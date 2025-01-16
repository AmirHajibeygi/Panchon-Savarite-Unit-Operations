tic;
%%% Inputs and Assumptions

% Antoine Coefficient for #Methanol and #Water
A_m=8.08097;        %methanol
B_m=1582.271;
C_m=239.73;
A_w=8.07131;        %water
B_w=1730.63;
C_w=233.43;

% Heat Capacity( %j/Kmole.K) and Latent Heat( %j/Kmole.K) for liquid and vapor #Methanol and #Water
cp_lm=81080;       %methanol
cp_vm=82656;      
L_m=33494000;
cp_lw=75240;       %water
cp_vw=41814;    
L_w=41652000; 

% The Equilibruim datas for #methanol and #water
X_methanol=[0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]; 
y_methanol=[0,0.261,0.405,0.571,0.66,0.72,0.772,0.82,0.867,0.912,0.955,1];

% Valus of excess enthalpy from Aspen+ for #methanol and #water mixture
x_m=[0.00107,0.0021,0.00412,0.00809,0.016,0.0641,0.1358,0.3172,0.3172,0.571,0.7516,0.8637,0.9276,0.9622,0.9805,0.99,0.9948,0.0318];
delta_h=[-7470,-14620,-28640,-55980,-110210,-403710,-698340,-886000,-880640,-748440,-544980,-351940,-203210,-112380,-60020,-31590,-16250,-213430];

% Inputs given for the problem
XD=0.99;
XB=0.01;
XF=0.5;
p1=760;     %mmHg(Total pressure)


%%% Calculations

% saturation T for #methanol and #water(Antoine Equation) and creating a
% matrix of tempuratures between them
T_sat1=B_m/(A_m-log10(760))-C_m;
T_sat2=B_w/(A_w-log10(760))-C_w;

a=T_sat1;
i=1;
while a<T_sat2 
    T(i)=a;
    a=a+2;
    i=i+1;
end
T(1,i)=T_sat2;

%Calculating saturated pressures(Antoine Equation)
x=zeros(1,19);
P_sat1=zeros(1,19);
P_sat2=zeros(1,19);
for i=1:19
    P_sat1(i)=10^(A_m-(B_m/(C_m+T(i))));
    P_sat2(i)=10^(A_w-(B_w/(C_w+T(i))));
    x(i)=(760-P_sat2(i))/(P_sat1(i)-P_sat2(i));
end

y1=zeros(1,19);
for i=1:19
    y1(i)=(x(i)*P_sat1(i))/p1;
end
y=double(y1);

% Writing the relation between excess enthalpy and mole fractions of
% methanol based on datas from Aspen+
p = polyfit(x_m,delta_h,5);
H_excess=@(x) p(1)*x^5+p(2)*x^4+p(3)*x^3+p(4)*x^2+p(5)*x+p(6);

H_excess1=zeros(1,19);
for i =1:19
    H_excess1(i)=H_excess(x(i));
end

% Writing the relation of liquid and vapor enthalpies based on mole
% fraction of #methanol
HL=zeros(1,19);
HV=zeros(1,19);
for i=1:19
    HL(i)=((cp_lm*x(i)+cp_lw*(1-x(i)))*(T(i)-25)+H_excess1(i));
    HV(i)=y(i)*(L_m+cp_vm*(T(i)-T_sat1))+(1-y(i))*(L_w+cp_vw*(T(i)-T_sat1));
end

b1=polyfit(x,HL,4);
H_L=@(x) b1(1)*x.^4+b1(2)*x.^3+b1(3)*x.^2+b1(4)*x+b1(5);

b2=polyfit(y,HV,3);
H_V=@(y) b2(1)*y.^3+b2(2)*y.^2+b2(3)*y+b2(4);

%Writing equilibrium equation (y and x for #methanol)
p = polyfit(X_methanol, y_methanol, 7);
y_equation= @(x) p(1)*x^7+p(2)*x^6+p(3)*x^5+p(4)*x^4+p(5)*x^3+p(6)*x^2+p(7)*x+p(8);

% Finding the location of important points (feed,min deltaD,deltaD,DeltaB)
y_feed=y_equation(XF);
H_vfeed=H_V(y_feed);
H_Lfeed=H_L(XF);
H_deltamin=((H_vfeed-H_Lfeed)/(y_feed-XF))*(XD-XF)+H_Lfeed;
Ratio=2*(H_deltamin-H_V(XD))/(H_V(XD)-H_L(XD));
H_deltaD=(H_V(XD)-H_L(XD))*Ratio+H_V(XD);
H_deltaB=((H_deltaD-H_Lfeed)/(XD-XF))*(XB-XF)+H_Lfeed;

%plot the H-x,y curves
plot(x,HL,y,HV)
grid on
hold on

%plot the tielines for the top section of the column
x_top1=XD;
y_top1=H_V(x_top1);
syms x
y=x_top1-y_equation(x)==0;
x_s = vpasolve(y, x );
x_top2=double(x_s(1));
i=0;
while x_top2>XF 
    y_top1=H_V(x_top1);
    syms x
    y=x_top1-y_equation(x)==0;
    x_s = vpasolve(y, x );
    x_top2=double(x_s(1));
    y_top2=H_L(x_top2);
    y_tie=@(x) ((H_deltaD-y_top2)/(XD-x_top2))*(x-XD)+H_deltaD;
    y2=H_V(x)-y_tie(x)==0;
    x_s1 = vpasolve(y2, x );
    x_top3=double(x_s1((1)));
    y_top3=H_V(x_top3);
    plot([x_top1 x_top2],[y_top1 y_top2],'r')
    plot([x_top2 XD],[y_top2 H_deltaD],'b')
    x_top1=x_top3;
    y_top1=y_top3;
    i=i+1;
end

%plot the tielines for the bottom section of the column
j=0;
x_bot1=x_top2;
while x_bot1>XB
    y_bot1=H_L(x_bot1);
    syms x
    y3= @(x) ((y_bot1-H_deltaB)/(x_bot1-XB))*(x-XB)+H_deltaB;
    y=H_V(x)-y3(x)==0;
    x_s = vpasolve(y, x );
    x_bot2=double(x_s(1));
    y_bot2=H_V(x_bot2);
    y4= x_bot2-y_equation(x)==0;
    x_s1 = vpasolve(y4, x );
    x_bot3=double(x_s1(1));
    y_bot3=H_L(x_bot3);
    plot([x_bot3 x_bot2],[y_bot3 y_bot2],'r')
    plot([x_bot2 XB],[y_bot2 H_deltaB],'b')
    x_bot1=x_bot3;
    y_bot1=y_bot3;
    j=j+1;
end
hold off


%%% Output

% The number of ideal Stages
disp(['The number of ideal stages = ',num2str(i+j)]);
toc;