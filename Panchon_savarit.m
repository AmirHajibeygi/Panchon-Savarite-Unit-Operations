tic;
%%% Generic Inputs (User-defined)
disp('Enter Antoine Coefficients for Component 1 (A, B, C):');
A1 = input('A1: '); B1 = input('B1: '); C1 = input('C1: ');
disp('Enter Antoine Coefficients for Component 2 (A, B, C):');
A2 = input('A2: '); B2 = input('B2: '); C2 = input('C2: ');

disp('Enter Heat Capacity (J/kmol.K) and Latent Heat (J/kmol) for both components:');
cp_l1 = input('Liquid Cp of Component 1: ');
cp_v1 = input('Vapor Cp of Component 1: ');
L1 = input('Latent Heat of Component 1: ');
cp_l2 = input('Liquid Cp of Component 2: ');
cp_v2 = input('Vapor Cp of Component 2: ');
L2 = input('Latent Heat of Component 2: ');

disp('Enter Equilibrium data for the components:');
X_comp1 = input('Mole fractions of Component 1 in liquid phase: ');
Y_comp1 = input('Mole fractions of Component 1 in vapor phase: ');

disp('Enter Excess Enthalpy data (x and delta_h values):');
x_excess = input('Mole fractions for Excess Enthalpy data: ');
delta_h = input('Excess Enthalpy data: ');

disp('Enter the feed composition and system pressure:');
XD = input('Distillate Composition: ');
XB = input('Bottom Composition: ');
XF = input('Feed Composition: ');
P = input('Total Pressure (mmHg): ');

%%% Calculations
T_sat1 = B1 / (A1 - log10(P)) - C1;
T_sat2 = B2 / (A2 - log10(P)) - C2;

% Generate temperature range between saturation temperatures
T = linspace(T_sat1, T_sat2, 100);

% Calculate saturation pressures
P_sat1 = 10.^(A1 - B1 ./ (C1 + T));
P_sat2 = 10.^(A2 - B2 ./ (C2 + T));
x = (P - P_sat2) ./ (P_sat1 - P_sat2);
y = (x .* P_sat1) / P;

% Fit Excess Enthalpy data
p = polyfit(x_excess, delta_h, 5);
H_excess = @(x) polyval(p, x);

% Calculate liquid and vapor enthalpies
HL = ((cp_l1 .* x + cp_l2 .* (1 - x)) .* (T - 25)) + arrayfun(H_excess, x);
HV = y .* (L1 + cp_v1 .* (T - T_sat1)) + (1 - y) .* (L2 + cp_v2 .* (T - T_sat2));

% Fit enthalpy data
H_L = @(x) polyval(polyfit(x, HL, 4), x);
H_V = @(y) polyval(polyfit(y, HV, 3), y);

% Fit equilibrium curve
y_equation = @(x) polyval(polyfit(X_comp1, Y_comp1, 7), x);

% Feed condition enthalpies
y_feed = y_equation(XF);
H_vfeed = H_V(y_feed);
H_Lfeed = H_L(XF);
H_deltamin = ((H_vfeed - H_Lfeed) / (y_feed - XF)) * (XD - XF) + H_Lfeed;
Ratio = 2 * (H_deltamin - H_V(XD)) / (H_V(XD) - H_L(XD));
H_deltaD = (H_V(XD) - H_L(XD)) * Ratio + H_V(XD);
H_deltaB = ((H_deltaD - H_Lfeed) / (XD - XF)) * (XB - XF) + H_Lfeed;

% Plot H-x,y curves
plot(x, HL, y, HV);
grid on;
hold on;

% Plot tielines (Top section)
x_top1 = XD;
while x_top1 > XF
    y_top1 = H_V(x_top1);
    x_top2 = fzero(@(x) x_top1 - y_equation(x), [0, x_top1]);
    y_top2 = H_L(x_top2);
    y_tie = @(x) ((H_deltaD - y_top2) / (XD - x_top2)) * (x - XD) + H_deltaD;
    x_top3 = fzero(@(x) H_V(x) - y_tie(x), [0, x_top1]);
    plot([x_top1 x_top2], [y_top1 y_top2], 'r');
    plot([x_top2 XD], [y_top2 H_deltaD], 'b');
    x_top1 = x_top3;
end

% Plot tielines (Bottom section)
x_bot1 = x_top2;
while x_bot1 > XB
    y_bot1 = H_L(x_bot1);
    y_tie = @(x) ((y_bot1 - H_deltaB) / (x_bot1 - XB)) * (x - XB) + H_deltaB;
    x_bot2 = fzero(@(x) H_V(x) - y_tie(x), [0, x_bot1]);
    x_bot3 = fzero(@(x) x_bot2 - y_equation(x), [0, x_bot2]);
    y_bot3 = H_L(x_bot3);
    plot([x_bot3 x_bot2], [y_bot3 y_bot1], 'r');
    plot([x_bot2 XB], [y_bot1 H_deltaB], 'b');
    x_bot1 = x_bot3;
end

hold off;

%%% Output
disp(['The number of ideal stages = ', num2str(i + j)]);
toc;