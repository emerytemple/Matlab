P = 101.325; % kPa
T = 293.15; % K
M = 0.2;
g = 9.81; % m/s

% gas side properties
[a, Cp, rho, h, gamma, k, mu, Pr] = refpropm('ACDHKLV^','T',293.15,'P',101.325,'air.ppf');
% a, m/s, Cp, J/kg.K, rho, kg/m^3, h, J/kg, gamma, -, k, W/m.K, mu, Pa.s, Pr, -

Po = P*Supersonic('Po/P','M',M,'g',gamma);
To = T*Supersonic('To/T','M',M,'g',gamma);

r = Pr^0.33; % 0.5 if lam, 0.33 if turb

top = 1+(r*((gamma-1)/2)*(M^2));
bot = 1+(((gamma-1)/2)*(M^2));
Taw = To*(top/bot);

v = M*a;
hg_est = (rho*v)^0.8;

eps = Supersonic('A/A*','M',M,'g',gamma);
Dt
cstar
R

top = (0.026/(Dt^0.2))*(((mu^0.2)*Cp)/(Pr^0.6))*(((Po*g)/cstar)^0.8)*((Dt/R)^0.1)*((1/eps)^0.9)
bot = 
hg = 


% coolant side properties
[a, Cp, rho, h, gamma, k, mu, Pr] = refpropm('ACDHKLV^','T',293.15,'P',101.325,'water.fld');

