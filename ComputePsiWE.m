function [PsiW,PsiE] = ComputePsiWE(X,T,Rho0,Constants)


Tbp = 100;
h0vap = Constants.h0vap;
[Properties] = ComputeMaterialProperties(X,T,Rho0,Constants);
% Pass in temp not theta

%% Read in properties
Xb = Properties.Xb;
Xfsp = Properties.Xfsp;
epsG = Properties.epsG;
RhoV = Properties.RhoV;
RhoA = Properties.RhoA;


hw = Properties.hw;
ha = Properties.ha;
hv = Properties.hv;
hs = Properties.hs;

% Enthalpy Terms
LiqSolEnthalpy = Rho0 .* (X .* hw + hs);
% LiqSolEnthalpy = Rho0 .* hs; % Rh0 * Cps * Tbp*Theta

GasEnthalpy = epsG .* (RhoV .* hv + RhoA .* ha);
LatentHeat = 0.4 * Rho0 .* Xb .* (hv - hw) .* (1 - Xb./Xfsp + 1/3*(Xb./Xfsp).^2);

% % Acumulation terms
PsiW = (Rho0 .* X + epsG .* RhoV);
PsiE = LiqSolEnthalpy + GasEnthalpy - LatentHeat;

% Acumulation terms
% PsiW = (Rho0 .* X) ./ Rho0;
% PsiE = (Rho0 .* Constants.CPs .* T) ./ (1400 * Rho0);
% PsiE = (Rho0 .* T) ./ (Rho0); / (Tbp * h0vap)

end