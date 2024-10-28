function [Properties] = ComputeMaterialProperties(X,T,Rho0,Constants)
%% Calculate material properties
Properties = struct;
% Constants

RhoS = Constants.RhoS; % kg / m^3
RhoW = Constants.RhoW; % kg / m^3
RhoFBar = Constants.RhoFBar*2; % kg / m^3
Patm = Constants.Patm;
CPw = Constants.CPw;
CPs = Constants.CPs;
CPa = Constants.CPa;
CPv = Constants.CPv;
h0vap = Constants.h0vap;
R = Constants.R;
Mv = Constants.Mv;
Ma = Constants.Ma;
Tr = Constants.Tr;
Ti = Constants.Ti;
DiffFact = Constants.DiffFact;
PermFact = Constants.PermFact;


Xfsp = 0.325 - 0.001 .* T;
XSat = (1 - Rho0/1530) * 1000 ./(Rho0) - 0*(Xfsp + 0);
% Important Properties

espX = 0.005;

% Smoothing of Xb
Xb = zeros(size(X));
for i = 1:length(X)

    Xp = X(i);
    X1 = Xfsp(i) - espX;
    X2 = Xfsp(i) + espX;


    if Xp < X1
        Xb(i) = Xp;
    elseif Xp < X2
        Xb(i) = -1/(4*espX) * (X1^2 - 2*X2*Xp + Xp^2);
    else
        Xb(i) = Xfsp(i);
    end
end

% epsS = Rho0 .* (RhoW + Rho0 .* Xb) ./ (RhoS .* (RhoW + Rho0 .* (Xb - Xfsp)));
epsS = (1 + RhoS.*Xb/RhoW) ./ (RhoS./Rho0 + (Xb - Xfsp)*RhoS./RhoW);
%
% Phi = 1 - Rho0 ./ RhoS;
Phi = 1 - epsS;
Sw = Rho0 .* (X - Xb) ./ (Phi*RhoW); % Saturation


%% Permeabilities
% absolute
a1 = [-0.82567033*10^-14,0.10275104*10^-10, 0.66577831*10^-8, 0.53733443*10^-6];
a2= [-0.14614421e-30,0.12115020e-27,0.46368212e-24,-0.71386999e-21,0.27514032e-18];

tv = polyval(a1,Rho0);
beta = polyval(a2,Rho0);


at = 50e-6;
al = 1.5e-3;
ar = 0.575e-4 - 0.375e-7*Rho0;

alpha = 1.8e24*(tv+2e-6);
gamma = 25 / at;
nt = gamma * (ar - 2*tv);
nr = gamma * (at - 2*tv);
nl = nt + nr;

% Final form
KR = PermFact * nr .* ar ./ (2*at.*al.*alpha);
KT = PermFact * nt .* at ./ (2*ar.*al.*alpha);
KL = PermFact * al ./ (at.*ar.*(alpha./nl + al./beta));

% Relative
kwR = Sw.^3;
kwT = Sw.^3;
kwL = Sw.^8;

kgT = 1 + (2*Sw - 3) .*Sw.^2;

%% Effective diffusivity

Dv = 2.26e-5*((T+Tr)/Tr).^1.81;

DeffT = 0.001 * kgT .* Dv;
DeffR = 2 * DeffT;
DeffL = 10 * DeffT;

%% Vapour pressure

Psat = exp(25.5058 - 5204.9 ./ (T + Tr));
BoundWaterFrac = (Xb./Xfsp);
Psi = 1 - exp(-0.76427*BoundWaterFrac - 3.36787*BoundWaterFrac.^2);

Pv = Psat .* Psi;


%% Capillary Pressure (water)

Pc = (0.0775 - 0.000185 * T) .* (3150 ./ ((Sw + 0.0001)) ...
    - (1047 + 3.368*Rho0)./((1.02 - Sw)) +...
    149.8 * Rho0 .* (1 - Sw) + 52350 + 168.4*Rho0 - 3150 / (1 +0.0001));
% Rho0 is at the CV
Pw = Patm - Pc;

%% Dynamic Viscosoity

muw = RhoW * exp(-19.143 + 1540./(T+Tr));

%% Air density
RhoABar = Ma*(Patm - Pv)./(R*(T+Tr));
%% Volume fractions

a = Ma.*Pv ./ (R*(T+Tr));
b = -(a .* Phi .* (1-Sw)+(RhoFBar+RhoABar));
c = Phi .* (1-Sw) * RhoFBar;

% final form
epsF = (-b - sqrt(b.^2 - 4*a.*c)) ./ (2*a);
epsW = Phi .* Sw;
% epsG = Phi - epsW - epsF;
epsG = Phi.*(1-Sw) - epsF;

% Evaluating DrayT at 500 kg m^3

epsSt = 500 .* (RhoW + 500 .* Xb) ./ (RhoS .* (RhoW + 500 .* (Xb - Xfsp)));
% epsSt = (1 + RhoS.*Xb/RhoW) ./ (RhoS./500 + (Xb - Xfsp)*RhoS./RhoW);
phit = 1 - epsSt;
% epsSt = 1-phit;
Xsatt = phit .* RhoW/500;
Swt = phit .* (X-Xb)./Xsatt; %is 0

a = Ma*Pv ./ (R*(T+Tr));
b = -(a .* phit .* (1-Swt)+(RhoFBar+RhoABar));
c = phit .* (1-Swt) .* RhoFBar;


epsFt = (-b - sqrt(b.^2 - 4*a.*c)) ./ (2*a);
epsGt = phit.*(1 - Sw) - epsFt;



%% Phase Densities
% need values for X, T and RhoABar

RhoV = Pv .* Mv ./ (R * (T+Tr));
RhoA = ((Patm-Pv)*Ma) ./(R * (T+Tr));
RhoG = RhoA + RhoV;


%% Mass fractions
% Air
wa = RhoA ./ RhoG;
% Vapour
wv = RhoV ./ RhoG;


%% Phase enthalpies

hw = CPw * T;
ha = CPa * T;
hs = CPs * T;
hv = h0vap + CPv*T;
hb = hw - 0.4 * (hv - hw) .* (1 - Xb./Xfsp).^2;

%% Bound Diffusivity

% Partial derivative dPv/dXb
DiffPv = Psat .* (0.76427./Xfsp + 2*3.6787 .* BoundWaterFrac) .*...
    (exp(-0.76427*BoundWaterFrac - 3.36787*BoundWaterFrac.^2));

% Compute pseudo micro diffusivities
DbMicro = RhoS .* exp(-12.8183993 + 10.8951601 * Xb - 4300./(T+Tr));
DvMicro = Mv ./ (R*(T + Tr)) * 2.26e-5 .* DiffPv;

Aseries = 1 ./ (epsS./DbMicro + epsG./DvMicro);

DbRMacro = Aseries .* (1 + 1.6*epsG.^1.8);
DbTMacro = 1.8*Aseries;

epsRay = 0.05;
DRayR = 0.1 * DvMicro;

DRayT = 1 ./ (epsSt ./ DbMicro + epsGt ./ DvMicro);

% Bound liquid pseudo-diffusivity
DbDiffR = (1 - epsRay)*DbRMacro + epsRay*DRayR;
DbDiffT =  1 ./ ((1-epsRay)./DbTMacro + epsRay ./ DRayT);
DbDiffL = epsS .* DbMicro + epsG .* DvMicro;

% Effective diffusivity
DbR = DiffFact*DbDiffR ./ Rho0;
DbT = DiffFact*DbDiffT ./ Rho0;
DbL = DiffFact*DbDiffL ./ Rho0;

%% Conduction time
LambdaAir = 0.023; % W/m/K
LambdaW = 0.6;
LambdaSPerp = 0.5; % W/m/K
n = 0.6;

LambdaSeries = 1./ (epsS./LambdaSPerp + epsG./LambdaAir);%
LambdaParallel = epsS.*LambdaSPerp + epsG*LambdaAir;%


% macroscopic conductivies
LambdaRMacro = 0.46*LambdaSeries + 0.54*LambdaParallel;%
LambdaTMacro = (epsG*(LambdaAir).^n + epsS*(LambdaSPerp).^n).^(1/n);%


% ray corrections
LambdaRayPerp =  (epsGt*LambdaAir.^n + epsSt*LambdaSPerp.^n).^(1/n);%
LambdaRayParallel = epsGt*LambdaAir + 2*epsSt*LambdaSPerp;%
%     LambdaRayParallel =(epsSt*LambdaAir.^n + epsGt.*LambdaParallel.^n).^(1/n);

% Effective conductivities
Lambda0R = (1-epsRay).*LambdaRMacro + epsRay*LambdaRayParallel;%
Lambda0T = 1 ./ ((1-epsRay)./LambdaTMacro + epsRay./LambdaRayPerp);%
Lambda0L = epsG*LambdaAir + epsS.*(2*LambdaSPerp);

KeffR = ((1-Sw).*(Lambda0R).^n + Sw.*(Phi*LambdaW^n + epsS*LambdaSPerp^n)).^(1/n);
KeffT = ((1-Sw).*(Lambda0T).^n + Sw.*(Phi*LambdaW^n + epsS*LambdaSPerp^n)).^(1/n);
KeffL = ((1-Sw).*(Lambda0L).^n + Sw.*(Phi*LambdaW + epsS*2*LambdaSPerp).^n).^(1/n);
% NumNodes = length(X);
% KeffR =  0.5.* ones(NumNodes,1);
% KeffT =  0.5.* ones(NumNodes,1);
% KeffL =  0.5.* ones(NumNodes,1);

%% Put into property Struct

Properties.Xb = Xb;
Properties.Xfsp = Xfsp;

Properties.Pv = Pv;
Properties.Pw = Pw;
Properties.muw = muw;

Properties.epsS = epsS;
Properties.epsW = epsW;
Properties.epsG = epsG;
Properties.epsF = epsF;

Properties.RhoV = RhoV;
Properties.RhoA = RhoA;
Properties.RhoG = RhoG;
Properties.RhoABar = RhoABar;

Properties.wv = wv;
Properties.wa = wa;

Properties.hw = hw;
Properties.ha = ha;
Properties.hb = hb;
Properties.hv = hv;
Properties.hs = hs;

Properties.DbR = DbR;
Properties.DbT = DbT;
Properties.DbL = DbL;

Properties.KeffR = KeffR;
Properties.KeffT = KeffT;
Properties.KeffL = KeffL;

Properties.DeffT = DeffT;
Properties.DeffR = DeffR;
Properties.DeffL = DeffL;

Properties.kwR = kwR;
Properties.kwT = kwT;
Properties.kwL = kwL;

Properties.KR = KR;
Properties.KT = KT;
Properties.KL = KL;

Properties.Sw = Sw;
Properties.XSat = XSat;
end