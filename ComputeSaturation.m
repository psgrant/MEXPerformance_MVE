function [Sw] = ComputeSaturation(X,T,Rho0,Constants)
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


Xfsp = 0.325 - 0.001 .* T;
XSat = (1 - Rho0/1530) * 1000 ./(Rho0) - (0.001);
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
