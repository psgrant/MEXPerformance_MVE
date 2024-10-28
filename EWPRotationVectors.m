function [RotationVectors] = EWPRotationVectors(Rotations,GA)
% precomputation of the rotation vectors so we're not computing them at
% every iteration


% Load Rotations
phi   = Rotations(:,1);
theta = Rotations(:,2);
psi   = Rotations(:,3);

RotationVectors = zeros(length(Rotations),18);

% Compute the common 'subfunctions'
R1 = (cos(GA).*(cos(psi).*cos(theta) + sin(phi).*sin(psi).*sin(theta)) - sin(GA).*(cos(theta).*sin(psi) - cos(psi).*sin(phi).*sin(theta)));

R2 = (cos(GA).*(cos(psi).*sin(theta) - cos(theta).*sin(phi).*sin(psi)) - sin(GA).*(sin(psi).*sin(theta) + cos(psi).*cos(theta).*sin(phi)));

R3 = (cos(GA).*cos(phi).*sin(psi) + sin(GA).*cos(phi).*cos(psi));

R4 = (cos(psi).*(cos(GA).*cos(theta) + sin(GA).*sin(phi).*sin(theta)) - sin(psi).*(sin(GA).*cos(theta) - cos(GA).*sin(phi).*sin(theta)));

R5 = (cos(psi).*(cos(GA).*sin(theta) - sin(GA).*cos(theta).*sin(phi)) - sin(psi).*(sin(GA).*sin(theta) + cos(GA).*cos(theta).*sin(phi)));

T1 = (cos(GA).*(cos(theta).*sin(psi) - cos(psi).*sin(phi).*sin(theta)) + sin(GA).*(cos(psi).*cos(theta) + sin(phi).*sin(psi).*sin(theta)));

T2 = (cos(GA).*(sin(psi).*sin(theta) + cos(psi).*cos(theta).*sin(phi)) + sin(GA).*(cos(psi).*sin(theta) - cos(theta).*sin(phi).*sin(psi)));

T3 = (cos(GA).*cos(phi).*cos(psi) - sin(GA).*cos(phi).*sin(psi));

T4 = (cos(psi).*(sin(GA).*cos(theta) - cos(GA).*sin(phi).*sin(theta)) + sin(psi).*(cos(GA).*cos(theta) + sin(GA).*sin(phi).*sin(theta)));

T5 = (cos(psi).*(sin(GA).*sin(theta) + cos(GA).*cos(theta).*sin(phi)) + sin(psi).*(cos(GA).*sin(theta) - sin(GA).*cos(theta).*sin(phi)));

L1 = cos(phi).^2 .* sin(theta).^2;

L2 = sin(phi); 

L3 = cos(phi).^2 .* cos(theta).^2;

L4 = cos(phi) .* sin(phi) .* sin(theta);

L5 = cos(phi).^2 .* cos(theta) .* sin(theta);

L6 = cos(phi) .* cos(theta) .* sin(phi);

% Rotation vectors (doesn't change for the duration of the simulation, do
% as precomputation step)
XX_R = R1.^2; XX_T = T1.^2; XX_L = L1;
YY_R = R3.^2; YY_T = T3.^2; YY_L = L2.^2;
ZZ_R = R2.^2; ZZ_T = T2.^2; ZZ_L = L3;

XY_R = R4 .* R3; XY_T = T4 .* T3; XY_L = L4;
XZ_R = R1 .* R2; XZ_T = T1 .* T2; XZ_L = L5;
YZ_R = R3 .* R5; YZ_T = T3 .* T5; YZ_L = L6;

% Stick into an array
RotationVectors(:,1) = XX_R;
RotationVectors(:,2) = XX_T;
RotationVectors(:,3) = XX_L;

RotationVectors(:,4) = YY_R;
RotationVectors(:,5) = YY_T;
RotationVectors(:,6) = YY_L;

RotationVectors(:,7) = ZZ_R;
RotationVectors(:,8) = ZZ_T;
RotationVectors(:,9) = ZZ_L;

RotationVectors(:,10) = XY_R;
RotationVectors(:,11) = XY_T;
RotationVectors(:,12) = XY_L;

RotationVectors(:,13) = XZ_R;
RotationVectors(:,14) = XZ_T;
RotationVectors(:,15) = XZ_L;

RotationVectors(:,16) = YZ_R;
RotationVectors(:,17) = YZ_T;
RotationVectors(:,18) = YZ_L;

end