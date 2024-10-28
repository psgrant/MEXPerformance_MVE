function [X,Y,Z] = EWPRotateAndTransform(Rotations,Translations,x,y,z)
phi   = Rotations(:,1);
theta = Rotations(:,2);
psi   = Rotations(:,3);

tx = Translations(:,1);
ty = Translations(:,2);
tz = Translations(:,3);

meanX = mean(minmax(x'));
meanY = mean(minmax(y'));
meanZ = mean(minmax(z'));


x = x - meanX;
y = y - meanY;
z = z - meanZ;

% Stupid dumb coord transfer
X = tx + x.*(cos(psi).*cos(theta) - sin(phi).*sin(psi).*sin(theta)) + y.*(cos(theta).*sin(psi) + cos(psi).*sin(phi).*sin(theta)) + z.*cos(phi).*sin(theta);
Y = ty - z.*sin(phi) + y.*cos(phi).*cos(psi) - x.*cos(phi).*sin(psi);
Z = tz - x.*(cos(psi).*sin(theta) + cos(theta).*sin(phi).*sin(psi)) - y.*(sin(psi).*sin(theta) - cos(psi).*cos(theta).*sin(phi)) + z.*cos(phi).*cos(theta);

% X = tx + x.*(cos(psi).*cos(theta) + sin(phi).*sin(psi).*sin(theta)) + y.*(cos(theta).*sin(psi) - cos(psi).*sin(phi).*sin(theta)) - z.*cos(phi).*sin(theta);
% Y = ty - z.*sin(phi) + y.*cos(phi).*cos(psi) - x.*cos(phi).*sin(psi);
% Z = tz + x.*(cos(psi).*sin(theta) - cos(theta).*sin(phi).*sin(psi)) + y.*(sin(psi).*sin(theta) + cos(psi).*cos(theta).*sin(phi)) + z.*cos(phi).*cos(theta);


X = X + meanX;
Y = Y + meanY;
Z = Z + meanZ;


X(abs(X) < sqrt(eps)) = 0;
Y(abs(Y) < sqrt(eps)) = 0;
Z(abs(Z) < sqrt(eps)) = 0;












end