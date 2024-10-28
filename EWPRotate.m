function [X,Y,Z] = EWPRotate(Rotations,x,y,z)


phi   = Rotations(:,1);
theta = Rotations(:,2);
psi   = Rotations(:,3);

% Stupid dumb coord transfer for middle board
% X = x.*(cos(psi).*cos(theta) + sin(phi).*sin(psi).*sin(theta)) + y.*(cos(theta).*sin(psi) - cos(psi).*sin(phi).*sin(theta)) - z.*cos(phi).*sin(theta);
% Y = - z.*sin(phi) + y.*cos(phi).*cos(psi) - x.*cos(phi).*sin(psi);
% Z = x.*(cos(psi).*sin(theta) - cos(theta).*sin(phi).*sin(psi)) + y.*(sin(psi).*sin(theta) + cos(psi).*cos(theta).*sin(phi)) + z.*cos(phi).*cos(theta);


X =  x.*(cos(psi).*cos(theta) - sin(phi).*sin(psi).*sin(theta)) + y.*(cos(theta).*sin(psi) + cos(psi).*sin(phi).*sin(theta)) + z.*cos(phi).*sin(theta);
Y = -z.*sin(phi) + y.*cos(phi).*cos(psi) - x.*cos(phi).*sin(psi);
Z = -x.*(cos(psi).*sin(theta) + cos(theta).*sin(phi).*sin(psi)) - y.*(sin(psi).*sin(theta) - cos(psi).*cos(theta).*sin(phi)) + z.*cos(phi).*cos(theta);


X(abs(X) < eps) = 0;
Y(abs(Y) < eps) = 0;
Z(abs(Z) < eps) = 0;
end