function [GradPhiX, GradPhiY, GradPhiZ] = TPEShapeFuncGradVectorised(A11, A12, ...
    A13, A21, A22, A23, A31, A32, A33, Phi, Elements)
% Computes the gradients of the shape functions

% Set up local corods array

LocalCoordsArray = [5/12 1/6 -0.5;
    5/12 5/12 -0.5;
    1/6 5/12 -0.5;
    5/12 1/6 0.5;
    5/12 5/12 0.5;
    1/6 5/12 0.5;
    5/24 5/24 0;
    7/12 5/24 0;
    5/24 7/12 0];

% Make element array where phi is at every index in element
PhiBigElement = Phi(Elements);

NineZeros = zeros(9,1);

% Create shape derivatives
dxi = 0.5*[-(1-LocalCoordsArray(:,3)), (1-LocalCoordsArray(:,3)), NineZeros, -(1+LocalCoordsArray(:,3)), (1+LocalCoordsArray(:,3)), NineZeros]';
deta = 0.5*[-(1-LocalCoordsArray(:,3)), NineZeros, (1-LocalCoordsArray(:,3)), -(1+LocalCoordsArray(:,3)), NineZeros, (1+LocalCoordsArray(:,3))]';
dzeta = 0.5* [-(1 - LocalCoordsArray(:,1) - LocalCoordsArray(:,2)), -LocalCoordsArray(:,1), -LocalCoordsArray(:,2),(1 - LocalCoordsArray(:,1) - LocalCoordsArray(:,2)), LocalCoordsArray(:,1), LocalCoordsArray(:,2)]';


% Create big derivative vector
dPhidxi = PhiBigElement * dxi;
dPhideta = PhiBigElement * deta;
dPhidzeta = PhiBigElement * dzeta;
% h = (eps);
% dPhidxi(abs(dPhidxi)<h) = 0;
% dPhideta(abs(dPhideta)<h) = 0;
% dPhidzeta(abs(dPhidzeta)<h) = 0;

% Big matrix multiply with JinV
GradPhiX = A11.* dPhidxi + A12.* dPhideta + A13.* dPhidzeta;
GradPhiY = A21.* dPhidxi + A22.* dPhideta + A23.* dPhidzeta;
GradPhiZ = A31.* dPhidxi + A32.* dPhideta + A33.* dPhidzeta;
end