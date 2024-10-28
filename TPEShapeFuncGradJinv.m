function [A11, A12, A13, A21, A22, A23, A31, A32, A33] = TPEShapeFuncGradJinv(NodePos, Elements)
%% This code generates all the values for the inverse jacobian for compute the shape function gradients
GradPhiX = zeros(size(Elements,1),9);
GradPhiY = GradPhiX;
GradPhiZ = GradPhiX;
LocalCoordsArray = [5/12 1/6 -0.5;
    5/12 5/12 -0.5;
    1/6 5/12 -0.5;
    5/12 1/6 0.5;
    5/12 5/12 0.5;
    1/6 5/12 0.5;
    5/24 5/24 0;
    7/12 5/24 0;
    5/24 7/12 0];


A11 = zeros(size(Elements,1),9);
A12 = zeros(size(Elements,1),9);
A13 = zeros(size(Elements,1),9);

A21 = zeros(size(Elements,1),9);
A22 = zeros(size(Elements,1),9);
A23 = zeros(size(Elements,1),9);

A31 = zeros(size(Elements,1),9);
A32 = zeros(size(Elements,1),9);
A33 = zeros(size(Elements,1),9);


NineZeros = zeros(9,1);

% 6 nodes, 9 IPs
dxi = 0.5*[-(1-LocalCoordsArray(:,3)), (1-LocalCoordsArray(:,3)), NineZeros, -(1+LocalCoordsArray(:,3)), (1+LocalCoordsArray(:,3)), NineZeros]';
deta = 0.5*[-(1-LocalCoordsArray(:,3)), NineZeros, (1-LocalCoordsArray(:,3)), -(1+LocalCoordsArray(:,3)), NineZeros, (1+LocalCoordsArray(:,3))]';
dzeta = 0.5* [-(1 - LocalCoordsArray(:,1) - LocalCoordsArray(:,2)), -LocalCoordsArray(:,1), -LocalCoordsArray(:,2),(1 - LocalCoordsArray(:,1) - LocalCoordsArray(:,2)), LocalCoordsArray(:,1), LocalCoordsArray(:,2)]';

for ele = 1:size(Elements,1)

    ElementIds = Elements(ele,:);
    GlobalCoords = NodePos(ElementIds,:);
    X = GlobalCoords(:,1)';
    Y = GlobalCoords(:,2)';
    Z = GlobalCoords(:,3)';


    dXdx = X * dxi;
    dYdx = Y * dxi;
    dZdx = Z * dxi;

    dXdn = X * deta;
    dYdn = Y * deta;
    dZdn = Z * deta;

    dXdz = X * dzeta;
    dYdz = Y * dzeta;
    dZdz = Z * dzeta;


    for i = 1:9
        J = [dXdx(i), dYdx(i), dZdx(i);
            dXdn(i), dYdn(i), dZdn(i);
            dXdz(i), dYdz(i), dZdz(i)];

        Jinv = inv(J);



        % only depend on spatial stuff, AXX is E x 9 in size gross but only
        % need to precompute once lol
        A11(ele,i) = Jinv(1,1);
        A12(ele,i) = Jinv(1,2);
        A13(ele,i) = Jinv(1,3);

        A21(ele,i) = Jinv(2,1);
        A22(ele,i) = Jinv(2,2);
        A23(ele,i) = Jinv(2,3);

        A31(ele,i) = Jinv(3,1);
        A32(ele,i) = Jinv(3,2);
        A33(ele,i) = Jinv(3,3);
        %%% Sparse matrix block for Jacobian backslash
        %%% for stability


        % A11(ele,:) = -(dYdn.*dZdz - dYdz.*dZdn)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        % A12(ele,:) = (dYdx.*dZdz - dYdz.*dZdx)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        % A13(ele,:) = (dYdn.*dZdx - dYdx.*dZdn)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        %
        % A21(ele,:) = (dXdn.*dZdz - dXdz.*dZdn)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        % A22(ele,:) = -(dXdx.*dZdz - dXdz.*dZdx)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        % A23(ele,:) = -(dXdn.*dZdx - dXdx.*dZdn)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        %
        % A31(ele,:) = -(dXdn.*dYdz - dXdz.*dYdn)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        % A32(ele,:) = (dXdx.*dYdz - dXdz.*dYdx)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
        % A33(ele,:) = (dXdn.*dYdx - dXdx.*dYdn)./(dXdn.*dYdx.*dZdz - dXdn.*dYdz.*dZdx - dXdx.*dYdn.*dZdz + dXdx.*dYdz.*dZdn + dXdz.*dYdn.*dZdx - dXdz.*dYdx.*dZdn);
    end
end
h = (eps);
% A11(abs(A11) < h) = 0;
% A12(abs(A12) < h) = 0;
% A13(abs(A13) < h) = 0;
% 
% A21(abs(A21) < h) = 0;
% A22(abs(A22) < h) = 0;
% A23(abs(A23) < h) = 0;
% 
% A31(abs(A31) < h) = 0;
% A32(abs(A32) < h) = 0;
% A33(abs(A33) < h) = 0;













