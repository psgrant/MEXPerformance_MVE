function [XX,YY,ZZ,XY,XZ,YZ] = ActionRotationVectors(RotationVectors,R,T,L)
% Applys the rotation vectors of the RTL (diagonal) tensors, this outputs a
% full symetric tensor elements.
% Pull out of array
XX_R = RotationVectors(:,1);
XX_T = RotationVectors(:,2);
XX_L = RotationVectors(:,3);

YY_R = RotationVectors(:,4);
YY_T = RotationVectors(:,5);
YY_L = RotationVectors(:,6);

ZZ_R = RotationVectors(:,7);
ZZ_T = RotationVectors(:,8);
ZZ_L = RotationVectors(:,9);

XY_R = RotationVectors(:,10);
XY_T = RotationVectors(:,11);
XY_L = RotationVectors(:,12);

XZ_R = RotationVectors(:,13);
XZ_T = RotationVectors(:,14);
XZ_L = RotationVectors(:,15);

YZ_R = RotationVectors(:,16);
YZ_T = RotationVectors(:,17);
YZ_L = RotationVectors(:,18);


% Compute rotated material props
XX = XX_R.*R + XX_T.*T + XX_L.*L;
YY = YY_R.*R + YY_T.*T + YY_L.*L;
ZZ = ZZ_R.*R + ZZ_T.*T + ZZ_L.*L;

XY =  XY_R.*R - XY_T.*T - XY_L.*L;
XZ = -XZ_R.*R - XZ_T.*T + XZ_L.*L;
YZ = -YZ_R.*R + YZ_T.*T - YZ_L.*L;


end

