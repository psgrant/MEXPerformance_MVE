function [SumXiXSatArea] = ComputeSurfaceFlux(Quads,QuadsIndex,BoundaryConds,QuadSCVArea,X,XSat)
SumXiXSatArea = 0;
for i = 1:length(Quads)

    Nodes = Quads(i,:);
    QuadFace = QuadsIndex(i);
    BCType = BoundaryConds(QuadFace);
    % WFlux = -qw * ((X(Nodes)-Properties.Xb(Nodes)) - XSat(Nodes));
    % Fix this hack later !!!!!!! add + FSP to XSat
    if BCType == 2 % Saturating
        SumXiXSatArea = SumXiXSatArea + sum((X(Nodes) - XSat(Nodes)) * 1 * QuadSCVArea(i));
    end

end