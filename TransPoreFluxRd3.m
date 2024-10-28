function [g,FluxDiff] = TransPoreFluxRd3(EWPMesh,Constants,XTheta)
% Extract stuff from struct
NodePos = EWPMesh.NodePos;
Elements = EWPMesh.Elements;
Triangles = EWPMesh.Triangles;
Quads = EWPMesh.Quads;
TrianglesIndex = EWPMesh.TrianglesSurfIndex;
QuadsIndex = EWPMesh.QuadsSurfIndex;
% TriangleNormals = EWPMesh.SurfTriNormalVectors;
% QuadNormals = EWPMesh.SurfTriNormalVectors;
% ScalingFacts = EWPMesh.ScaleFacts;
NormalVectorsX = EWPMesh.NormalVecsX;
NormalVectorsY = EWPMesh.NormalVecsY;
NormalVectorsZ = EWPMesh.NormalVecsZ;
% Alpha2 = EWPMesh.Alpha2;
% Beta2 = EWPMesh.Beta2;
% AlphaBeta = EWPMesh.AlphaBeta;
% CenterPos = EWPMesh.CenterPos;
% SubVolumes = EWPMesh.SubVolumes;
SubAreas = EWPMesh.SubAreas;
InterlacedFullVolumes = EWPMesh.InterlacedFullVolumes;
FullVolumes = EWPMesh.FullVolumes;
TriSCVArea = EWPMesh.SurfTriArea;
QuadSCVArea = EWPMesh.SurfQuadArea;
NodeRho0 = EWPMesh.NodeRho0;
ElementRho0 = EWPMesh.ElementRho0;
NodalGrainAngle = EWPMesh.NodalGrainAngle;
ShapeFuncArray = EWPMesh.ShapeFuncArray;
ShapeFuncArrayT = EWPMesh.ShapeFuncArray';
XSat = EWPMesh.XSat;
GlueInteractions = EWPMesh.GlueInteractions;
NodeRotations = EWPMesh.NodeRotations;
UseMex = EWPMesh.UseMex;
% BC stuff
SaturatingParam = EWPMesh.sigma_w;
HglueX = EWPMesh.HglueX;
HglueT = EWPMesh.HglueT;
pTdb = EWPMesh.pTdb;
pRelHum =  EWPMesh.pRelHum;
% B = EWPMesh.b;


BoundaryConds = EWPMesh.BoundaryConds;
% ElementRotations = EWPMesh.ElementRotations;

% Extract inv jacobian terms for grad shape funcs
A11 = EWPMesh.A11;
A12 = EWPMesh.A12;
A13 = EWPMesh.A13;
A21 = EWPMesh.A21;
A22 = EWPMesh.A22;
A23 = EWPMesh.A23;
A31 = EWPMesh.A31;
A32 = EWPMesh.A32;
A33 = EWPMesh.A33;


% NumMidNodes = EWPMesh.NumMidNodes;
NumMidNodes = 0;

% get the number of nodes and elements
NumNodes = length(NodeRho0);
NumNodes2 = NumNodes * 2;
NumElements = length(ElementRho0);

X = XTheta(1:2:NumNodes2);
Theta = XTheta(2:2:NumNodes2);

% Specify constants
% Constants.RhoS = 1530;           % kg/m^3    - Solid phase density
RhoW = Constants.RhoW;             % kg/m^3    - Density of water
% Constants.RhoFBar = 0.002;       % kg/mol    - Superficial density (fixed phase)
Patm = Constants.Patm;             % Pa        - Atmospheric pressure
% Constants.CPw = 4187;            % J/(kg K)  - Specific heat of water
% Constants.CPs = 1400;            % J/(kg K)  - Specific heat of solid
% Constants.CPa = 1005.683;        % J/(kg K)  - Specific heat of air
% Constants.CPv = 1900;            % J/(kg K)  - Specific heat of vapor
% Constants.h0vap = 2.503e6;       % J/kg      - Latent heat of vaporization
R = Constants.R;                   % J/(mol K) - Universal gas constant
Mv = Constants.Mv;                 % kg/mol    - Molecular weight of vapor
% Constants.Ma = 0.028952;         % kg/mol    - Molecular weight of air
Tr = Constants.Tr;                 % K         - Reference temperature
q = Constants.q;                   %           - Mass transfer coeff
km = Constants.km;                 %           - heat transfer coeff
Constants.Ti = 25;                 % C         - Initial temperature
Tbp = 100;                         % C         - Boiling point of water
Gravity = 9.81;
T = Theta * Tbp;
RhoWGrav = (RhoW*Gravity);
Time = Constants.Time;
MeshVolume = Constants.MeshVolume;
% Compute the material properties
[Properties] = ComputeMaterialProperties(X,T,NodeRho0,Constants);



Xb = Properties.Xb; % Bound Water
% Xfsp = Properties.Xfsp; % Fibre Saturation Point
Pv = Properties.Pv; % Vapour pressure
Pw = Properties.Pw; % Water pressure
muw = Properties.muw; % Dynamic Viscosity
% epsS = Properties.epsS; % Solid mass fraction
% epsW = Properties.epsW; % Water mass fraction
% epsG = Properties.epsG; % Gas mass Fraction
% epsF = Properties.epsF; % Fixed gas phase fraction
% RhoV = Properties.RhoV; % Vapour density
% RhoA = Properties.RhoA; % Air phase Density
RhoG = Properties.RhoG; % Gas phase Density
% RhoABar = Properties.RhoABar; % Air Density
wv = Properties.wv; % Vapour mass fraction
% wa = Properties.wa; % Air mass fraction
hw = Properties.hw; % Water Enthalpy-Temperature Relation
ha = Properties.ha; % Air Enthalpy-Temperature Relation
hb = Properties.hb; % Bound Water Enthalpy-Temperature Relation
hv = Properties.hv; % Vapour Enthalpy-Temperature Relation
% hs = Properties.hs; % Solid Enthalpy-Temperature Relation
XSat = Properties.XSat;
muw = muw(Elements)';
RhoG = RhoG(Elements)';

% eWa = eWa(Elements)';
ha = ha(Elements)';
hb = hb(Elements)';
hvN = hv;
hv = hv(Elements)';
hwN = hw;
hw = hw(Elements)';

%% Applying rotation for RTL -> XYZ
% Work out a better way for hadling rotations in a more general way. This
% is hard coded for mid board X <-> Z rotation. Keeping the normal vecotors
% the same sign.


% Bound Diffusivity
[DbX,DbY,DbZ,DbXY,DbXZ,DbYZ] = ActionRotationVectors(EWPMesh.RotationVectors, ...
    Properties.DbR,Properties.DbT,Properties.DbL);

DbX = DbX(Elements)'; DbY = DbY(Elements)'; DbZ = DbZ(Elements)';
DbXY = DbXY(Elements)'; DbXZ = DbXZ(Elements)'; DbYZ = DbYZ(Elements)';
DBXc = DbX(Elements);

% % Effective Thermal Conductivity
[KeffX,KeffY,KeffZ,KeffXY,KeffXZ,KeffYZ] = ActionRotationVectors(EWPMesh.RotationVectors, ...
    Properties.KeffR,Properties.KeffT,Properties.KeffL);
KeffX = KeffX(Elements)'; KeffY = KeffY(Elements)'; KeffZ = KeffZ(Elements)';
KeffXY = KeffXY(Elements)'; KeffXZ = KeffXZ(Elements)'; KeffYZ = KeffYZ(Elements)';

% % Effective vapour diffusivity
[DeffX,DeffY,DeffZ,DeffXY,DeffXZ,DeffYZ] = ActionRotationVectors(EWPMesh.RotationVectors, ...
    Properties.DeffR,Properties.DeffT,Properties.DeffL);
DeffX = DeffX(Elements)'; DeffY = DeffY(Elements)'; DeffZ = DeffZ(Elements)';
DeffXY = DeffXY(Elements)'; DeffXZ = DeffXZ(Elements)'; DeffYZ = DeffYZ(Elements)';

% % Relative permeabilities;
% [kwX, kwY, kwXY] = ComputeRotationRT2XY(Properties.kwR, Properties.kwT, NodalGrainAngle);
% [kwX,kwY,kwZ,kwXY,kwXZ,kwYZ] = EWPRotateSandwich(NodeRotations,kwX,kwY,Properties.kwL,kwXY);
[kwX,kwY,kwZ,kwXY,kwXZ,kwYZ] = ActionRotationVectors(EWPMesh.RotationVectors, ...
    Properties.kwR, Properties.kwT,Properties.kwL);
kwX = kwX(Elements)'; kwY = kwY(Elements)'; kwZ = kwZ(Elements)';
kwXY = kwXY(Elements)'; kwXZ = kwXZ(Elements)'; kwYZ = kwYZ(Elements)';

% % Absolute Permeabilities
[KX,KY,KZ,KXY,KXZ,KYZ] = ActionRotationVectors(EWPMesh.RotationVectors, ...
    Properties.KR,Properties.KT,Properties.KL);
KX = KX(Elements)'; KY = KY(Elements)'; KZ = KZ(Elements)';
KXY = KXY(Elements)'; KXZ = KXZ(Elements)'; KYZ = KYZ(Elements)';



%% Gradient computation
% Compute the gradient of the bound water for each element



[GradXbX, GradXbY, GradXbZ] = TPEShapeFuncGradVectorised(A11, A12, A13, ...
    A21, A22, A23, A31, A32, A33, Xb, Elements);
[GradWvX, GradWvY, GradWvZ] = TPEShapeFuncGradVectorised(A11, A12, A13, ...
    A21, A22, A23, A31, A32, A33, wv, Elements);
PressureHead = Pw./RhoWGrav + NodePos(:,3); % Gravity check
[GradPwX, GradPwY, GradPwZ] = TPEShapeFuncGradVectorised(A11, A12, A13, ...
    A21, A22, A23, A31, A32, A33, PressureHead, Elements);
% gravity check
% Total head is H = Pw./RhoWGrav - NodePos(:,3);
[GradTX, GradTY, GradTZ] = TPEShapeFuncGradVectorised(A11, A12, A13, ...
    A21, A22, A23, A31, A32, A33, T, Elements);

wv = wv(Elements)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   Internal Flux Loop    %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Do we wanna use mex (usually yes)? using the codegen at the start of the
% EEM script
if UseMex == 1

    [g] = FVMTPETransPoreElementLoop_mex(NumNodes, NumElements, Elements, ...
        NormalVectorsX, NormalVectorsY, NormalVectorsZ, ElementRho0, ...
        ShapeFuncArrayT, DbX, DbY, DbZ, DbXY, DbXZ, DbYZ, DeffX, DeffY, ...
        DeffZ, DeffXY, DeffXZ, DeffYZ, KeffX, KeffY, KeffZ, KeffXY, KeffXZ, ...
        KeffYZ, KX, KY, KZ, KXY, KXZ, KYZ, kwX, kwY, kwZ, kwXY, kwXZ, kwYZ, ...
        muw, RhoG, wv, ha, hb, hv, hw, PressureHead, RhoWGrav, RhoW, GradPwX, ...
        GradPwY, GradPwZ, GradXbX, GradXbY, GradXbZ, GradWvX, GradWvY, GradWvZ, ...
        GradTX, GradTY, GradTZ, SubAreas);
else
    [g] = FVMTPETransPoreElementLoop(NumNodes, NumElements, Elements, ...
        NormalVectorsX, NormalVectorsY, NormalVectorsZ, ElementRho0, ...
        ShapeFuncArrayT, DbX, DbY, DbZ, DbXY, DbXZ, DbYZ, DeffX, DeffY, ...
        DeffZ, DeffXY, DeffXZ, DeffYZ, KeffX, KeffY, KeffZ, KeffXY, KeffXZ, ...
        KeffYZ, KX, KY, KZ, KXY, KXZ, KYZ, kwX, kwY, kwZ, kwXY, kwXZ, kwYZ, ...
        muw, RhoG, wv, ha, hb, hv, hw, PressureHead, RhoWGrav, RhoW, GradPwX, ...
        GradPwY, GradPwZ, GradXbX, GradXbY, GradXbZ, GradWvX, GradWvY, GradWvZ, ...
        GradTX, GradTY, GradTZ, SubAreas);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   Boundary Conditions    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Psychometric

% tend = 4 * 24 * 60 * 60; % Days x hours x minutes x seconds
% x = Time/3600/24;
Pbar = 101325;
Patm = Pbar;
% Tdb = 26.64355 + 5.166129*sin(Time./(86400/2/pi)); % 30;
a = pTdb(1); b = pTdb(2); d = pTdb(3);
Tdb = a + b*sin((Time./3600)*(2*pi)/24-d);
% Tdb = 25 + 5*sin(Time./(86400/2/pi)); % 30;
a = pRelHum(1); b = pRelHum(2); d = pRelHum(3);
RelHum = a + b*sin((Time./3600)*(2*pi)/24-d);  % Needs to be a percentage here
% % % Twb = 25;% + 5*sin(Time./(86400/2/pi)); % 28;
SatVapPres = @(T) exp(-5.8002206e3 ./ T + 1.3914993 - 4.8640239e-2 .* T...
    + 4.1764768e-5 .* T.^2 - 1.4452093e-8 .* T.^3 + 6.5459673 .* log(T));
% SatHumRatio = @(SVP)0.621945 .* SVP ./ (Patm - SVP);
SatVapPres_db = SatVapPres(Tdb + Tr);
% Given relative humidity, calculate Twb
VapPres = (RelHum/100) * SatVapPres_db;
% HumRatio = 0.621945 * VapPres / (Patm - VapPres);
% % From https://journals.ametsoc.org/view/journals/apme/50/11/jamc-d-11-0143.1.xml
Twb = Tdb * atan(0.151977 * (RelHum + 8.313659)^0.5) + atan(Tdb + RelHum)...
    - atan(RelHum - 1.676331) + 0.00391838 * RelHum^1.5 * atan(0.023101 * RelHum)...
    - 4.686035;

%% Tris


c = Patm / (R*(Tdb + Tr));
% Triangle Faces
for i = 1:length(Triangles)
    Nodes = Triangles(i,:);
    TriFace = TrianglesIndex(i);
    BCType = BoundaryConds(TriFace);

    if BCType == 0 % No flux
        WatFlux = 0;
        TFlux = 0;
    elseif BCType == 1 % Atmospheric
        xv = Pv(Nodes) ./ Patm;
        xvinf = VapPres / Patm;
        WLogTerm = log((1 - xvinf) ./ (1 - xv));
        WatFlux = -km * c * Mv * WLogTerm;

        TFlux = -q* (T(Nodes) - Tdb) + 1*hvN(Nodes) .*WatFlux;

    elseif BCType == 2 % Saturating

        WatFlux = -SaturatingParam.*(X(Nodes) - XSat(Nodes));
        TFlux = -HglueT* (T(Nodes) - Twb);

    else
        WatFlux = 0;
        TFlux = 0;

    end

    switch TriFace
        case 1 % X = 0
            % Remember TriSCVArea is 1/3 the area of the triangle

            g(Nodes*2-1) = g(Nodes*2-1) + WatFlux * TriSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFlux * TriSCVArea(i);
        case 2 % X = Xmax
            g(Nodes*2-1) = g(Nodes*2-1) + WatFlux * TriSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFlux * TriSCVArea(i);
        case 3 % Y = 0
            g(Nodes*2-1) = g(Nodes*2-1) + WatFlux * TriSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFlux * TriSCVArea(i);
        case 4 % Y = Ymax
            g(Nodes*2-1) = g(Nodes*2-1) + WatFlux * TriSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFlux * TriSCVArea(i);

        case 5 % Z = 0
            g(Nodes*2-1) = g(Nodes*2-1) + WatFlux * TriSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFlux * TriSCVArea(i);

        case 6 % Z = Zmax
            g(Nodes*2-1) = g(Nodes*2-1) + WatFlux * TriSCVArea(i);
            % TFlux = -1* (T(Nodes) - Twb);
            g(Nodes*2) = g(Nodes*2) + TFlux * TriSCVArea(i);
            %
            case 7
                % Glue line interaction layer, implement this later
                error("Glue line interaction on triangular face! Not good (implemented)!")

        otherwise
            error("Tris surface index not good")
    end


end
%
%
%
% %% %% getting surface transfer coeff sigma_w
% SumXiXSatArea = 0;
% QArea = 0;
% % PhiAvgRho0(n) = sum(NodeRho0 .* X .* FullVolumes) / sum(NodeRho0 .* FullVolumes);
% XDiff = zeros(length(Quads),1);
% SAreas = zeros(length(Quads),1);
% SRho0 = zeros(length(Quads),1);
% SVols = zeros(length(Quads),1);
% for i = 1:length(Quads)
% 
%     Nodes = Quads(i,:);
%     QuadFace = QuadsIndex(i);
%     BCType = BoundaryConds(QuadFace);
%     % WFlux = -qw * ((X(Nodes)-Properties.Xb(Nodes)) - XSat(Nodes));
%     % Fix this hack later !!!!!!! add + FSP to XSat
%     if BCType == 2 % Saturating
% 
%         % SumXiXSatArea = SumXiXSatArea + sum((X(Nodes) - XSat(Nodes)) * 4 * QuadSCVArea(i));
%         SumXiXSatArea = SumXiXSatArea + sum((Properties.Sw(Nodes) - 1) * 4 * QuadSCVArea(i));
%         % QArea = QArea + 4 * QuadSCVArea(i);
%         % XDiff(i) = mean(X(Nodes) - XSat(Nodes));
%         % SRho0(i) = mean(NodeRho0(Nodes));
%         % SAreas(i) = 4*QuadSCVArea(i);
%     end
% 
% end
% % SumXiXSatArea2 = sum(SRho0 .* XDiff .* SAreas) / sum(SRho0 .* SAreas)*sum(SAreas);
% % % Values from curve fitter
% % XAvgdt = 8.28e-2./((Time/3600)+1.0331)/3600;
% % XAvgdt = 6.8050e-05/(1.0000e-04*Time + 1.9544)^2;
% % XAvgdt2 = 0.7589/((0.5965*(Time/3600) + 3.6638)^2);
% % XAvgdt = 3.6364e-04/(2.0000e-04*Time + 5.2220)^2;
% XAvgdt = 0.0829/(Time + 3.7531e+03);
% Rho0Vol = sum(NodeRho0 .* FullVolumes);
% SaturatingParam = min((SaturatingParam),0.01);
FluxDiff = SaturatingParam;
% sigma_w = 0*min(sigma_w,50);
%% Quad faces

for i = 1:length(Quads)

    Nodes = Quads(i,:);

    QuadFace = QuadsIndex(i);
    BCType = BoundaryConds(QuadFace);

    if BCType == 0 % No flux
        
        WatFluxQ = 0;
        TFluxQ = 0;
    elseif BCType == 1 % Atmospheric
        xvQ = Pv(Nodes) ./ Patm;
        xvQinf = VapPres / Patm;
        WLogTermQ = log((1 - xvQinf) ./ (1 - xvQ));
        WatFluxQ = -km * c * Mv * WLogTermQ;

        TFluxQ = -q* (T(Nodes) - Tdb) + hvN(Nodes) .*WatFluxQ;

    elseif BCType == 2 % Saturating
        % WatFlux = min(-SaturatingParam*(X(Nodes) - XSat(Nodes)), SaturatingParam);

        WatFluxQ = -1*SaturatingParam.*(X(Nodes) - XSat(Nodes));
        TFluxQ = -HglueT* (T(Nodes) - Twb);
    else
        WatFluxQ = 0;
        TFluxQ = 0;

    end
    switch QuadFace
        case 1 % X = 0
            % Remember QuadSCVArea is 1/4 the area of the quad
            g(Nodes*2-1) = g(Nodes*2-1) + WatFluxQ * QuadSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFluxQ * QuadSCVArea(i);

        case 2 % X = XMAX
            g(Nodes*2-1) = g(Nodes*2-1) + WatFluxQ * QuadSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFluxQ * QuadSCVArea(i);

        case 3 % Y = 0
            g(Nodes*2-1) = g(Nodes*2-1) + WatFluxQ * QuadSCVArea(i);
            g(Nodes*2) = g(Nodes*2) +  (-q*(T(Nodes) - Twb)) * QuadSCVArea(i);

        case 4 % Y = YMAX
            g(Nodes*2-1) = g(Nodes*2-1) + WatFluxQ * QuadSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFluxQ * QuadSCVArea(i);

        case 5 % Z = 0
            g(Nodes*2-1) = g(Nodes*2-1) + WatFluxQ * QuadSCVArea(i);
            g(Nodes*2) = g(Nodes*2) + TFluxQ * QuadSCVArea(i);

        case 6 % Z = ZMAX
            % g(Nodes*2-1) = g(Nodes*2-1) + 0.00001 * QuadSCVArea(i);
            % g(Nodes*2) = g(Nodes*2)  + 2 * QuadSCVArea(i);
            g(Nodes*2-1) = g(Nodes*2-1) + WatFluxQ * QuadSCVArea(i);
            g(Nodes*2) = g(Nodes*2)  + TFluxQ * QuadSCVArea(i);

        case 7
            % Glue line interaction layer
            for node = 1:4 % for each node on the otherside of the quad
                % try
                    OtherQuadNodes = Quads(GlueInteractions(Nodes(node),1),:);
                % catch
                %     1
                % end
                XQuads = X(OtherQuadNodes);
                XSatQuads = XSat(OtherQuadNodes);
                TQuads = T(OtherQuadNodes);

                QNode = Nodes(node);

                z = GlueInteractions(QNode,2);
                n = GlueInteractions(QNode,3);

                XInterp = 1/4*(XQuads(1)*(1-z)*(1-n) + XQuads(2)*(1+z)*(1-n) ...
                    + XQuads(3)*(1+z)*(1+n) + XQuads(4)*(1-z)*(1+n));

                XSatInterp = 1/4*(XSatQuads(1)*(1-z)*(1-n) + XSatQuads(2)*(1+z)*(1-n) ...
                    + XSatQuads(3)*(1+z)*(1+n) + XSatQuads(4)*(1-z)*(1+n));

                TInterp = 1/4*(TQuads(1)*(1-z)*(1-n) + TQuads(2)*(1+z)*(1-n) ...
                    + TQuads(3)*(1+z)*(1+n) + TQuads(4)*(1-z)*(1+n));

               % Ensuring no oversaturation
                Xci = min(X(QNode),XSatInterp);
                Xcab = min(XInterp,XSat(QNode));
                WatFluxG = -HglueX*(Xci - Xcab);


                TFluxG = -HglueT*(T(QNode) - TInterp) - hwN(QNode) .* WatFluxG;

                g(QNode*2-1) = g(QNode*2-1) + WatFluxG  * QuadSCVArea(i);
                g(QNode*2) = g(QNode*2) + TFluxG  * QuadSCVArea(i);
            end

        otherwise
            error("Quads surface index not good")
    end


end



%% Finite difference approximation

h = sqrt(eps) .* norm(XTheta);
% Compute Accumulation Terms and their shifted values for X and T
[PsiW,PsiE] = ComputePsiWE(X,Theta*Tbp,NodeRho0,Constants); % Standard
[PsiWX,PsiEX] = ComputePsiWE(X+h,Theta*Tbp,NodeRho0,Constants);
[PsiWT,PsiET] = ComputePsiWE(X,(Theta+h)*Tbp,NodeRho0,Constants);

% compute derivatives
dPsiWdX = (PsiWX - PsiW)/h;
dPsiWdT = (PsiWT - PsiW)/h;

dPsiEdX = (PsiEX - PsiE)/h;
dPsiEdT = (PsiET - PsiE)/h;

JinvFactor = 1 ./ (dPsiWdX .* dPsiEdT - dPsiWdT .* dPsiEdX);


% Change of variables for EEM
for node = 1:NumNodes % X; Theta
    g(node*2-1) = JinvFactor(node) * (dPsiEdT(node) * g(node*2-1) - dPsiWdT(node) * g(node*2));
    g(node*2) = JinvFactor(node) * (-dPsiEdX(node) * g(node*2-1) + dPsiWdX(node) * g(node*2));

end

% FVM CV scaling
g = g ./ (InterlacedFullVolumes);