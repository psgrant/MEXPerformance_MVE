%% Defining offsets and roations
% We're going to be defining everything in terms of the bottom board. GMSH
% extrudes in the Z direction so the z = 0 face is the orginal 2d mesh
% we've generated.
% clear all
clearvars -except JobCounter
ITime = tic;



% Change from QST1.mat - QST6.mat
MeshStrCell = {'QST3.mat'};


UpdateMesh = 1; 
if exist('EWPMesh','var') == 0 || UpdateMesh == 1
    [EWPMesh] = OverlappingSurfacePoints(MeshStrCell,EWPOffset,EWPRotation);
end


% 1 for mex (function will need to be built if it is not already)
EWPMesh.UseMex = 0;
BTime = tic;
Bounded = 1; % unbounded (infxtypeof) = 0, bounded = 1
if EWPMesh.UseMex == 1
        fprintf('Building Mex file for the element loop!\n')
    BuildFVMElementLoop(EWPMesh,Bounded)
end
BTime = toc(BTime);


ConsolePrint = 1; % display in command window progress
PlottingInterval = 10; % When to update the plot

%% Stuff for parameter fitting not important rn
bdiff = linspace(40,250,16); % Diff
S = logspace(-2,1,16); % Perm
% SatParamArray = (logspace(-4,-2,12));
SatParamArray = 1;%(logspace(-4,-2,12));
[B, S] = meshgrid(bdiff,S);
[r,c] = ind2sub(size(B),JobCounter);
%% set up saving dirs for multiple sims
% homeDir = getenv('HOME');  % Get the home directory path
% % homeDir = pwd; % Get the home directory path
% % DirStr = fullfile(homeDir, ['/Codes/QSIngress_Face_LinFunc/Sim_',...
% %     num2str(JobCounter),'_SatParam=',num2str(SatParamArray(JobCounter),'%.2e')]);  % Specify the folder name
% DirStr = fullfile(homeDir, ['/Codes/QSF_LargeLog/Sim_',...
%     sprintf('%04.0f_B%g_S%g',JobCounter,r,c)]);  % Specify the folder name
%
% if ~exist(DirStr, 'dir')
%        mkdir(DirStr)
% end

%% Build mesh
% Define offests
% These are the offsets for the board, note that the rotation axis are
% centred at 0 so we need to counteract that with a translation eithqdel er
% one of the two other directions.
MidOffset = [0 0.065/2 -0.065/2];
% MidOffset = [-0.04/2 0.040/2 0];
BotOffset = [0 0.045/2 0.045/2];
EWPOffset = [MidOffset; BotOffset];

% Rotations are in radians and are clockwise about each axis
% centered at X/Y/Z = 0
MidRotation = [0.5 0.5 0.5]*pi; % clockwise rotation about % Y axis
% MidRotation = [0 0 0.5]*pi; % clockwise rotation about % Y axis
BotRotation = [0 0 0];
EWPRotation = [MidRotation; BotRotation];


%
% clf
% map = flipud(copper(256));
% map = (map(1:end-30,:));
% ElementCol = map(round(rescale(EWPMesh.NodeRho0)*(length(map)-1))+1,:);
% Plot_Mesh3D(NodePos,Elements,ElementCol,1,1)
% colormap(map)
% clim(minmax(EWPMesh.ElementRho0'))
% colorbar
% view([-45 25])
% xlabel('X Direction [Metres]')
% ylabel('Y Direction [Metres]')
% zlabel('Z Direction [Metres]')
% % axis([0 0.09 0 0.09 0 0.09])
% % pbaspect([1 1 1])
% axis equal
% drawnow
%


% Set bounbdary conditions for each external face
% 0 = no flux, 1 = atmospheric, 2 = saturating, 3 = Glue
XminBC = 0; % idx = 1
XmaxBC = 0; % idx = 2
YminBC = 0; % idx = 3
YmaxBC = 0; % idx = 4
ZminBC = 0; % idx = 5
ZmaxBC = 2; % idx = 6

GlueBC = 3; % idx = 7
EWPMesh.BoundaryConds = [XminBC, XmaxBC, YminBC, YmaxBC, ZminBC, ZmaxBC, GlueBC];
BoundaryConds = [XminBC, XmaxBC, YminBC, YmaxBC, ZminBC, ZmaxBC, GlueBC];
% y1 = rand([1,18])*0.9;
% y2 = rand([1,18])*0.9;
% y1y2 = [y1;y2];


% EWPMesh.SaturatingParam = S(JobCounter);
% EWPMesh.y1y2 = y1y2(:,JobCounter);
EWPMesh.HglueX = 1e-4;
EWPMesh.HglueT = 50;

% Atmoshpheric conditions
% p(1) + p(2)*sin((x/3600)*(2*pi)/24-p(3))
% Set constant with [X 0 0]
% EWPMesh.pTdb = [26.6394, 2.4622, 20.8259];
% EWPMesh.pRelHum = [62.2095, 12.1179 5.0825];
EWPMesh.pTdb = [25, 0, 0];
EWPMesh.pRelHum = [62, 0 0];

NodePos = EWPMesh.NodePos;
Elements = EWPMesh.Elements;
SurfTriangles = EWPMesh.Triangles;
SurfQuads = EWPMesh.Quads;
ScalingFacts = EWPMesh.ScaleFacts;
NormalVectorsX = EWPMesh.NormalVecsX;
NormalVectorsY = EWPMesh.NormalVecsY;
NormalVectorsZ = EWPMesh.NormalVecsZ;
CenterPos = EWPMesh.CenterPos;
SubVolumes = EWPMesh.SubVolumes;
SubAreas = EWPMesh.SubAreas;
FullVolumes = EWPMesh.FullVolumes;
SurfTriArea = EWPMesh.SurfTriArea;
SurfQuadArea = EWPMesh.SurfQuadArea;
TrianglesIndex = EWPMesh.TrianglesSurfIndex;
Quads = EWPMesh.Quads;
QuadsIndex = EWPMesh.QuadsSurfIndex;
NodeRho0 = EWPMesh.NodeRho0;
ElementRho0 = EWPMesh.ElementRho0;

NodalGrainAngle = EWPMesh.NodalGrainAngle;
%%
EWPMesh.ShapeFuncArray = ComputeShapeFuncArray();
EWPMesh.RotationVectors = EWPRotationVectors(EWPMesh.NodeRotations, ...
    NodalGrainAngle);

%%

% EWPMesh.NodeRotations = zeros(length(EWPMesh.Elements),3);

% Compute inverse Jacobian elements for shape functions
tic
% [IP1, IP2, IP3, IP4, IP5, IP6, IP7, IP8, IP9] = TPEShapeFuncGradJ(NodePos, Elements);
[A11, A12, A13, A21, A22, A23, A31, A32, A33] = TPEShapeFuncGradJinv(NodePos, Elements);
toc
% Load in mesh data
EWPMesh.A11 = A11;
EWPMesh.A12 = A12;
EWPMesh.A13 = A13;
EWPMesh.A21 = A21;
EWPMesh.A22 = A22;
EWPMesh.A23 = A23;
EWPMesh.A31 = A31;
EWPMesh.A32 = A32;
EWPMesh.A33 = A33;

EWPMesh.SigmaSave = 0;

SigmaSave = zeros(1000,1);

% EWPMesh.IP1 = IP1;
% EWPMesh.IP2 = IP2;
% EWPMesh.IP3 = IP3;
% EWPMesh.IP4 = IP4;
% EWPMesh.IP5 = IP5;
% EWPMesh.IP6 = IP6;
% EWPMesh.IP7 = IP7;
% EWPMesh.IP8 = IP8;
% EWPMesh.IP9 = IP9;

%% Specify constants
Constants.RhoS = 1530;           % kg/m^3    - Solid phase density
Constants.RhoW = 1000;           % kg/m^3    - Density of water
Constants.RhoFBar = 0.002;       % kg/mol    - Superficial density (fixed phase)
Constants.Patm = 101325;         % Pa -      - Atmospheric pressure
Constants.CPw = 4187;            % J/(kg K)  - Specific heat of water
Constants.CPs = 1400;            % J/(kg K)  - Specific heat of solid
Constants.CPa = 1005.683;        % J/(kg K)  - Specific heat of air
Constants.CPv = 1900;            % J/(kg K)  - Specific heat of vapor
Constants.h0vap = 2.503e6;       % J/kg      - Latent heat of vaporization
Constants.R = 8.314;             % J/(mol K) - Universal gas constant
Constants.Mv = 0.018016;         % kg/mol    - Molecular weight of vapor
Constants.Ma = 0.028952;         % kg/mol    - Molecular weight of air
Constants.Tr = 273.15;           % K         - Reference temperature
Constants.Ti = 25;               % K         - Initial temperature
Constants.km = 0.025;            %           - mass transfer coeff
Constants.q = 25;
Constants.MeshVolume = prod(max(NodePos));
% Bound water diffusivity and permeabilities
% scaling factors
Constants.DiffFact = B(JobCounter);
Constants.PermFact = S(JobCounter);
Tbp = 100;

% Face and edge alpha for plotting
FAlpha = 1;
EAlpha = 1;
PlottingTriangles = SurfTriangles(:,1:3);
PlottingQuads = SurfQuads(:,1:4);

% Setting time params for simu
t = 0;
Tn = 0.1;
TnMax = 3600*100;
TMax = 3600*96; % Days x hours x minutes x seconds
SaveTimeInterval = 10*60;
CurrentSaveTime = 0;

Constants.Time = t;

VapourAverage = zeros(1000,1);
PhiAvgRho0 = zeros(1000,1);
PhiAvgRho0T = zeros(1000,1);
KrylovArray = zeros(1000,1);
TnArray = zeros(1000,1);


% get the number of nodes and elements
NumNodes = length(NodeRho0);
NumNodes2 = NumNodes * 2;
NumElements = length(ElementRho0);


% Phi Pade path
folderName = 'Pade Code';  % Replace 'subfolder' with the name of your folder
currentDir = pwd;  % Get the current directory
folderPath = fullfile(currentDir, folderName);  % Create the full folder path
addpath(folderPath);  % Add the folder to the MATLAB path

% % Set X and T
XFSP = 0.3;

P = 1 - NodeRho0/1530;
% X0 = (P .* 0.75 * Constants.RhoW) ./ (NodeRho0);
% XSat = (1 - NodeRho0/1530) * 1000 ./(NodeRho0);
%

X0 = 0.12 * ones(NumNodes,1); % +0.3*NodePos(:,1)/0.09;
% X0(1:EWPMesh.NumMidNodes) = X0(1:EWPMesh.NumMidNodes)  + 1;
% X0 = X0  + 0.25*NodePos(:,2)/0.09 + 0.15*NodePos(:,3)/0.09;
X = X0;
T0 = 26.5*ones(NumNodes,1);
% T0(1:EWPMesh.NumMidNodes) = T0(1:EWPMesh.NumMidNodes)  + 10;
T0 = T0/Tbp;
Theta = T0;
[Properties] = ComputeMaterialProperties(X,Theta*Tbp,NodeRho0,Constants);
EWPMesh.XSat = Properties.XSat;
XSat = Properties.XSat;

Rho0Vol = sum(NodeRho0.*FullVolumes);

% XAvgdt = 0.0829/(t + 3.7531e+03);
% a = 0.4538;
% b = -1.7351e+04;
% a = 0.0829;
% b = -3.7531e+03;
% a = 0.5718;
% b = -2.3059e+04;
% XAvgdt = a/(t - b);

% b = 9.7621e+04;
% c = 6.1910e+04;
b = 7.7849e+03;
c = 2.2358e+04;
% b = 1.116e+04;
% c = 1.75e+04;
XAvgdt = b./(c + t).^2;
[SumXiXSatArea] = ComputeSurfaceFlux(Quads,QuadsIndex, ...
    BoundaryConds,SurfQuadArea,X,XSat);
EWPMesh.sigma_w = -Rho0Vol * XAvgdt / SumXiXSatArea;
SigmaSave(1) = EWPMesh.sigma_w ;
%%
XTheta = zeros(NumNodes2,1);
XTheta(1:2:NumNodes2) = X;
XTheta(2:2:NumNodes2) = Theta;
XTheta = real(XTheta);



% Weighting Vector tolerances
Atol = 1e-8;
Rtol = 1e-8;
C = 0.25;
q = 2;
Eta = 0.9;
n = 1;
Acounter = 0;

% XSat = (1 - NodeRho0/1530) * 1000 ./(NodeRho0);
% XSat = (1 - NodeRho0/1530) * 1000 ./(NodeRho0) - (Xfsp + 0.001);
% X = XSat-0.001;

epsG = Properties.epsG;
RhoV = Properties.RhoV;
PhiAvgRho0(n) = sum(NodeRho0 .* X0 .* FullVolumes) / sum(NodeRho0 .* FullVolumes);
PhiAvgRho0T(n) = sum(NodeRho0 .* T0 .* FullVolumes * Tbp) / sum(NodeRho0 .* FullVolumes);
VapourAverage(n) = sum(epsG .* RhoV .* FullVolumes) / sum(NodeRho0 .* FullVolumes);



%% Creating simulation plotting enviroment
% Fig 1 simulations
% close all
% figure('Name','Simulation Plots')
% pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');% #ok this code might not work in future
% % You can remove it though as it just maximised the figure window
% set(frame_h,'Maximized',1);
% pause(0.01)

% Save Interval
SaveCounter = 0;
SaveCSVCounter = 1;
SaveInterval = 100;



%% Optional data from previous sim

LoadData = 0; % Calls in mesh
if LoadData == 1
    clear all %#ok
    load([pwd,'\\SimulationCSVResults\\FixedGLSatBC','\\SimulationMat.mat'])
else
    SaveCounter = 1;
    SaveCSVCounter = 1;
    SaveInterval = 100;
    TimeArraySave  = zeros(1,SaveInterval);
    TimeArray = zeros(1000,1);
    XArray = zeros(length(XTheta),SaveInterval);
    X0T0 = [X0'; T0'];
    XArray(:,1) = X0T0(:);
    CurrentSaveTime = CurrentSaveTime + SaveTimeInterval;
    AlphaLimit = 1000 ;
    % Perterbation amount order e-8
    Pert = sqrt(eps);
    MaxKrylov = 500;
    MinKrylov = 1;
    % Weighting Vector tolerances
    Atol = 1e-4;
    Rtol = 1e-4;
    C = 0.25;
    q = 2;
    Eta = 0.9;
    n = 1;

end
% save([DirStr,'/PreSimMeshData.mat'])
initTime = toc(ITime);
initTime = initTime - BTime;
%%
% Start Loop
Stime = tic;
while t < TMax

    % for n = 1:2500
    T1 = tic;

    % Compute Vm and Hm using Arnoldi's Method Modified Gram-Schmidt
    [g0, ~] = TransPoreFluxRd3(EWPMesh,Constants,XTheta);


    % error('hehe')
    r0 = g0;
    beta0 = norm(g0,2);

    % Initalise the V and H matricies
    V = zeros(NumNodes2,MaxKrylov + 1);
    H = zeros(MaxKrylov + 1);


    % First iteration in V
    if beta0 == 0
        V(:,1) = r0;
    else
        V(:,1) = r0 / beta0;
    end
    NeedMoreKrylov = 1;
    m = 1;

    %     Arnoldi's Method Modified Gram-Schmidt
    while NeedMoreKrylov == 1
        % Apply perturbation to compute first part of Vm+1
        % VPert = Pert * V(:,m);
        % [g1, FluxDiff] = TransPoreFluxRd3(EWPMesh, Constants, XTheta + Pert * V(:,m));
        [g1, FluxDiff] = TransPoreFluxRd3(EWPMesh, Constants, XTheta + Pert * V(:,m));
        V(:,m+1) = (g1 - g0)/ Pert;
        % Update hessenberg
        for j = 1:m
            H(j,m) = V(:,j)' * V(:,m+1);
            V(:,m+1) = V(:,m+1) - H(j,m) * V(:,j);
        end
        betam = norm(V(:,m+1),2);
        H(m+1,m) = betam;
        if abs(betam) < 1e-14
            % fprintf("Enough Krylov Basis Vectors Found at m = %g \n",m)
            break
        else % Compute Vm+1
            Vm1= V(:,m+1) / H(m+1,m);
            V(:,m+1) = Vm1;
        end

        % Matrix expoential
        PhiExpoFull = phipade(Tn * H(1:m,1:m), 1);

        % Selection of Krylov dimension
        pm = Tn * beta0 * betam * PhiExpoFull(m,1) * Vm1;
        WeightVec = Rtol * abs(pm) + Atol;
        pmTnNWRMS = Tn*(1 / NumNodes * sum((pm ./ WeightVec).^2))^0.5;

        % Check to see if the new additons of spaces is needed
        if pmTnNWRMS < 1 && m >= MinKrylov
            NeedMoreKrylov = 0; % Got enough Krylov
        else
            m = m + 1; % More Krylov
        end

        if m >= (MaxKrylov) % Half timestep if we reach the max krylov iters
            m = m - 1;
            while pmTnNWRMS > 1
                Tn = 0.5 * Tn;
                % Matrix expoential
                PhiExpoFull = phipade(Tn * H(1:m,1:m), 1);

                % Selection of Krylov dimension
                pm = Tn * beta0 * betam * PhiExpoFull(m,1) * Vm1;
                WeightVec = Rtol * abs(pm) + Atol;
                pmTnNWRMS = Tn*(1 / NumNodes * sum((pm ./ WeightVec).^2))^0.5;
                % fprintf(2,"Max Krylov dimension reached - half time step Tn\n")
                if ConsolePrint == 1
                    ffprintf(2,"Max Krylov dimension reached - half time step Tn\n")
                end
            end
            NeedMoreKrylov = 0;

        end

    end
    Vm = V(:,1:m);

    e1 = [1; zeros(m-1,1)];

    % Make sure solution stays physical, reduce timestep if not


    PhiExpoHalf = phipade(0.5 * Tn * H(1:m,1:m), 1);%(0.5 * Tn * Hm) \ (expm(0.5 * Tn * Hm) - eye(m));

    % Computing the 3 steps of the EEM (1 full and 2 half steps)


    % Compute the 3 steps
    EEMFull = XTheta + Tn * beta0 * Vm * PhiExpoFull(1:m,1:m) * e1;
    EEMHalf = XTheta + 0.5 * Tn * beta0 * Vm * PhiExpoHalf(1:m,1:m) * e1;



    [SwHalf] = ComputeSaturation(EEMHalf(1:2:end),EEMHalf(2:2:end),NodeRho0,Constants);
    [SwFull] = ComputeSaturation(EEMFull(1:2:end),EEMFull(2:2:end),NodeRho0,Constants);

    % make sure solution stays physical
    if (any(SwHalf > 1) || any(SwFull > 1)) ||...
            (any(isnan(SwHalf)) || any(isnan(SwFull)))
        Tn = 0.5 * Tn;
        if ConsolePrint == 1
            fprintf(2,"Solution non-physical half time step Tn\n")
        end
        %
        continue
    end

    [gHalf,~] = TransPoreFluxRd3(EWPMesh,Constants,EEMHalf);

    % make sure solution stays physical
    if any(~isreal(gHalf))

        Tn = 0.5 * Tn;
        if ConsolePrint == 1
            fprintf(2,"Solution non-physical half time step Tn\n")
        end

        continue
    end
    EEMFull2Step = EEMHalf + 0.5 * Tn * Vm * (PhiExpoHalf(1:m,1:m) * Vm' * gHalf);
    % Calculate error
    DeltaN = EEMFull - EEMFull2Step;
    WeightVec = Rtol * abs(XTheta) + Atol;
    DeltaNWRMS = (1 / NumNodes * sum((DeltaN ./ WeightVec).^2))^0.5;

    % Determine stepsize control factor
    Alpha = Eta * (1 / DeltaNWRMS) ^ (1/(q+1));


    % Successfull step
    if C * DeltaNWRMS <= 1 % add in check for non physicality

        T1End = toc(T1);
        if ConsolePrint == 1
            fprintf("Step = %g : m = %g: Alpha = %.2f : Tn = %.4fs : t = %.4fh : Step Time = %d Mins and %.3f Seconds\n"...
                ,n,m,Alpha,Tn,t/3600,floor(T1End/60), rem(T1End,60))
        end
        % Limit alpha
        Alpha = min(Alpha,AlphaLimit);
        TimeArray(n) = t;


        t = min(t + Tn,TMax);


        % Compute sigma_w
        XAvgdt = b./(c + t).^2;
        [SumXiXSatArea] = ComputeSurfaceFlux(Quads,QuadsIndex, ...
            BoundaryConds,SurfQuadArea,X,XSat);
        EWPMesh.sigma_w = -Rho0Vol * XAvgdt / SumXiXSatArea;
        Constants.Time = t;

        % Update function
        XTheta = EEMFull2Step;
        KrylovArray(n) = m;
        TnArray(n) = Tn;

        X = XTheta(1:2:NumNodes2);
        Theta = XTheta(2:2:NumNodes2);
        if t + Tn > CurrentSaveTime
            SaveCounter = SaveCounter + 1;
            TimeArraySave(SaveCounter) = t;
            XArray(:,SaveCounter) = XTheta;
            CurrentSaveTime = ceil((Tn + t - CurrentSaveTime)/SaveTimeInterval)...
                * SaveTimeInterval + CurrentSaveTime;
        end
        % Saving for HPC parameter testing not needed right now
        % if SaveCounter == SaveInterval% Save data and stuff
        %     fprintf("Writting data to CSV File :::::: SimulationData_%g.csv\n",SaveCSVCounter)
        %     SaveStr = [DirStr, sprintf('/SimulationData_%g.csv',SaveCSVCounter)];
        %     save([DirStr,'/SimulationWorkspace.mat'])
        %     % TimeIndexing = (SaveCSVCounter-1)*SaveInterval + 1:SaveCSVCounter*SaveInterval;
        %     SaveCSVCounter = SaveCSVCounter + 1;
        %     % writematrix([TimeArray(TimeIndexing)';XArray],SaveStr);
        %     % writematrix([TimeArraySave;XArray],SaveStr);
        %     XArray = zeros(length(XTheta),SaveInterval);
        %     TimeArraySave  = zeros(1,SaveInterval);
        %     SaveCounter = 0;
        % end

        if mod(n-1,PlottingInterval) == 0
            PlotXandT
            drawnow
        end
        n = n + 1;

        % Save error and Average Phi

        % PhiAvg(n) = sum(FullVolumes .* X) / sum(NodeRho0.*FullVolumes);
        [Properties] = ComputeMaterialProperties(X,Theta,NodeRho0,Constants);
        epsG = Properties.epsG;
        RhoV = Properties.RhoV;
        PhiAvgRho0(n) = sum(NodeRho0 .* X .* FullVolumes) / sum(NodeRho0 .* FullVolumes);
        PhiAvgRho0T(n) = sum(NodeRho0 .* Theta .* Tbp .* FullVolumes) / sum(NodeRho0 .* FullVolumes);
        VapourAverage(n) = sum(epsG .* RhoV .* FullVolumes) / sum(NodeRho0 .* FullVolumes);
        SigmaSave(n) = FluxDiff;
        % if n >= 2
        %     DiffAvg(n) = sum(NodeRho0 .* FullVolumes)*((PhiAvgRho0(n) - PhiAvgRho0(n-1))...
        %         + (VapourAverage(n) - VapourAverage(n-1)))- FluxDiff * Tn;
        %     % fprintf('DiffAvg = %g\n',abs(DiffAvg(n)))
        % end
        % Step forward in time
        Tn = Tn * Alpha;
        Tn = min(Tn,TnMax);
        % Tn = min(TMax-t,Tn);

        %     frame = getframe(h);
        %     writeVideo(vidfile,frame);

    else
        Acounter = Acounter + 1;
        if ConsolePrint == 1
            fprintf(2,"Alpha heuristic not met, reducing timestep without advancing solution\n")
        end

        Tn = Alpha * Tn;
        continue
    end


end
%% save time results
SimTime = toc(Stime);
FinalData = [initTime,BTime,SimTime,Acounter]
n = n - 1;
PlotXandT
% if t + Tn > CurrentSaveTime
% SaveCounter = SaveCounter + 1;
% TimeArraySave(SaveCounter) = t;
% XArray(:,SaveCounter) = XTheta;
% % while CurrentSaveTime < t + Tn % make sure next save step is above current step
% %     CurrentSaveTime = CurrentSaveTime + SaveTimeInterval;
% % end
% % end
% % if SaveCounter == SaveInterval% Save data and stuff
% fprintf("Writting data to CSV File :::::: SimulationData_%g.csv\n",SaveCSVCounter)
% SaveStr = [DirStr, sprintf('/SimulationData_%g.csv',SaveCSVCounter)];
% save([DirStr,'/SimulationWorkspace.mat'])
% % TimeIndexing = (SaveCSVCounter-1)*SaveInterval + 1:SaveCSVCounter*SaveInterval;
% SaveCSVCounter = SaveCSVCounter + 1;
% % writematrix([TimeArray(TimeIndexing)';XArray],SaveStr);
% writematrix([TimeArraySave;XArray],SaveStr);
% XArray = zeros(length(XTheta),SaveInterval);
% SaveCounter = 0;
% end