function [EWPMesh] = OverlappingSurfacePoints(MeshStr,EWPTranslation,EWPRotation)
% we want to combine the meshes of multiple baords together, this includes
% rotating and translating them as needed
%
% Establish the convention that arrays are constructed in a top to bottom
% fashion. That is, array = [Top; Mid; Bot]!!!
% fprintf('Loading in mesh data\n')
% load mesh data
MeshStr1 = MeshStr{1};
NumBoards = numel(MeshStr);
if NumBoards == 1 %#ok
    load(MeshStr1,"MeshGeoPropsMid")

    MNodePos = MeshGeoPropsMid.NodePos;

    MElements = MeshGeoPropsMid.Elements;

    NumMidElements = length(MElements);

    NumElements = NumMidElements;

    ElementRotations = zeros(NumElements,3);
    ElementRotations(1:NumMidElements,:) = EWPRotation(1,:).*ones(NumMidElements,1);

    NumMidNodes = length(MeshGeoPropsMid.NodePos);
    NumNodes = NumMidNodes;

    NodeRotations = zeros(NumNodes,3);
    NodeRotations(1:NumMidNodes,:) = EWPRotation(1,:).*ones(NumMidNodes,1);

    Translations = zeros(NumNodes,3);
    Translations(1:NumMidNodes,:) = EWPTranslation(1,:).*ones(NumMidNodes,1);


    x = [MNodePos(:,1)];
    y = [MNodePos(:,2)];
    z = [MNodePos(:,3)];

    Elements = MElements;
    NumElements = length(Elements);

    NormalsX = [MeshGeoPropsMid.NormalVecsX];
    NormalsY = [MeshGeoPropsMid.NormalVecsY];
    NormalsZ = [MeshGeoPropsMid.NormalVecsZ];

    Triangles = [MeshGeoPropsMid.Triangles];
    Quads =  [MeshGeoPropsMid.Quads];

elseif NumBoards == 2

    MeshStr2 = MeshStr{2};
    load(MeshStr1,"MeshGeoPropsMid")
    load(MeshStr2,"MeshGeoPropsBot")

    MNodePos = MeshGeoPropsMid.NodePos;
    BNodePos = MeshGeoPropsBot.NodePos;
    MElements = MeshGeoPropsMid.Elements;
    BElements = MeshGeoPropsBot.Elements;
    NumMidElements = length(MElements);
    NumBotElements = length(BElements);
    NumElements = NumMidElements + NumBotElements;

    ElementRotations = zeros(NumElements,3);
    ElementRotations(1:NumMidElements,:) = EWPRotation(1,:).*ones(NumMidElements,1);
    ElementRotations(NumMidElements+1:end,:) = EWPRotation(2,:).*ones(NumBotElements,1);

    NumMidNodes = length(MeshGeoPropsMid.NodePos);
    NumBotNodes = length(MeshGeoPropsBot.NodePos);
    NumNodes = NumMidNodes + NumBotNodes;

    NodeRotations = zeros(NumNodes,3);
    NodeRotations(1:NumMidNodes,:) = EWPRotation(1,:).*ones(NumMidNodes,1);
    NodeRotations(NumMidNodes+1:end,:) = EWPRotation(2,:).*ones(NumBotNodes,1);

    Translations = zeros(NumNodes,3);
    Translations(1:NumMidNodes,:) = EWPTranslation(1,:).*ones(NumMidNodes,1);
    Translations(NumMidNodes+1:end,:) = EWPTranslation(2,:).*ones(NumBotNodes,1);


    x = [MNodePos(:,1); BNodePos(:,1)];
    y = [MNodePos(:,2); BNodePos(:,2)];
    z = [MNodePos(:,3); BNodePos(:,3)];

    BElements = MeshGeoPropsBot.Elements;
    Elements = [MElements; BElements + length(MNodePos)];
    NumElements = length(Elements);

    NormalsX = [MeshGeoPropsMid.NormalVecsX; MeshGeoPropsBot.NormalVecsX];
    NormalsY = [MeshGeoPropsMid.NormalVecsY; MeshGeoPropsBot.NormalVecsY];
    NormalsZ = [MeshGeoPropsMid.NormalVecsZ; MeshGeoPropsBot.NormalVecsZ];

    Triangles = [MeshGeoPropsMid.Triangles; MeshGeoPropsBot.Triangles + length(MNodePos)];
    Quads =  [MeshGeoPropsMid.Quads; MeshGeoPropsBot.Quads + length(MNodePos)];

end


[X,Y,Z] = EWPRotateAndTransform(NodeRotations,Translations,x,y,z);

Xmin = min(X); Xmax = max(X);
Ymin = min(Y); Ymax = max(Y);
Zmin = min(Z); Zmax = max(Z);

NodePos = zeros(NumNodes,3);
NodePos(:,1) = X;
NodePos(:,2) = Y;
NodePos(:,3) = Z;
% clf
% hold on

% Plot_Mesh3D(NodePos,Elements,([92 211 217])/255.*ones(NumElements,1),1,1)
% view([-45 25])
% % view([0 0])
% xlabel('X Direction [Metres]')
% ylabel('Y Direction [Metres]')
% zlabel('Z Direction [Metres]')
% 
% plot3(x(1),y(1),z(1),'k.','markersize',30);
% plot3(X(1),Y(1),Z(1),'r.','markersize',30);
% plot3(x(5312),y(5312),z(5312),'m.','markersize',30);
% plot3(X(5312),Y(5312),Z(5312),'y.','markersize',30);
% axis equal
%% Recompute mid mesh props

%

Xn = zeros(NumElements,9);
Yn = zeros(NumElements,9);
Zn = zeros(NumElements,9);
for i = 1:9



    x = NormalsX(:,i);
    y = NormalsY(:,i);
    z = NormalsZ(:,i);

    [Xn(:,i),Yn(:,i),Zn(:,i)] = EWPRotate(ElementRotations,x,y,z);
end

%% Reassign tri and quad element boundary Ids EXTERNAL ONLY

h = sqrt(eps);
TriSurfIndex = zeros(length(Triangles),1);
QuadSurfIndex = zeros(length(Quads),1);

TriX = X(Triangles); QuadX = X(Quads);
TriY = Y(Triangles); QuadY = Y(Quads);
TriZ = Z(Triangles); QuadZ = Z(Quads);
%%%%% X DIRECTION %%%%%
% X = X_min --- Idx = 1
idxT = all(abs(TriX - Xmin) < h,2); idxQ = all(abs(QuadX - Xmin)< h,2);
TriSurfIndex(idxT) = 1; QuadSurfIndex(idxQ) = 1;

% X = X_max --- Idx = 2
idxT = all(abs(TriX - Xmax)< h,2); idxQ = all(abs(QuadX - Xmax)< h,2);
TriSurfIndex(idxT) = 2; QuadSurfIndex(idxQ) = 2;

%%%%% Y DIRECTION %%%%%
% Y = Y_min --- Idx = 3
idxT = all(abs(TriY - Ymin)< h,2); idxQ = all(abs(QuadY - Ymin)< h,2);
TriSurfIndex(idxT) = 3; QuadSurfIndex(idxQ) = 3;

% Y = Y_max --- Idx = 4
idxT = all(abs(TriY - Ymax)< h,2); idxQ = all(abs(QuadY - Ymax)< h,2);
TriSurfIndex(idxT) = 4; QuadSurfIndex(idxQ) = 4;

%%%%% Z DIRECTION %%%%%
% Z = 0 --- Idx = 5
idxT = all(abs(TriZ - Zmin)< h,2); idxQ = all(abs(QuadZ - Zmin)< h,2);
TriSurfIndex(idxT) = 5; QuadSurfIndex(idxQ) = 5;

% Z = Z_max --- Idx = 6
idxT = all(abs(TriZ - Zmax)< h,2); idxQ = all(abs(QuadZ - Zmax)< h,2);
TriSurfIndex(idxT) = 6; QuadSurfIndex(idxQ) = 6;

% Glue Line % TEMP FIX
if NumBoards == 2
    idxQ = all((abs(QuadZ - 0.045) < h),2);
    QuadSurfIndex(idxQ) = 7;
end

%% Isolate the glue line quads
% We're gonna look at the pairwise distance between each quad for each
% direction seperatley. That is if meshes lie along the same plane thier
% pairwise distance will be 0.
GlueInteractions = zeros(NumNodes,3);
if NumBoards == 2
    NumMidQuads = length(MeshGeoPropsMid.Quads);

    XDist = squareform(pdist(QuadX)) < h;
    YDist = squareform(pdist(QuadY)) < h;
    ZDist = squareform(pdist(QuadZ)) < h;

    % w only need to look at the entries where the two meshes interact, that is
    % the bottom left section of the matrix
    ZDist(1:NumMidQuads,:) = 0;
    ZDist(:,NumMidQuads:end) = 0;

    % Find the find the indicies of the quads that are overlapping Zr is
    % bot->mid and Zc is mid -> bot
    [Zr,Zc] = find(ZDist);
    Zr = unique(Zr); Zc = unique(Zc);
    ZrQuads = Quads(unique(Zr),:); ZcQuads = Quads(unique(Zc),:);
    ZrNodes = unique(Quads(unique(Zr),:)); ZcNodes = unique(Quads(unique(Zc),:));
    %% We now need to find which nodes correspond with which quads in the other mesh
    % clf
    GlueInteractions = zeros(NumNodes,3);
    n = 1;

    ZrNodePos = NodePos(ZrNodes,:);
    ZcNodePos = NodePos(ZcNodes,:);
    Xr = ZrNodePos(:,1);
    Xc = ZcNodePos(:,1);
    Yr = ZrNodePos(:,2);
    Yc = ZcNodePos(:,2);
    Zr = ZrNodePos(:,3);
    Zc = ZcNodePos(:,3);

    MidNodePos = NodePos(1:NumMidNodes,:);
    BotNodePos = NodePos(NumMidNodes+1:end,:);

    % Loop over each quad in one of the meshes and see which points in the
    % other lie inside of it, if any
    for i = 1:length(ZcQuads)

        XQuad = NodePos(ZcQuads(i,:),1);
        YQuad = NodePos(ZcQuads(i,:),2);
        in = inpolygon(Xr,Yr,XQuad,YQuad); % logical index of points inside
        %%Plotting
        % plot([XQuad; XQuad],[YQuad; YQuad]) % polygon
        % hold on
        % plot(Xr(in),Yr(in),'r+') % points inside
        % plot(Xr(~in),Yr(~in),'bo') % points outside
        % plot(Xc,Yc,'k.')
        % hold off
        % grid on
        % grid minor
        % pbaspect([1 1 1])
        % drawnow

        if any(in)
            innodes = find(in);

            for node = 1:length(innodes) % Loop over each node found in the quad
                QNode = innodes(node);
                [z,n] = QuadShapeFuncLocalFromGlobal(Xr(QNode),Yr(QNode),XQuad,YQuad);
                % this is really bad code I'm sure there is a faster way to get
                % the global node Ids but this is easy. Simply search over the
                % global node list
                GlobalNodeIdxCol = ismember(BotNodePos,...
                    [Xr(QNode),Yr(QNode),Zr(QNode)],'rows');
                GlobalNodeIdx = find(GlobalNodeIdxCol)+NumMidNodes;
                QuadIdx = find(ismember(Quads,ZcQuads(i,:),'rows'));

                try
                    GlueInteractions(GlobalNodeIdx(1),:) = [QuadIdx,z,n];
                catch
                    1
                end
            end
        end

    end
    % Other Dir
    for i = 1:length(ZrQuads)

        XQuad = NodePos(ZrQuads(i,:),1);
        YQuad = NodePos(ZrQuads(i,:),2);
        in = inpolygon(Xc,Yc,XQuad,YQuad); % logical index of points inside
        %Plotting
        % % plot([XQuad; XQuad],[YQuad; YQuad]) % polygon
        % % hold on
        % % plot(Xc(in),Yc(in),'r+') % points inside
        % % plot(Xc(~in),Yc(~in),'bo') % points outside
        % % plot(Xr,Yr,'k.')
        % % hold off
        % % grid on
        % % grid minor
        % % pbaspect([1 1 1])
        % % drawnow
        % %
        if any(in)
            innodes = find(in);

            for node = 1:length(innodes) % Loop over each node found in the quad
                QNode = innodes(node);
                [z,n] = QuadShapeFuncLocalFromGlobal(Xc(QNode),Yc(QNode),XQuad,YQuad);
                % this is really bad code I'm sure there is a faster way to get
                % the global node Ids but this is easy. Simply search over the
                % global node list
                GlobalNodeIdx = ismember(MidNodePos,...
                    [Xc(QNode),Yc(QNode),Zc(QNode)],'rows');
                GlobalNodeIdx = find(GlobalNodeIdx);
                QuadIdx = find(ismember(Quads,ZrQuads(i,:),'rows'));

                try
                    GlueInteractions(GlobalNodeIdx(1),:) = [QuadIdx,z,n];
                catch
                    1
                end
            end
        end

    end
end

%%


% MNodePos = MNodePos(:,MidDimVec); % Apply X - Y roation

% BNodePos = BNodePos(:,[3 2 1]);

% Plot
% clf
% Plot_Mesh3D(NodePos,Elements,([92 211 217])/255.*ones(NumElements,1),1,1)
% view([-45 25])
% % view([0 0])
% xlabel('X Direction [Metres]')
% ylabel('Y Direction [Metres]')
% zlabel('Z Direction [Metres]')
% % pbaspect([1 1 1])



% fprintf('Isolating glue line interaction nodes\n')


%% Create new struct for big mesh
% fprintf('Creating EWP mesh struct\n')
EWPMesh = struct;
% NormalX = [MeshGeoPropsMid.NormalVecsX;
%             MeshGeoPropsBot.NormalVecsX];
% NormalZ = [MeshGeoPropsMid.NormalVecsZ;
%             MeshGeoPropsBot.NormalVecsZ];
MeshNames = fieldnames(MeshGeoPropsMid);
NumMidNodes = length(MNodePos);
NumBotNodes = length(MNodePos);
% Concatinate arrays, we have some special cases to consider
for k=1:numel(MeshNames)
    if MeshNames{k} == "Elements"
        EWPMesh.Elements = Elements;
    elseif MeshNames{k} == "NodePos"
        EWPMesh.NodePos = NodePos;
    elseif MeshNames{k} == "Triangles"
        EWPMesh.Triangles = Triangles;
    elseif MeshNames{k} == "Quads"
        EWPMesh.Quads = Quads;
    elseif MeshNames{k} == "NormalVecsX"
        EWPMesh.NormalVecsX = Xn;
    elseif MeshNames{k} == "NormalVecsY"
        EWPMesh.NormalVecsY = Yn;
    elseif MeshNames{k} == "NormalVecsZ"
        EWPMesh.NormalVecsZ = Zn;
    elseif MeshNames{k} == "QuadsSurfIndex"
        EWPMesh.QuadsSurfIndex = QuadSurfIndex;
    elseif MeshNames{k} == "TrianglesSurfIndex"
        EWPMesh.TrianglesSurfIndex = TriSurfIndex;
    elseif MeshNames{k} == "ScaleFacts"
        if NumBoards == 2
        EWPMesh.ScaleFacts = 0.5*(MeshGeoPropsMid.ScaleFacts ...
            + MeshGeoPropsBot.ScaleFacts);
        elseif NumBoards == 1 %#ok
            EWPMesh.ScaleFacts = MeshGeoPropsMid.ScaleFacts;
        end

    else % Just throw everything else together
        if NumBoards == 2
        EWPMesh.(MeshNames{k}) = [MeshGeoPropsMid.(MeshNames{k});
            MeshGeoPropsBot.(MeshNames{k})];
        elseif NumBoards == 1 %#ok
            EWPMesh.(MeshNames{k}) = [MeshGeoPropsMid.(MeshNames{k})];
        end
        
    end
end
EWPMesh.NodeRotations = NodeRotations;
EWPMesh.NodeTranslations = Translations;
EWPMesh.ElementRotations = ElementRotations;

EWPMesh.GlueInteractions = GlueInteractions;
EWPMesh.NumMidNodes = NumMidNodes;
EWPMesh.NumMidElements = NumMidElements;


% %% Plot
% % fprintf('Plotting\n')
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


%% Format for plotting interface
% 
% GlueQuadsIds = find(QuadSurfIndex == 7);
% MidQuads = Quads(GlueQuadsIds(GlueQuadsIds <= NumMidQuads),:);
% if NumBoards == 2
%     BotQuads = Quads(GlueQuadsIds(GlueQuadsIds >  NumMidQuads),:);
% end

%
% [NodeCol, CMin, CMax] = ReturnColsFromData(Nod, map);
%
% patch('Faces', PlottingQuads, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
%     'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);



end





function [z,n] = QuadShapeFuncLocalFromGlobal(xg,yg,xp,yp)
% This funtion is pretty slow but at least it works
% Perform a single nonlinear solve on the functions fx and fy where we feed
% in global coords of the quad verticies and the point we're interpolating
% too. It performs a single newton step as the functions are bilinear so it
% converges exactly within this step.

% Shape functions to return the global x and y coords when given local
% coords (z,n)
fx = @(z,n) 0.25*(xp(1)*(1-z).*(1-n) + xp(2)*(1+z).*(1-n) + xp(3)*(1+z).*(1+n) ...
    + xp(4)*(1-z).*(1+n)) - 1*xg;
fy = @(z,n) 0.25*(yp(1)*(1-z).*(1-n) + yp(2)*(1+z).*(1-n) + yp(3)*(1+z).*(1+n) ...
    + yp(4)*(1-z).*(1+n)) - 1*yg;

% derivative terms
Jxz = @(z,n) 0.25*(-xp(1)*(1-n) + xp(2)*(1-n) + xp(3)*(1+n) - xp(4)*(1+n));
Jyz = @(z,n) 0.25*(-yp(1)*(1-n) + yp(2)*(1-n) + yp(3)*(1+n) - yp(4)*(1+n));
Jxn = @(z,n) 0.25*(-xp(1)*(1-z) - xp(2)*(1+z) + xp(3)*(1+z) + xp(4)*(1-z));
Jyn = @(z,n) 0.25*(-yp(1)*(1-z) - yp(2)*(1+z) + yp(3)*(1+z) + yp(4)*(1-z));
% Intial guess of the position in the middle
z = 0;
n = 0;

% Setup newton step
b = -[fx(z,n);fy(z,n)]; % F(f(x))
J = [Jxz(z,n),Jxn(z,n);
    Jyz(z,n), Jyn(z,n)]; % jacobian
dx = J\b; % solve
z = dx(1); % Iterate solution (set equal as inital guess is 0 for z and n)
n = dx(2);




end