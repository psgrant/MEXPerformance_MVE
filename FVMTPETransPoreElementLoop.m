function [g] = FVMTPETransPoreElementLoop(NumNodes, NumElements, Elements, NormalVectorsX, NormalVectorsY, NormalVectorsZ, ElementRho0, ShapeFuncArrayT, DbX, DbY, DbZ, DbXY, DbXZ, DbYZ, DeffX, DeffY, DeffZ, DeffXY, DeffXZ, DeffYZ, KeffX, KeffY, KeffZ, KeffXY, KeffXZ, KeffYZ, KX, KY, KZ, KXY, KXZ, KYZ, kwX, kwY, kwZ, kwXY, kwXZ, kwYZ, muw, RhoG, wv, ha, hb, hv, hw, PressureHead, RhoWGrav, RhoW, GradPwX, GradPwY, GradPwZ, GradXbX, GradXbY, GradXbZ, GradWvX, GradWvY, GradWvZ, GradTX, GradTY, GradTZ, SubAreas)%#codegen


NodesFrom = [1 2 3 4 5 6 1 2 3];
NodesTo = [2 3 1 5 6 4 4 5 6];
% Loop over each element
g = zeros(NumNodes*2,1);
for ele = 1:NumElements

 eNodes = Elements(ele,:);

 nx = NormalVectorsX(ele,:);
 ny = NormalVectorsY(ele,:);
 nz = NormalVectorsZ(ele,:);

 eRho0 = ElementRho0(ele);

 eDbX = ShapeFuncArrayT * DbX(:,ele);
 eDbY = ShapeFuncArrayT * DbY(:,ele); % update all with the matrix multiply
 eDbZ = ShapeFuncArrayT * DbZ(:, ele);
 eDbXY = ShapeFuncArrayT * DbXY(:, ele);
 eDbXZ = ShapeFuncArrayT * DbXZ(:, ele);
 eDbYZ = ShapeFuncArrayT * DbYZ(:, ele);




 eDeffX = ShapeFuncArrayT * DeffX(:, ele);
 eDeffY = ShapeFuncArrayT * DeffY(:, ele);
 eDeffZ = ShapeFuncArrayT * DeffZ(:, ele);
 eDeffXY = ShapeFuncArrayT * DeffXY(:, ele);
 eDeffXZ = ShapeFuncArrayT * DeffXZ(:, ele);
 eDeffYZ = ShapeFuncArrayT * DeffYZ(:, ele);

 eKeffX = ShapeFuncArrayT * KeffX(:, ele);
 eKeffY = ShapeFuncArrayT * KeffY(:, ele);
 eKeffZ = ShapeFuncArrayT * KeffZ(:, ele);
 eKeffXY = ShapeFuncArrayT * KeffXY(:, ele);
 eKeffXZ = ShapeFuncArrayT * KeffXZ(:, ele);
 eKeffYZ = ShapeFuncArrayT * KeffYZ(:, ele);


 % Upwinding total head water pressure + gravity
 eKX = ShapeFuncArrayT * KX(:, ele);
 eKY = ShapeFuncArrayT * KY(:, ele);
 eKZ = ShapeFuncArrayT * KZ(:, ele);
 eKXY = ShapeFuncArrayT * KXY(:, ele);
 eKXZ = ShapeFuncArrayT * KXZ(:, ele);
 eKYZ = ShapeFuncArrayT * KYZ(:, ele);

 ekwX = ShapeFuncArrayT * kwX(:, ele);
 ekwY = ShapeFuncArrayT * kwY(:, ele);
 ekwZ = ShapeFuncArrayT * kwZ(:, ele);
 ekwXY = ShapeFuncArrayT * kwXY(:, ele);
 ekwXZ = ShapeFuncArrayT * kwXZ(:, ele);
 ekwYZ = ShapeFuncArrayT * kwYZ(:, ele);



 emuw = ShapeFuncArrayT * muw(:, ele);

 eRhoG = ShapeFuncArrayT * RhoG(:, ele);

 eWv = ShapeFuncArrayT * wv(:, ele);

 eha = ShapeFuncArrayT * ha(:, ele);
 ehb = ShapeFuncArrayT * hb(:, ele);
 ehv = ShapeFuncArrayT * hv(:, ele);
 ehw = ShapeFuncArrayT * hw(:, ele);

 ePressureHead = PressureHead(eNodes);


 % Loop over the 9 flux-passing faces
 for ip = 1:9

  MyNode = NodesFrom(ip);
  OtherNode = NodesTo(ip);

  Node1 = eNodes(MyNode);
  Node2 = eNodes(OtherNode);

  % Upwinding pressure stuffs
  % if pressure head at my node is greater than the other node we
  % take kw to our value
  if ePressureHead(MyNode) > PressureHead(OtherNode)
   ekwX_UP = ekwX(MyNode);
   ekwY_UP = ekwY(MyNode);
   ekwZ_UP = ekwZ(MyNode);

   ekwXY_UP = ekwXY(MyNode);
   ekwXZ_UP = ekwXZ(MyNode);
   ekwYZ_UP = ekwYZ(MyNode);

   emuw_UP = emuw(MyNode);

   ehw_UP = ehw(MyNode);
  else
   ekwX_UP = ekwX(OtherNode);
   ekwY_UP = ekwY(OtherNode);
   ekwZ_UP = ekwZ(OtherNode);

   ekwXY_UP = ekwXY(OtherNode);
   ekwXZ_UP = ekwXZ(OtherNode);
   ekwYZ_UP = ekwYZ(OtherNode);

   emuw_UP = emuw(OtherNode);

   ehw_UP = ehw(OtherNode);
  end

  % LiquidAdvection
  LiquidAdvection = RhoWGrav * RhoW / emuw_UP * [  eKX(ip) * ekwX_UP * GradPwX(ele,ip)  +  eKXY(ip) * ekwXY_UP * GradPwY(ele,ip)  +  eKXZ(ip) * ekwXZ_UP * GradPwZ(ele,ip);
   eKXY(ip) * ekwXY_UP * GradPwX(ele,ip)  + eKY(ip) * ekwY_UP * GradPwY(ele,ip)  +  eKYZ(ip) * ekwYZ_UP * GradPwZ(ele,ip)
   eKXZ(ip) * ekwXZ_UP * GradPwX(ele,ip)  +  eKYZ(ip) * ekwYZ_UP * GradPwY(ele,ip)  +   eKZ(ip) * ekwZ_UP * GradPwZ(ele,ip)];

  % Bound liquid (only when Xb < XFSP)
  BoundLiquidFlux = eRho0 * [eDbX(ip) * GradXbX(ele,ip)   +  eDbXY(ip) * GradXbY(ele,ip)  +  eDbXZ(ip) * GradXbZ(ele,ip);
   eDbXY(ip) * GradXbX(ele,ip)  +   eDbY(ip) * GradXbY(ele,ip)  +  eDbYZ(ip) * GradXbZ(ele,ip);
   eDbXZ(ip) * GradXbX(ele,ip)  +  eDbYZ(ip) * GradXbY(ele,ip)  +   eDbZ(ip) * GradXbZ(ele,ip)];

  % Vapour Diffusion
  RhoGDeffGradWv = eRhoG(ip) *[ eDeffX(ip) * GradWvX(ele,ip)  +  eDeffXY(ip) * GradWvY(ele,ip)  +  eDeffXZ(ip) * GradWvZ(ele,ip);
   eDeffXY(ip) * GradWvX(ele,ip)  +   eDeffY(ip) * GradWvY(ele,ip)  +  eDeffYZ(ip) * GradWvZ(ele,ip);
   eDeffXZ(ip) * GradWvX(ele,ip)  +  eDeffYZ(ip) * GradWvY(ele,ip)  +   eDeffZ(ip) * GradWvZ(ele,ip)];


  VapourDiffusion =  RhoGDeffGradWv ./ (1 - eWv(ip));

  LiquidAdvectionEnthalpy = ehw_UP * LiquidAdvection;


  AirDiffusionEnthalpy = -eha(ip) * RhoGDeffGradWv ./ (eWv(ip));

  % Bound water diffusion Enthalpy
  BoundLiquidEnthalpyFlux = ehb(ip) * BoundLiquidFlux;

  % Vapour Diffusion Enthalpy
  VapourDiffusionEnthalpy = ehv(ip) * VapourDiffusion;

  % Temperature Conduction
  TempConduction =  1*[ eKeffX(ip) * GradTX(ele,ip)  +  eKeffXY(ip) * GradTY(ele,ip)  +  eKeffXZ(ip) * GradTZ(ele,ip);
   eKeffXY(ip) * GradTX(ele,ip)  +   eKeffY(ip) * GradTY(ele,ip)  +  eKeffYZ(ip) * GradTZ(ele,ip);
   eKeffXZ(ip) * GradTX(ele,ip)  +  eKeffYZ(ip) * GradTY(ele,ip)  +   eKeffZ(ip) * GradTZ(ele,ip)];



  % Create water (liquid + vapour) flux vector
  q_w = LiquidAdvection + BoundLiquidFlux + VapourDiffusion;

  % Create Theta (scaled temp) flux vector
  q_e  =  LiquidAdvectionEnthalpy + BoundLiquidEnthalpyFlux...
   + VapourDiffusionEnthalpy + AirDiffusionEnthalpy + TempConduction;

  % q_e = TempConduction;
  % Compute fluxes
  wFlux = SubAreas(ele,ip) * dot([nx(ip); ny(ip); nz(ip)] , q_w);
  eFlux = SubAreas(ele,ip) * dot([nx(ip); ny(ip); nz(ip)] , q_e);

  % Add to flux vector (remeber the negative flux for the other node)

  g(Node1*2-1) = g(Node1*2-1) + wFlux;
  g(Node2*2-1) = g(Node2*2-1) - wFlux;


  g(Node1*2) = g(Node1*2) + eFlux;
  g(Node2*2) = g(Node2*2) - eFlux;

 end


end

end


