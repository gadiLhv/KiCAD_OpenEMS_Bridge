clc
clear
close all

load('simParams.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize openEMS environment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_CSX = 'openEMS_sim.xml';

load(fullfile(Sim_Path,'model_data.mat'));
load(fullfile(Sim_Path,'model_ports.mat'));

pkg load geometry
pkg load miscellaneous

addpath(modelerPath);
addpath(miscPath);
addpath(openEMSpath);

physical_constants;

CSX = InitCSX();

% Setup initial simulation conditions
if exist('timeStep','var')
   FDTD = InitFDTD('NrTS',  maxTimeSteps,'TimeStep',timeStep,'EndCriteria',10^(endSimulationEnergy_dB/10));
else
   FDTD = InitFDTD('NrTS',  maxTimeSteps,'EndCriteria',10^(endSimulationEnergy_dB/10));
end


FDTD = SetGaussExcite( FDTD, f_center, f_cutoff );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate radiation addition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Currently add wl/2
fmin = f_center - f_cutoff;
wl_max = c0/fmin;

bBoxAdd = wl_max*airbox_pad_wl/unit;

%%%%%%%%%%%%%%%%%%%%
% Define materials %
%%%%%%%%%%%%%%%%%%%%

conductorNames = {};
for metalIdx = 1:numel(metals)
   cMetal = metals(metalIdx);

   % Find thickness of these metals in the layers list
   isThisMetal = @(s) strcmp(s.material,cMetal.name);
   metalFound = arrayfun(isThisMetal,layers);


   if ~max(metalFound)
      fprintf('\nNo metal named ''%s'' found in layers. Assuming PEC\n',cMetal.name);
      CSX = AddMetal(CSX,cMetal.name);
      conductorNames(numel(conductorNames) + 1) = {cMetal.name};
      continue;
   end

   % Transform thickness to [m]
   cThick = layers(find(metalFound,1)).thickness*1e-3;
   cSigma = cMetal.sigma;


   metalFoundIdxs = find(metalFound);
   for mIdx = metalFoundIdxs(:).'
      cThick = layers(mIdx).thickness*unit;
      cSigma = cMetal.sigma;
      conductorNames(numel(conductorNames) + 1) = {[cMetal.name '-' layers(mIdx).name]};
      CSX = AddConductingSheet(CSX,conductorNames{end},cSigma,cThick);
   end

end

% Add dielectrics
for matIdx = 1:numel(materials)

   cMaterial = materials(matIdx);

   % Calculate sigma from loss tangent
   epp = cMaterial.epsR*cMaterial.tanD;

   % Estimate conductivity at center frequency.
   cKappa = epp*2*pi*f_center*EPS0;

   CSX = AddMaterial(CSX,cMaterial.name);
   CSX = SetMaterialProperty( CSX, cMaterial.name, 'Epsilon',cMaterial.epsR,'Kappa',cKappa);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load layers into simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detect all metal layers in simulation
layerIsMetal = @(s) ~isempty(strfind(s.type,'conductor'));
%arrayfun(layerIsMetal,layers,'UniformOutput',false)
binMetalLayers = arrayfun(layerIsMetal,layers);
metalLayers = find(binMetalLayers);

layersNoThick = layers;

for lIdx = foundLayers(~~foundLayers).';
   layersNoThick(lIdx).thickness = 0;
end
[elTable,subThickTable] = layerStruct_2_tables(layersNoThick);


for lIdx = 1:numel(elTable)
   isMetal = binMetalLayers(lIdx);

   foundLayer = find(foundLayers(:) == lIdx);

   % If no layer was found, reference to 0 layer (edge cuts)
   if isempty(foundLayer)
      if isMetal
         continue;
      end
      foundLayer = find(foundLayers(:) == 0);
   end

   cPolygons = layerPolygons{foundLayer};

   % Material name per layer
   if isMetal
      matName = [layers(lIdx).material '-' layers(lIdx).name];
   else
      matName = layers(lIdx).material;
   end
   % Load all polygons into workspace
   for pIdx = 1:numel(cPolygons)
      cPoly = cPolygons{pIdx};
      p = [cPoly.x cPoly.y].';
      CSX = AddLinPoly(...
            CSX,...                             % Handle for the openEMS environment
            matName, ...                        % Material name
            isMetal*10 + (~isMetal)*1, ...      % Priority. 10 For metal, 2 for the rest
            'z', ...                            % Normal direciton
            elTable(lIdx), ...                  % Elevation - TBD
            p, ...                              % Elevation (To be fixed)
            subThickTable(lIdx)*(~isMetal));    % Currently forcing zero thickness on all metalic layers
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vias into simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reverse order, because KiCAD considers top layer #1
metalLayers = flipud(metalLayers);

for viaIdx = 1:size(viaLoc,1)

   cLoc = viaLoc(viaIdx,:);
   cD = drillD(viaIdx);
   cPair = layerPairs(viaIdx,:);

   stopEl = elTable(metalLayers(cPair(1)));
   startEl = elTable(metalLayers(cPair(2))) + subThickTable(metalLayers(cPair(2)));

   CSX = AddBox(CSX,...
            viaMaterial,...
            10,...
            [cLoc 0] + [-cD/sqrt(2) -cD/sqrt(2) startEl],...
            [cLoc 0] + [cD/sqrt(2) cD/sqrt(2) stopEl]);

%   CSX = AddCylinder(...
%            CSX,...
%            viaMaterial,...      % Vias are always metalic
%            1,...                % Metals are always top priority
%            [cLoc startEl],...   % Start coordinates
%            [cLoc stopEl],...    % End coordinates
%            cD*0.5);             % Radius
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add ports: Currently only lumped ports %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

portPoints = [];
portObjs = {};
for pIdx = 1:numel(ports)
   cPort = ports(pIdx);

   mid1 = mean(cPort.edge1);
   mid2 = mean(cPort.edge2);
   dl = (mid1 - mid2)/sqrt(sum((mid1 - mid2).^2,2));

   [CSX newPort] = AddLumpedPort(...
      CSX, ...
      15, ...           % Ports are priority 15
      pIdx, ...         % Port number
      portZ0, ...       % Impedance
      cPort.edge1(1,:), ...         % Start point
      cPort.edge2(2,:), ...         % Stop point
      dl, ...           % Direction
      true);            % To excite?


   portPoints = [portPoints ; cPort.edge1 ; cPort.edge2];

   portObjs{numel(portObjs) + 1} = newPort;
end



%%%%%%%%%%%
% Mesh!!! %
%%%%%%%%%%%

allMesh = [];
for metalIdx = 1:numel(conductorNames);
   cMesh = DetectEdges(CSX, [], 'SetProperty',conductorNames{metalIdx});
   if isempty(allMesh)
      allMesh = cMesh;
   else
      allMesh.x = [allMesh.x cMesh.x];
      allMesh.y = [allMesh.y cMesh.y];
      allMesh.z = [allMesh.z cMesh.z];
   end
end

for materialIdx = 1:numel(materials);
   cMesh = DetectEdges(CSX, [], 'SetProperty',materials(materialIdx).name);
   allMesh.x = [allMesh.x cMesh.x];
   allMesh.y = [allMesh.y cMesh.y];
   allMesh.z = [allMesh.z cMesh.z];
end

minMax = @(x) [min(x(:)) max(x(:))];
bBox = [minMax(allMesh.x) minMax(allMesh.y) minMax(allMesh.z)];

% Define original mesh
allMesh.x = unique([allMesh.x portPoints(:,1).']);
allMesh.y = unique([allMesh.y portPoints(:,2).']);
allMesh.z = unique([allMesh.z portPoints(:,3).']);

% These points are considered "immovable"
hardMesh = allMesh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate largest mesh step allowed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmax = f_center + f_cutoff;

maxDk = 0;
for mIdx = 1:numel(materials)
   maxDk =  (maxDk < materials(mIdx).epsR)*materials(mIdx).epsR + ...
            (maxDk >= materials(mIdx).epsR)*maxDk;
end

step_max = (1/unit)*maxStep_wl*(c0/fmax)/sqrt(maxDk);

%%%%%%%%%%%%%%%%%%%%%%%
% Generate sub-meshes %
%%%%%%%%%%%%%%%%%%%%%%%

allMesh.x = AutoSmoothMeshLines(allMesh.x,step_max(1), gradRatio,'symmetric',0,'homogeneous',0,'algorithm',[ 3]);
allMesh.y = AutoSmoothMeshLines(allMesh.y,step_max(2), gradRatio,'symmetric',0,'homogeneous',0,'algorithm',[ 3]);
allMesh.z = AutoSmoothMeshLines(allMesh.z,step_max(3), gradRatio,'symmetric',0,'homogeneous',0,'algorithm',[ 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth mesh in areas where step size is smaller than minimum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine mesh limitations in material
minMeshCell = [min(diff(allMesh.x(:))) min(diff(allMesh.y(:))) min(diff(allMesh.z(:)))];
dt_smallest = minMeshCell/(c0/sqrt(maxDk));

% Stupid approximation
dt_min_approx = dt_smallest/2e4;

axNames = 'xyz';
for axIdx = 1:3
   cMin = minCellSize_abs(axIdx);

   if ~cMin
      continue;
   end

   cAxVals = getfield(allMesh,axNames(axIdx));
   cAxMesh = [];



   cDiff = diff(cAxVals(:));

   while max(cDiff < cMin)
      % Find first small mesh cell
      i1 = find(cDiff < cMin,1);
      i2 = i1;
      while cDiff(i2) < cMin
         i2 = i2 + 1;
      end

      % Take one cell forward and one cell backwards
      meshSnip = [cAxVals(i1-1) cAxVals(i2+1)];

      minNparts = floor((meshSnip(2) - meshSnip(1))/cMin);
      smoothSnip = linspace(meshSnip(1),meshSnip(2),minNparts + 1);

      cAxVals = [cAxVals(1:(i1 - 2)) smoothSnip cAxVals((i2 + 2):end)];
      cDiff = diff(cAxVals(:));
   end

   eval(['allMesh.' axNames(axIdx) ' = cAxVals;']);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Add air box and mesh %
%%%%%%%%%%%%%%%%%%%%%%%%
step_max_air = (1/unit)*airboxMaxStep_wl*(c0/fmax);

meshLeft = AutoSmoothMeshLines([(allMesh.x(1) - [bBoxAdd 0]) allMesh.x(2)],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshRight = AutoSmoothMeshLines([allMesh.x(end-1) (allMesh.x(end) + [0 bBoxAdd])],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshDown = AutoSmoothMeshLines([(allMesh.y(1) - [bBoxAdd 0]) allMesh.y(2)],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshUp = AutoSmoothMeshLines([allMesh.y(end-1) (allMesh.y(end) + [0 bBoxAdd])],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshBot = AutoSmoothMeshLines([(allMesh.z(1) - [bBoxAdd 0]) allMesh.z(2)],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);
meshTop = AutoSmoothMeshLines([allMesh.z(end-1) (allMesh.z(end) + [0 bBoxAdd])],step_max_air(1), gradRatio,'symmetric',0,'homogeneous',0);

allMesh.x = [meshLeft(1:(end-2)) allMesh.x meshRight(3:end)];
allMesh.y = [meshDown(1:(end-2)) allMesh.y meshUp(3:end)];
allMesh.z = [meshBot(1:(end-2)) allMesh.z meshTop(3:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine mesh limitations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CSX = DefineRectGrid( CSX, unit, allMesh );

% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
start = [allMesh.x(4)     allMesh.y(4)     allMesh.z(4)];
stop  = [allMesh.x(end-3) allMesh.y(end-3) allMesh.z(end-3)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX,'--dump-statistics');

save(...
   fullfile(Sim_Path,'sim_data.mat'),...
   '-mat7-binary',...
   'f_center',...
   'f_cutoff',...
   'maxTimeSteps',...
   'portZ0',...
   'unit',...
   'portObjs');
