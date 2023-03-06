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
load(fullfile(Sim_Path,'openEMS_data.mat'));

pkg load geometry
pkg load miscellaneous

addpath(modelerPath);
addpath(miscPath);
addpath(openEMSpath);

physical_constants;

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

step_min = (1/unit)*minStep_wl*(c0/fmax)/sqrt(maxDk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force sub minimal count %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

axVect = 'xyz';
for axIdx = 1:3
   cMesh = getfield(allMesh,axVect(axIdx));

   cDiff = diff(cMesh);
   binBadDiv = ((ceil(cDiff/step_max(axIdx)) - 1) < minEdge2EdgeCount_abs(axIdx));
   badDivIdxs = find(binBadDiv);

   startIdx = 1;
   newMesh = [];
   for dIdx = badDivIdxs

      newDiv = linspace(cMesh(dIdx),cMesh(dIdx + 1),minEdge2EdgeCount_abs(axIdx) + 2);

      newMesh = [newMesh cMesh(startIdx:dIdx)];
      % Also, verify absolute minimum step restraint
      % Bad bypass. Need to re-calculate the difference
%      if (newDiv(2) - newDiv(1)) >= step_min(axIdx)
         newMesh = [newMesh newDiv(2:(end-1))];
%      end
      startIdx = dIdx + 1;
   end

   newMesh = [newMesh cMesh(startIdx:end)];

%   figure;
%   plot(newMesh(:),ones(size(newMesh(:))),'.r','markersize',20);
%   hold on;
%   plot(cMesh(:),ones(size(cMesh(:))),'.b','markersize',10);
%   hold off;

   % Stock back after done
   eval(['allMesh.' axVect(axIdx) ' = newMesh;']);

end


%%%%%%%%%%%%%%%%%%%%%%%
% Generate sub-meshes %
%%%%%%%%%%%%%%%%%%%%%%%

allMesh.x = AutoSmoothMeshLines(allMesh.x,step_max(1), gradRatio,'symmetric',0,'homogeneous',0,'algorithm',[ 3]);
allMesh.y = AutoSmoothMeshLines(allMesh.y,step_max(2), gradRatio,'symmetric',0,'homogeneous',0,'algorithm',[ 3]);
allMesh.z = AutoSmoothMeshLines(allMesh.z,step_max(3), gradRatio,'symmetric',0,'homogeneous',0,'algorithm',[ 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove duplicate lines %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here very close lines will be removed. This can happen for various reasons,
% but mostly slightly crooked lines. This is a common artifact from KiCAD
allMesh.x = removeSameMeshPoints(allMesh.x,samePointRes);
allMesh.y = removeSameMeshPoints(allMesh.y,samePointRes);
allMesh.z = removeSameMeshPoints(allMesh.z,samePointRes);

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

save( fullfile(Sim_Path,'mesh_data.mat'),...
      'allMesh',...
      'hardMesh',...
      'step_max',...
      'maxDk',...
      'fmax');


