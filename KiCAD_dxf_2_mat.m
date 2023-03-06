% KiCAD_dxf_2_mat.m

clc
clear
close all

pkg load geometry
pkg load miscellaneous

load('simParams.mat');

if ~exist(Sim_Path,'dir')
   [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
end

%%%%%%%%%%%
% Predefs %
%%%%%%%%%%%

addpath(modelerPath);
addpath(miscPath);
addpath(openEMSpath);

%%%%%%%%%%%%%%%%%%%%%%%
% Parse and read DRLs %
%%%%%%%%%%%%%%%%%%%%%%%
drlFiles = dir([dxfPath '\*.drl']);

layerPairs = [];
drillD = [];
viaLoc = [];
for fIdx = 1:numel(drlFiles)
   [cLayerPairs,cD,cLoc,~,~] = parseDrillFile(fullfile(dxfPath,drlFiles(fIdx).name));
   layerPairs = [layerPairs ; cLayerPairs];
   drillD = [drillD ; cD];
   viaLoc = [viaLoc ; cLoc];
end

%%%%%%%%%%%%%%%%%%%%%%%
% Parse and read DXFs %
%%%%%%%%%%%%%%%%%%%%%%%

% Convert layers to height map
[elTable,subThickTable] = layerStruct_2_tables(layers);

%% Find all DXFs first
dxfFiles = dir([dxfPath '\*.dxf']);

layerPolygons = cell([numel(dxfFiles) 1]);
origPolygons = layerPolygons;
foundLayers = zeros(size(layerPolygons));

for fileIdx = 1:numel(dxfFiles)

   isLayer = @(s) strfind(dxfFiles(fileIdx).name,s.name);
   foundLayer = find(~cellfun(@isempty,arrayfun(isLayer,layers,'UniformOutput',false)));
   if isempty(foundLayer) foundLayer = 0; end;

   % Current DXF imported in mm
   [header,tables,entities] = mod2D_importDXF(fullfile(dxfPath,dxfFiles(fileIdx).name));

   % Convert the read DXF file to mod2D "objects"
   [dxfPolygons, ignoredEnts] = mod2D_convertEntityListToPolygons(entities);
   fprintf(1,'Total of %d entities ignored while converting to polygons\n',numel(ignoredEnts));

   % Subtract overlaps. This is to replace a manual macro.
   dxfPolygons = mod2D_subtractOvelapingPolygons(dxfPolygons);

   % Store polygons, and which layer was found
   layerPolygons{fileIdx} = dxfPolygons;
   foundLayers(fileIdx) = foundLayer;
end

%%%%%%%%%%%%%%%%%%%
% Feature removal %
%%%%%%%%%%%%%%%%%%%

for polsIdx = 1:numel(layerPolygons)
   cPolygons = layerPolygons{polsIdx};
   for polIdx = 1:numel(cPolygons)
      cPoly = cPolygons{polIdx};

      cPoly = removeSmallFeatures(cPoly,minPointRes);

      % If exist, remove any duplicate points
      % [uNode,uIdx,u2oIdx] = mod2D_uniqueNodeByDist(node,distTH)
      [uNodes,~,~] = mod2D_uniqueNodeByDist([cPoly.x cPoly.y],samePointRes);
      if size(uNodes,1) < numel(cPoly.x)
         cPoly.x = uNodes(:,1);
         cPoly.y = uNodes(:,2);
      end

      cPolygons{polIdx} = cPoly;

   end

   layerPolygons{polsIdx} = cPolygons;



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split for APPCSXCAD readability %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

origPolygons = layerPolygons;

for polsIdx = 1:numel(layerPolygons)
   cPolygons = layerPolygons{polsIdx};

   % Before splitting, store all coordinates, for later use with mesh
   cPolygons = splitSelfIntersectingBodies(cPolygons);

   layerPolygons{polsIdx} = cPolygons;

end

%%%%%%%%%%%%%%
% Store data %
%%%%%%%%%%%%%%

save(...
   fullfile(Sim_Path,'model_data.mat'),...
   '-mat7-binary',...
   'layerPolygons',...
   'origPolygons',...
   'foundLayers',...
   'layers',...
   'materials',...
   'metals',...
   'layerPairs',...
   'drillD',...
   'viaLoc');


rmpath(modelerPath);
rmpath(miscPath);
rmpath(openEMSpath);
