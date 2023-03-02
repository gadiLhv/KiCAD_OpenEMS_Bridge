% KiCAD_mat_define_port_faces.m

clc
clear
close all

pkg load geometry
pkg load miscellaneous

modelerPath = 'C:\Repositories\EM2D_Solver\Modeler_2D';
miscPath = 'C:\Repositories\EM2D_Solver\misc';

Sim_Path = 'Corrected_PIFA_1p6mm';

load('simParams.mat');

%%%%%%%%%%%
% Predefs %
%%%%%%%%%%%

addpath(modelerPath);
addpath(miscPath);

load(fullfile(Sim_Path,'model_data.mat'));

%%%%%%%%%%%%%%
% Start menu %
%%%%%%%%%%%%%%

% Detect all metal layers in simulation
layerIsMetal = @(s) ~isempty(strfind(s.type,'conductor'));

metalLayers = find(arrayfun(layerIsMetal,layers));

layersNoThick = layers;

for lIdx = foundLayers(~~foundLayers).';
   layersNoThick(lIdx).thickness = 0;
end
[elTable,subThickTable] = layerStruct_2_tables(layersNoThick);

% Container for port edges, including elevation
portEdges = [];

while true
   %%%%%%%%%%%%%%%%%%%%
   % Promt layer menu %
   %%%%%%%%%%%%%%%%%%%%

   fprintf('Choose layer to place ports on:\n0) Exit\n');
   for lIdx = 1:numel(metalLayers)
      fprintf('%d) %s\n',lIdx,layers(metalLayers(lIdx)).name);
   end
   editLayer = input('Choice: ');

   if (editLayer > 0) && (editLayer <= numel(metalLayers))
      fprintf('\nLoading layer #%d, ''%s''...\n',editLayer,layers(metalLayers(editLayer)).name);
   elseif editLayer == 0
      break;
   else
      fprintf('Wrong choice. Try again!');
   end

   % Pull out correct set of polygons
   foundLayer = find(foundLayers(:) == metalLayers(editLayer));
   cPolygons = origPolygons{foundLayer};

   %%%%%%%%%%%%%%%
   % Start "GUI" %
   %%%%%%%%%%%%%%%
   figHdl = figure;
   allX = [];
   allY = [];
   for pIdx = 1:numel(cPolygons)
      cPolygon = cPolygons{pIdx};
      allX = [allX ; cPolygon.x];
      allY = [allY ; cPolygon.y];
      mod2D_showPolygon(gca,cPolygon,[hex2dec('B8') hex2dec('73') hex2dec('33')]/255,[0 0 0]);
   end
   hold off;

   % Cleare calculate bounding box and centroid
   allX(isnan(allX)) = [];
   allY(isnan(allY)) = [];
   bBox = [min(allX(:)) max(allX(:)) min(allY(:)) max(allY(:))];
   bBoxSize = [(bBox(2) - bBox(1)) (bBox(4) - bBox(3))];
   cent = [mean(bBox(1:2).') mean(bBox(3:4).')];

   % Make nice axes
   axis([cent(1) + 0.5*bBoxSize(1)*[-1.1 1.1] cent(2) + 0.5*bBoxSize(2)*[-1.1 1.1]]);
   axis('square');

   %%%%%%%%%%%%%%%%%%%%%%%%
   % Wait for first click %
   %%%%%%%%%%%%%%%%%%%%%%%%

   % Initialize edge layers as same layer
   e1Layer = editLayer;
   e2Layer = editLayer;

   while true

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Wait for second click\keyboard press %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      title('Click first edge on this layer');

      [x1,y1,in1] = ginput(1);
      if in1 <= 5
         % Now look for closest edge
         minDistPoly = 0;
         minDist = inf;
         minDistEdge = [0 0 ; 0 0 ];
         for pIdx = 1:numel(cPolygons)
            cPolygon = cPolygons{pIdx};

            v1 = [cPolygon.x(1:(end-1)) cPolygon.y(1:(end-1))];
            v2 = [cPolygon.x(2:end) cPolygon.y(2:end)];

            [ds,ts] = dist2edge(v2,v1,[x1 y1]);
            ds((ts < 0) | (ts > 1)) = inf;
            [cMinDist,cMinIdx] = min(ds.');
            if (cMinDist < minDist)
               minDist = cMinDist;
               minDistPoly = pIdx;
               minDistEdge = [v1(cMinIdx,:) ; v2(cMinIdx,:)];
            end
         end
         hold on;
         plot(minDistEdge(:,1),minDistEdge(:,2),'-r','linewidth',4);
         txtHdl = text(x1,y1,'Edge #1','fontsize',16);
         hold off;

         portEdge1 = minDistEdge;

         break;
      end
   end

   while true
      % Check if upper and lower layers are available
      msgCells = {};
      if editLayer < numel(metalLayers)
         msgCells(1) = sprintf('For layer ''%s'', press ''u''',layers(metalLayers(editLayer + 1)).name);
      end
      if editLayer > 1
         msgCells(numel(msgCells) + 1) = sprintf('For layer ''%s'', press ''l''',layers(metalLayers(editLayer - 1)).name);
      end
      msgCells(numel(msgCells) + 1) = 'For another edge on same layer, click mouse';

      title(msgCells);

      [x2,y2,in2] = ginput(1);
      if in2 <= 5
         % Now look for closest edge
         minDistPoly = 0;
         minDist = inf;
         minDistEdge = [0 0 ; 0 0 ];
         for pIdx = 1:numel(cPolygons)
            cPolygon = cPolygons{pIdx};

            % If this is a self intersecting polygon, don't close the last loop
            if max(isnan(cPolygon.x))
               v1 = [cPolygon.x(1:(end-1)) cPolygon.y(1:(end-1))];
               v2 = [cPolygon.x(2:end) cPolygon.y(2:end)];
            else
               v1 = [cPolygon.x cPolygon.y];
               v2 = [[cPolygon.x(2:end) ; cPolygon.x(1)] [cPolygon.y(2:end) ; cPolygon.y(1)]];
            end
            [ds,ts] = dist2edge(v2,v1,[x2 y2]);
            ds((ts < 0) | (ts > 1)) = inf;
            [cMinDist,cMinIdx] = min(ds.');
            if (cMinDist < minDist)
               minDist = cMinDist;
               minDistPoly = pIdx;
               minDistEdge = [v1(cMinIdx,:) ; v2(cMinIdx,:)];
            end
         end
         hold on;
         plot(minDistEdge(:,1),minDistEdge(:,2),'-r','linewidth',4);
         text(x2,y2,'Edge #2','fontsize',16);
         hold off;

         portEdge2 = minDistEdge;

         % Create smallest possible rectangle:
         [~,ts] = dist2edge(portEdge1(2,:),portEdge1(1,:),portEdge2);

         [ts,sortIdxs] = sort(ts);
         % Order them according to 'ts'
         portEdge2 = portEdge2(sortIdxs,:);

         oldEdge1 = portEdge1;
         oldEdge2 = portEdge2;
         % Project to see which one is larger
         dl = oldEdge1(2,:) - oldEdge1(1,:);
         tau = dl/sqrt(sum(dl(:).^2));
         if ts(1) < 0
            portEdge2(1,:) = oldEdge2(1,:) - ts(1)*dl;
         else
            portEdge1(1,:) = oldEdge1(1,:) + ts(1)*dl;
         end

         % Same with 2nd point
         if ts(2) > 1
            portEdge2(2,:) = oldEdge2(2,:) - (ts(2) - 1)*dl;
         else
            portEdge1(2,:) = oldEdge1(2,:) + (ts(2) - 1)*dl;
         end

         % Show user the created port
         portPol = mod2D_createPolygonStruct([portEdge1(:,1) ; portEdge2([2 1],1)],[portEdge1(:,2) ; portEdge2([2 1],2)]);

         hold on;
         mod2D_showPolygon(gca,portPol,[1 0 0],'none');
         % quiver(mean(portEdge1(:,1)),mean(portEdge1(:,2)),mean(portEdge2(:,1)) - mean(portEdge1(:,1)),mean(portEdge2(:,2)) - mean(portEdge1(:,2)),'filled');
         hold off;

         break;
      else
         portEdge2 = portEdge1;
         switch tolower(in2)
            case 'u'
               e2Layer = e1Layer + 1;

               if editLayer == numel(metalLayers)
                  error('No upper layer available');
               end

               set(txtHdl,'string','Edge#1 - Edge#2 @ upper layer');
               break;
            case 'l'
               e2Layer = e1Layer - 1;

               if editLayer == 1
                  error('No lower layer available');
               end

               set(txtHdl,'string','Edge#1 - Edge#2 @ lower layer');
               break;
         end

         title('Wrong choice, try again.');
         pause(1.5);
      end
   end

   pause(1.5);

   close(gcf);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Store edges in container %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if e1Layer == e2Layer
      zCoor1 = elTable(metalLayers(e1Layer));
      zCoor2 = zCoor1;
   elseif e1Layer > e2Layer
      zCoor1 = elTable(metalLayers(e1Layer));
      zCoor2 = elTable(metalLayers(e2Layer)) + subThickTable(metalLayers(e2Layer));
   else
      zCoor1 = elTable(metalLayers(e1Layer)) + subThickTable(metalLayers(e1Layer));;
      zCoor2 = elTable(metalLayers(e2Layer));
   end

   portEdges = cat(3,portEdges,...
      [[portEdge1 ; portEdge2] [ones([2 1])*zCoor1 ; ones([2 1])*zCoor2]]);

end

% Store all ports


%ports = repmat([size(portEdges,3) 1],struct('portNum',[],'edge1',[],'edge2',[]));
ports(1).portNum = 0;
for pIdx = 1:size(portEdges,3)
   ports(pIdx).portNum = pIdx;
   ports(pIdx).edge1 = portEdges(1:2,:,pIdx);
   ports(pIdx).edge2 = portEdges(3:4,:,pIdx);
end

fprintf('Storing in ''model_ports.mat''\n');
save( fullfile(Sim_Path,'model_ports.mat'),...
      '-mat7-binary',...
      'ports');
