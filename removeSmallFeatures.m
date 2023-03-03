function p = removeSmallFeatures(p,th)

oP = p;
p.x = [nan ; p.x ; nan];
p.y = [nan ; p.y ; nan];

partEdgeIdxs = find(isnan(p.x));

smallFeatures = {};
for partIdx = 1:(numel(partEdgeIdxs) - 1)
   cRange = (partEdgeIdxs(partIdx) + 1):(partEdgeIdxs(partIdx + 1) - 1);
   cPart = [p.x(cRange) p.y(cRange)];

   % Check if part is NOT closed
   if ~min(cPart(1,:).' == cPart(end,:).')
      cPart = [cPart ; cPart(1,:)];
   end

   % Calculate each edges length
   edgeL = sqrt(sum((cPart(2:end,:) - cPart(1:(end-1),:)).^2,2));

   % Are there small features here?
   cBinSmallFeatures = edgeL < th;

   % Edge case 1: No small features
   if ~sum(cBinSmallFeatures)
      break;
   end

   cSmallFeatureList = cat(3,cPart(cBinSmallFeatures,:),cPart(find(cBinSmallFeatures) + 1,:));

   % Edge case 2: Entire list is a set of small features
   if sum(cBinSmallFeatures) == numel(cRange)
      smallFeatures{numel(smallFeatures) + 1} = cSmallFeatureList;
      break;
   end

   % Detect clusters and store them in sub-lists
   clusters = {};
   cCluster = [];
   for i = 1:numel(cBinSmallFeatures)
      if cBinSmallFeatures(i)
         cCluster = [cCluster ; cat(3,cPart(i,:),cPart(i + 1,:))];
      else
         if ~isempty(cCluster)
            clusters{numel(clusters) + 1} = cCluster;
            cCluster = [];
         end
      end
   end

   % Edge case 3: Last cluster connects to first cluster
   if cBinSmallFeatures(end)
      % Connect clusters and delete last one
      clusters{1} = [cCluster ; clusters{1}];
   end

   smallFeatures((numel(smallFeatures) + 1):(numel(smallFeatures) + numel(clusters))) = clusters;

end

p = oP;

% Iterate through features, query their removal
badFeatureCtr = 0;
for fIdx = 1:numel(smallFeatures)
   axs = [0,0];

   cFeat = smallFeatures{fIdx};

   bBox = [ min(reshape(cFeat(:,1,:),[],1)) max(reshape(cFeat(:,1,:),[],1)) ...
            min(reshape(cFeat(:,2,:),[],1)) max(reshape(cFeat(:,2,:),[],1))];

   bBoxDims = [(bBox(2) - bBox(1)) (bBox(4) - bBox(3))];


   % Sanity check: Does feature have volume?
   if ~bBoxDims(1) || ~bBoxDims(2)
      badFeatureCtr = badFeatureCtr + 1;
      continue;
   end

   clipPart = mod2D_createRectangleStruct([bBox(1) bBox(3)],[bBox(2) bBox(4)]);

   axisBox = bBox + [bBoxDims(1)*[-1 1]*2 bBoxDims(2)*[-1 1]*2];
   showPart = mod2D_createRectangleStruct([axisBox(1) axisBox(3)],[axisBox(2) axisBox(4)]);
   figure;
   cPos = get(gcf,'position');
   for plotIdx = 1:2
      subplot(1,2,plotIdx);

      mod2D_showPolygon(gca,p,[0 1 0],[0 0 0]);

      % Draw red lines
      cLine = [cFeat(:,:,1) ; cFeat(end,:,2)];
      hold 'on';
      plot(cLine(:,1),cLine(:,2),'-r','linewidth',4);

      if plotIdx == 1
         mod2D_showPolygon(gca,showPart,'none',[0 0 0]);
      else
         axis(axisBox);
         rectHdl = mod2D_showPolygon(gca,clipPart,'none',[0 0 0]);
      end

      hold off;
   end

   titleStr = 'Press ''a'' for ''Add'', ''c'' for ''Clip''';
   if size(cLine,1) > 2
      titleStr = [titleStr ', ''r'' for ''remove'''];
   end
   title([titleStr ', ''q'' to ignore']);

   hold 'off';

   while true
      [~,~,keyPressed] = ginput(1);

      % Make sure this isn't a mouse press
      if (keyPressed <= 5)
         continue;
      end

      keyPressed = tolower(keyPressed);

      boolActionTaken = false;
      switch keyPressed
         case 'a'
            boolActionTaken = true;
            tempP = mod2D_booleanOperation(p,clipPart,'add');
            break;
         case 'c'
            boolActionTaken = true;
            tempP = mod2D_booleanOperation(p,clipPart,'subtract');
            break;
         case 'r'
            boolActionTaken = true;
            tempP = p;
            rho = [tempP.x tempP.y];
            [featureIdxs,~] = find((p.x == cLine(:,1).') & (p.y == cLine(:,2).'));

            tempP.x(featureIdxs(2:(end-1))) = [];
            tempP.y(featureIdxs(2:(end-1))) = [];
            break;
         case 'q'
            break;
         otherwise
            continue;
      end



   end

   if boolActionTaken
      cla;
      mod2D_showPolygon(gca,tempP,[0 1 0],[0 0 0]);
      title('Happy? ''y'' for ''yes'', ''n'' for ''no''');


      while true
         [~,~,keyPressed] = ginput(1);

         % Make sure this isn't a mouse press
         if (keyPressed <= 5)
            continue;
         end

         keyPressed = tolower(keyPressed);

         switch keyPressed
            case 'y'
               p = tempP;
               break;
            case 'n'
               error('Unhappy with result, breaking...');
               break;
            otherwise
               continue;
         end

      end

   end
   close(gcf);

end

fprintf(1,'Could not remove %d features due to bad bounding boxes\n',badFeatureCtr);

% Function end
end
