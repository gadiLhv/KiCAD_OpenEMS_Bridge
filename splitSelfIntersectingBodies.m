function pols = splitSelfIntersectingBodies(pols)
% Splits polygons with self intesections (holes) so they can be accepted by the
% CSXCAD environment


pIdx = 1;
while pIdx <= numel(pols)
   cPol = pols{pIdx};

   % If does not containt self-intersections
   if ~max(isnan(cPol.x))
      pIdx = pIdx + 1;
      continue;
   end

   cX = cPol.x;
   cY = cPol.y;

   cX = [cX ; nan];
   cY = [cY ; nan];

   % Always deal with one self intersection at a time
   siIdxs = find(isnan(cX));

   % Determine bounding box of polygon
   bBox = [min(cX) max(cX) min(cY) max(cY)];

   cents = zeros([numel(siIdxs)-1 2]);
   centsTR = cents;
   for holeIdx = 1:(numel(siIdxs) - 1)
      cents(holeIdx,:) = mean([cX((siIdxs(holeIdx)+1):(siIdxs(holeIdx + 1)-1)) cY((siIdxs(holeIdx)+1):(siIdxs(holeIdx + 1)-1))]);
   end

   % Find closest centeroid to bottom left corner
   [~,minIdx] = min((bBox(1) - cents(:,1)).^2 + (bBox(3) - cents(:,2)).^2);

%   warning off;
   % Check weather this is a hole or just a 2nd part
   testPol = mod2D_createPolygonStruct(...
               cX((siIdxs(minIdx) + 1):(siIdxs(minIdx+1) - 1)),...
               cY((siIdxs(minIdx) + 1):(siIdxs(minIdx+1) - 1)));
   subPol = mod2D_booleanOperation(cPol,testPol,'Subtract');
%   warning on;

   if cPol.nParts == subPol.nParts;
      % If the number of loops didn't change, this is a hole. In this case clip a part

      % Create a rectangle between bottom left corner and closest centroid, then
      % clip it out of this polygon.
%      clipBox = mod2D_createRectangleStruct([bBox(1) bBox(3)],cents(minIdx,:));
      clipBox = mod2D_createRectangleStruct([bBox(1) bBox(3)],[cents(minIdx,1) bBox(4)]);

      % Intersect the two bodies
%      warning off;
      pol2 = mod2D_booleanOperation(cPol,clipBox,'Intersect');
      cPol = mod2D_booleanOperation(cPol,clipBox,'Subtract');
%      warning on;

   else
      % The other option is that the number of loops is now smaller, meaning part of the object was clipped out

      cPol = subPol;
      pol2 = testPol;

   end

   pols{pIdx} = cPol;
   % Store in new container
   newPols = cell([1 numel(pols)+1]);
   newPols(1:pIdx) = pols(1:pIdx);
   newPols(pIdx+1) = pol2;
   newPols((pIdx+2):end) = pols((pIdx+1):end);
   pols = newPols;

%   fprintf('Polygon %d/%d disected\n',pIdx,numel(pols));

%   figure;
%   ax = gca;
%   mod2D_showPolygon(ax,cPol,[0 0 1],[0 0 1]);
%   mod2D_showPolygon(ax,pol2,[1 0 0],[1 0 0]);
%   close(gcf);


end

