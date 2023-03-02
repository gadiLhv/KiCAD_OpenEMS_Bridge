function [layerPair,D,loc,plating,type] = parseDrillFile(fileName)

fHdl = fopen(fileName);
% First, look for the layer pair line
while true
   cLine = fgetl(fHdl);
   [vals, pos, errMsg] = textscan(cLine,'; #@! TF.%[^,],%[^,],%d,%d,%s');

   command = testEmptyCell(vals{1});

   % Try and parse line
   plating = testEmptyCell(vals{2});

   layerTop = vals{3};
   layerBot = vals{4};

   type = testEmptyCell(vals{5});

   if strcmp(command,'FileFunction') && ~isempty(plating) && layerTop && layerBot && ~isempty(type)
      break;
   end

end

% Now look for "METRIC" line
while ~strcmp(cLine,'METRIC')
   cLine = fgetl(fHdl);
end

% If this is a multiple layer drill, make list of drills
pairSet = (layerTop:(layerBot - 1)).';
pairSet = [pairSet (pairSet + 1)];

% First column is drill ID, 2nd is drill diameter, in mm
drillD = [];
cLine = fgetl(fHdl); % Initialize first line. Instead of doing do\while
while ~strcmp(cLine,'%')
   vals = sscanf(cLine,'T%dC%f');

   drillD = [drillD ; vals.'];
   cLine = fgetl(fHdl);
end


% Now look for drill titles and drills
D = [];
layerPair = [];
loc = [];

currentlyReading = false;
cLine = fgetl(fHdl);
cD = 0;  % Current drill size
while ~feof(fHdl) && ~isempty(drillD)


   if ~currentlyReading
      % If still looking for title
      [val, count, errmsg, pos] = sscanf(cLine,'T%d');
      if (count == 1) && isempty(errmsg) && (pos == (length(cLine) + 1))
         currentlyReading = true;
         cD = drillD(find(drillD(:,1) == val),2);
      end
   else
      % If found title and reading pairs
      [val, count, errmsg, pos] = sscanf(cLine,'X%fY%f');
      if (count == 2) && isempty(errmsg) && (pos == (length(cLine) + 1))
         % Find proper drill ID

         % Put drill and pair(s) in list
         D = [D ; ones(size(pairSet,1))*cD];
         loc = [loc ; [ones(size(pairSet,1))*val(1) ones(size(pairSet,1))*val(2)]];
         layerPair = [layerPair ; pairSet];
      else
         % Mismatch means in the next title
         currentlyReading = false;
         cD = 0;
         continue;
      end
   end
   cLine = fgetl(fHdl);
end

fclose(fHdl);


end

function c = testEmptyCell(c)

if isempty(c)
   c = [];
else
   c = c{1};
end

end

