function [elTable,subThickTable] = layerStruct_2_tables(layers)

% Convert layers to height map
elTable = zeros(size(layers));
subThickTable = zeros(size(layers));
for lIdx = 1:numel(layers)
   cLayer = layers(lIdx);
   subThickTable(lIdx) = subThickTable(lIdx) + cLayer.thickness;

   % If not last layer
   if lIdx < numel(layers)
      elTable((lIdx + 1):end) = elTable((lIdx + 1):end) + cLayer.thickness;
   end

   switch (cLayer.type)
      case 'conductor-top'
         if lIdx < numel(layers)
            elTable((lIdx + 1):end) = elTable((lIdx + 1):end) - cLayer.thickness;
            subThickTable(lIdx + 1) = subThickTable(lIdx + 1) + cLayer.thickness;
         end

      case 'conductor-bottom'
         if lIdx > 1
            subThickTable(lIdx - 1) = subThickTable(lIdx - 1) + cLayer.thickness;
         end
   end
end

end
