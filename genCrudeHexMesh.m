function m = genCrudeHexMesh(hardM,r,Lm)

% Start left to right

diffPairs = [hardM(1:(end-1)).' hardM(2:end).'];
diffPairIdxs = [(2:numel(hardM)).' (1:(numel(hardM)-1)).'];

Ds = diffPairs(:,1) - diffPairs(:,2);

m = [];

% This generates equally spaced meshes
for dIdx = 1:size(diffPairs,1)
   [cM q] = AutoSmoothMeshLines(diffPairs(dIdx,:),Lm,r);
   m = [m cM];

end


end
