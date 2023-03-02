function [d,t] = dist2edge(node2,node1,vert)

% 1. Line difference vector
dr = node2 - node1;

% 2. Line coefficients according to line equation Ax + By + C = 0:
% dy*x-dx*y+dx*y1-dy*x1 = 0
A = dr(:,2);
B = -dr(:,1);
C = dr(:,1).*node1(:,2) - dr(:,2).*node1(:,1);

% 3. Determine distances between all vertices and all lines:
% |Ax + By + C|/sqrt(A^2+B^2)
d = bsxfun(@plus, bsxfun(@times,A.',vert(:,1)) + ...
                  bsxfun(@times,B.',vert(:,2)),...
                  C.');
d = bsxfun(@times,abs(d),1./sqrt((A.^2).' + (B.^2).'));

% Check normalized distance from node #1
% 1. Determine edge lengths
L = sqrt(sum(dr.^2,2));

% 2. Normalize direction vector
dl = bsxfun(@times,dr,1./L);

% 3. Project (normalized coordinates, between 0 and 1)
node_offset = bsxfun(@minus,permute(vert,[1 3 2]),permute(node1,[3 1 2]));
t = sum(bsxfun(@times,node_offset,permute(dl,[3 1 2])),3);
t = bsxfun(@times,t,1./(L.'));

end
