function m = refineHexMesh(m,r,Lm,hardM)

% m      - Mesh in the given direction
% r      - Required ratio
% Lm     - Maximal sub-division size allowed
% hardM  - Hard mesh points
% First go from left to right

% Start from left to right
while true;

   D = diff(m(:));

   % This has limited accuracy, check within 5 per-cent
   binDivRatio =  D(2:end)./D(1:(end - 1)) > r*1.05;

%   fprintf(1,'Number of anomalies %d\n',sum(binDivRatio(:)));
%   pause(0.1);
   if ~max(binDivRatio)
      break;
   end

   % Start from left to right
   i1 = find(binDivRatio,1);

   % Find first hard mesh point that might brake the smoothing
   mHard = hardM(find(hardM > m(i1 + 1),1));
   iHard = find(m == mHard);

   % Take segment to smooth
   cM = [m(i1) m(i1 + 1) m(iHard)];
   i3 = iHard;

   cD = cM(2:end) - cM(1:(end-1));
   cR = cD(2:end)./cD(1:(end-1));

   % Refine this bit of mesh
   nM = AutoSmoothMeshLines(cM,Lm,r,'symmetric',0,'homogeneous',0);

   nD = nM(2:end) - nM(1:(end-1));
   nR = nD(2:end)./nD(1:(end-1));

%   figure('position',[1747    223   1420    811]);
%   subplot(2,1,1);
%   plot(cM,ones(size(cM)),'.r','markersize',14);
%   hold on;
%   plot(m,ones(size(m)),'.b','markersize',10);
%   hold off;
%   axis([(min(m(:))-5) (max(m(:))+5) 0.5 1.5]);
%
%   subplot(2,1,2);
%   plot(cM,ones(size(cM)),'.b','markersize',20);
%   hold on;
%   plot(nM,ones(size(nM)),'.r','markersize',10);
%   hold off
%   axis([(min(cM)-5) (max(cM)+5) 0.5 1.5]);
%
%   close(gcf);

   % Stick this snippet in
   m1 = m(1:(i1 - 1));
   m3 = m((i3 + 1):end);

   m = [m1 nM m3];

end



end
