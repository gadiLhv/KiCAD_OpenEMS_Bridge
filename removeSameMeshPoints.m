function x = removeSameMeshPoints(x,th)

   binIsSame = true;
   while max(binIsSame = (diff(x(:)) < th))
      sameIdx = find(binIsSame,1);
      x(sameIdx) = mean(sameIdx + [0 ; 1]);
      x(sameIdx + 1) = [];
   end
end
