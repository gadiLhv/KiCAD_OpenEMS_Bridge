function x = reducePtsRes(x,minD)

   i0 = 1;
   i1 = 2;
   x0 = x(1);
   binRemove = false(size(x));
   while max(i0,i1) <= numel(x)
      x0 = x(i0);
      x1 = x(i1);

      if (x1 - x0) < minD
         binRemove(i1) = true;
         i1 = i1 + 1;
      else
         i0 = i1;
         i1 = i0 + 1;
      end
   end

   x(binRemove) = [];
end
