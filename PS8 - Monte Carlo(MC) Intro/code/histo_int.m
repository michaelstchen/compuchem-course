% Problem 1 part (ii)
function [ h ] = histo_int( h, dat )

if (h.count==0)     
  h.numbins = h.range(2) - h.range(1) + 1;
  h.hist = zeros(1,h.numbins);
  h.vals = 1:h.numbins;
  h.vals = h.range(1) + h.vals - 1;
end

if (dat >= h.range(1) && dat <= h.range(2))
  bin = dat - h.range(1) + 1;
  h.hist(bin) = h.hist(bin) + 1;
  h.count = h.count+1;

end

