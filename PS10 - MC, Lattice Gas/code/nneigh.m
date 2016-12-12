function nneigh = nneigh(L,n,x,y)

left = x-1;
if (left==0)
  left = L;
end

right = x+1;
if (right==L+1)
  right=1;
end

up = y+1;
if (up==L+1)
  up=1;
end

down = y-1;
if (down==0)
  down=L;
end

nneigh = n(left,y) + n(right,y) + n(x,up) + n(x,down);