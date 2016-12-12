function [ u ] = LJPotential( r, rcut )
    if (r > rcut)
        u = 0;
    else
       u = 4 * ((r^-12 - r^-6) - (rcut^-12 - rcut^-6)); 
    end
end

