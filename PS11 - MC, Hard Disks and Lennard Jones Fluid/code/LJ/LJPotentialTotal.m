function [ U ] = LJPotentialTotal( N, L, r )

    U = 0;
    for i = 1:N
        for j = (i+1):N
            deltar = r(:,i)-r(:,j);
            deltar = deltar - L*round(deltar/L);
            deltar2 = deltar'*deltar;
            dr = sqrt(deltar2);
            
            U = U + LJPotential(dr, 2.5);
        end
    end

end

