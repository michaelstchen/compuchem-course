function [ F, U, h ] = LJForcePot( r, rcut, L, h )

    U = 0;
    F = zeros(2, length(r));
    for i = 1:length(r)
        for j = (i+1):length(r)
            deltar = r(:,i)-r(:,j);
            deltar = deltar - L*round(deltar/L);
            deltar2 = deltar'*deltar;

            dr = sqrt(deltar2);
            h = histo(h, dr);
            if (sqrt(deltar2) <= rcut)
                dr2inv = 1/deltar2;
                dr4inv = dr2inv^2;
                dr6inv = dr2inv*dr4inv;
                dr8inv = dr6inv*dr2inv;
                dr12inv = dr6inv*dr6inv;
                dr14inv = dr12inv*dr2inv;

                uprimer = 48*(dr14inv - 0.5*dr8inv);

                F(:,i) = F(:,i) + uprimer*deltar;
                F(:,j) = F(:,j) - uprimer*deltar;

                U = U + 4*(dr12inv - dr6inv);
            end
        end
    end

end

