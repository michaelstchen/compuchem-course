% Calculates the harmonic component of the matrix element
% for two gaussian bases centered at XA and XB, with width ALPHA
function f = f(xA,xB,alpha)

deltax = xA-xB;
deltax2 = deltax^2;

sumx = xA+xB;
sumx2 = sumx^2;

alpha2 = alpha^2;

f = 0.5*s(xA,xB,alpha)*(alpha - alpha2*deltax2 ...
            + 0.25*(1/alpha + sumx2));

