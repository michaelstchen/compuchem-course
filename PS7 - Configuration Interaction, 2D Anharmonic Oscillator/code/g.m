% Calculates the matrix element associated with the anharmonic
% term for two gaussian bases centered at XA and XB, with width ALPHA
function g = g(xA,xB,alpha)

deltax = xA-xB;
deltax2 = deltax^2;

xP = 0.5*(xA+xB);
xP2 = xP^2;
xP4 = xP2^2;

alpha2 = alpha^2;
g = s(xA,xB,alpha)*( 3/(16*alpha2) + 1.5*xP2/alpha + xP4);
