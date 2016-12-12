% Calculates overlap matrix element for two gaussian
% bases with peaks at XA and XB, and width/stdev ALPHA
function s = s(xA,xB,alpha)

deltax = xA-xB;
deltax2 = deltax^2;

s = sqrt(pi/(2*alpha))*exp(-0.5*alpha*deltax2);