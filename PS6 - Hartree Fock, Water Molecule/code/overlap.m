function overlap = overlap(alpha,beta,RA,RB)

apb = alpha+beta;
abfac = alpha*beta/apb;
dRAB = RA-RB;
dRAB2 = dRAB*dRAB';

overlap = ((pi/apb)^(3/2))*exp(-abfac*dRAB2);