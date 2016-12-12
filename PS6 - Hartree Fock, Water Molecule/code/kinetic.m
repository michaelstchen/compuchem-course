function kinetic = kinetic(alpha,beta,RA,RB)

apb = alpha+beta;
abfac = alpha*beta/apb;
dRAB = RA-RB;
dRAB2 = dRAB*dRAB';

kinetic = ((pi/apb)^(3/2))*exp(-abfac*dRAB2)*abfac* ...
    (3-2*abfac*dRAB2);