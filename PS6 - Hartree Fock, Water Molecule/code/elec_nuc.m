function elec_nuc = elec_nuc(alpha,beta,RA,RB,RC)

apb = alpha+beta;
abfac = alpha*beta/apb;
dRAB = RA-RB;
dRAB2 = dRAB*dRAB';

RP = (alpha*RA+beta*RB)/apb;
dRPC = RP-RC;
dRPC2 = dRPC*dRPC';

elec_nuc = -(2*pi/apb)*exp(-abfac*dRAB2)*F0(apb*dRPC2);