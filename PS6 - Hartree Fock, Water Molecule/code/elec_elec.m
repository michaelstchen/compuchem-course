function elec_elec = elec_elec(al,bet,gam,del,RA,RB,RC,RD)

apb=al+bet;
abfac = al*bet/apb;
gpd=gam+del;
gdfac = gam*del/gpd;

dRAB = RA-RB;
dRAB2 = dRAB*dRAB';
dRCD = RC-RD;
dRCD2 = dRCD*dRCD';

RP = (al*RA+bet*RB)/apb;
RQ = (gam*RC+del*RD)/gpd;
dRPQ = RP-RQ;
dRPQ2 = dRPQ*dRPQ';

elec_elec = ((2*pi^(5/2))/(apb*gpd*sqrt(apb+gpd))) * ...
  exp(-abfac*dRAB2-gdfac*dRCD2) * ...
  F0((apb*gpd/(apb+gpd))*dRPQ2);