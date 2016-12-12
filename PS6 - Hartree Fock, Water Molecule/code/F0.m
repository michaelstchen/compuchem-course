function F0 = F0(x)

if (x<0.00001)
    F0=1-x/3;
else
    sqx = sqrt(x);
    F0 = (sqrt(pi)/2)*erf(sqx)/sqx;
end