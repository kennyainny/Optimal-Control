function[g,dg]=fundata(u);
a=3; b=1; c=-2; d=4;
g=a*u^3+b*u^2+c*u+d;
dg=3*a*u^2+2*b*u+c;
