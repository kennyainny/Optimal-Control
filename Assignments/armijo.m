function[g_armijo,dg_armijo]=armijo(Uk);
%initialize variables
Q = [ 5 0 8 -1 -3; 0 10 9 7 11; 8 9 25 0 6; -1 7 0 19 5; -3 11 6 5 18];
B = [-2; 1; -1; 3; 1];

g_armijo = 0.5 * transpose(Uk) * Q * Uk + transpose(B)*Uk;
dg_armijo = transpose(Uk)*Q + transpose(B);

%a=3; b=1; c=-2; d=4;
%g=a*u^3+b*u^2+c*u+d;
%dg=3*a*u^2+2*b*u+c;
