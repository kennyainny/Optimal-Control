clc;
clear;
syms u(t) v(t)
ode1 = diff(u) == v;
ode2  = diff(v) == -u;
odes = [ode1; ode2]

% odes(t) =
%  diff(u(t), t) == v(t)
%  diff(v(t), t) == -u(t)
 
 
S = dsolve(odes)
uSol(t) = S.u
vSol(t) = S.v