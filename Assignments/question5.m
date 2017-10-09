clc;
clear all;

syms u(t) v(t)
ode1 = diff(u) == v;
ode2  = diff(v) == -u -1;
odes = [ode1; ode2]

% odes(t) =
%  diff(u(t), t) == v(t)
%  diff(v(t), t) == -u(t)
 
 
S = dsolve(odes)
uSol(t) = S.u
vSol(t) = S.v

%cond1 = u(1) == 0;
%cond2 = v(1) == 1;
%conds = [cond1; cond2];
%[uSol(t), vSol(t)] = dsolve(odes,conds)