function gamma_min = gammaMin(u);
%initialize variables
b=[-2;1;-1;3;1];
Q=[5 0 8 -1 -3;0 10 9 7 11;8 9 25 0 6;-1 7 0 19 5;-3 11 6 5 18];

%g=0.5*(u.'*Q*u) + b.'*u;

dg= u.'*Q+b.';
gamma_min = dg'*dg*inv(dg*Q*dg');

