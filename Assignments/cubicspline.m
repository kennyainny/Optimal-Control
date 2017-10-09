% Cubic splines over [0,1];

% Initial and final condition coefficients
% M*c=v
Mx10=[0,0,0,1];
Mx11=[1/6,-1/2,1,1];
Mx20=[0,0,1,0];
Mx21=[1/2,-1,1,0];
Ml10=[1,0,0,0];
Ml11=[1,0,0,0];
Ml20=[0,1,0,0];
Ml21=[-1,1,0,0];

cs=2;

if (cs==1);
% If all boundary points (in x) are given
    M=[Mx10;Mx11;Mx20;Mx21];
    v=[0;1;-1;0];
elseif (cs==2);
% If only end positions but no velocities are given
    M=[Mx10;Mx11;Ml20;Ml21];
    v=[0;1;0;0];
elseif (cs==3);
% If initial position and velocity is given and final position
    M=[Mx10;Mx11;Mx20;Ml21];
    v=[0;2;1;0];
end;
% Compute the coefficients
c=inv(M)*v;

t=0:0.01:1;
p=c(1)/6.*t.^3-c(2)/2.*t.^2+c(3).*t+c(4);

plot(t,p);






