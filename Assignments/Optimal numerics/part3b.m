clear all;
clc;

% Numerics
t=0; dt=0.01; tf=1;
NL=[]; k=250; JJ=[];
delta=0.001;
% Initial guess
z0=[0;0;0;0.1;0.1;0.1];
G=100;
G_STOR=[];
% Gradient descent step size
gamma=0.01;

while (G>delta)
    Z=[]; J=0; i=0; z=z0;
    while (t<=tf);  % solve for z forward in time
        Z=[Z,z]; 
        u1=-z(4)+z(6)*z(2);
        u2=-z(5)-z(6)*z(1);
        dz=[u1;u2;z(1)*u2-z(2)*u1;-u2*z(6);u1*z(6);0];
        %meu=[-lambda(3)*u2;lambda(3)*u1;0;0;0;lambda(4)*u2+lambda(5)*u1]
        %meu=meu+dt.*dmeu
        z=z+dt.*dz; % forward Euler
        t=t+dt; i=i+1;
    end;
    G=(z(4))^2+(z(5))^2+(z(6)+2*z(3))^2;
    t=tf; 
%     meu=[0;0;-z(3);0.25;0.25;0.2];
%     meu=[0;0;z(3)*2+z(6);z(4);z(5);z(6)];
    meu=[0;0;z(3)*2+z(6);z(4);z(5);z(3)*2+z(6)];
    while (t>=0);   % solve for meu backwards in time
        u1=-z(4)+z(6)*z(2);
        u2=-z(5)-z(6)*z(1);
        dmeu=[-meu(3)*u2;meu(3)*u1;0;0;0;meu(4)*u2+meu(5)*u1];
        meu=meu-dt.*dmeu;% backwards Euler
        t=t-dt; 
        z=Z(:,i);
        i=i-1;
    end;
    z0=z0-gamma.*meu;
    z0=[0;0;0;z0(4:6)]
    G_STOR = [G_STOR G];
end;

plot(G_STOR);
