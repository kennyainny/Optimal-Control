% Gradient descent step size
gamma=0.01;
rho=2;
% Numerics
t=0; dt=0.01; tf=1;
NL=[]; k=200; JJ=[];
delta=0.001;
% Initial guess
z0=[0;0;0;0.1;0.1;0.1];
G=100;
G_STOR=[];
while (G>delta)
    Z=[]; J=0; i=0; z=z0;
    %lambda=[0;0;0;0;0;0];
    while (t<=tf);  % solve for x forward in time
        Z=[Z,z]; 
        %J=J+dt*(x'*Q*x-q);
        u1=-z(4)+z(6)*z(2);
        u2=-z(5)-z(6)*z(1);
        dz=[u1;u2;z(1)*u2-z(2)*u1;-u2*z(6);u1*z(6);0];
        %dlambda=[-lambda(3)*u2;lambda(3)*u1;0;0;0;lambda(4)*u2+lambda(5)*u1]
        %lambda=lambda+dt.*dlambda
        z=z+dt.*dz; % forward Euler
        t=t+dt; i=i+1;
    end;
    G=(z(4))^2+(z(5))^2+(z(6)+rho*z(3))^2;
    t=tf; 
%     lambda=[0;0;-z(3);0.25;0.25;0.2];
%     lambda=[0;0;z(3)*rho+z(6);z(4);z(5);z(6)];
    lambda=[0;0;z(3)*rho+z(6);z(4);z(5);z(3)*rho+z(6)];
    while (t>=0);   % solve for lambda backward in time
        u1=-z(4)+z(6)*z(2);
        u2=-z(5)-z(6)*z(1);
        dlambda=[-lambda(3)*u2;lambda(3)*u1;0;0;0;lambda(4)*u2+lambda(5)*u1];
        lambda=lambda-dt.*dlambda;% backwards Euler
        t=t-dt; 
        z=Z(:,i);
        i=i-1;
    end;
    z0=z0-gamma.*lambda;
    z0=[0;0;0;z0(4:6)]
    G_STOR = [G_STOR G];
    %NL=[NL;norm(lambda)];
    %JJ=[JJ;J];
end;

%hold off;
%figure(1);
%subplot(2,1,1);
%plot(JJ);
%subplot(2,1,2);
%plot(NL);
%figure(2);
plot(G_STOR);
