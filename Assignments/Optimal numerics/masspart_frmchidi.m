%% MASSPARTICLE                  %%
clc;
clear all;
close all;
%h=0.04;

% Numerical parameters
eps=0.01;     % test shooting
delta=0.031;   % cutoff for optimization
gamma=0.1;    % descent step size
p = 2;

T=1; dt=0.01; 

% let x=[x,y,vx,vy]'
x0=zeros(3,1);   % initial condition

lambda0=[0.1;0.1;0.1]; % initial lambda0 guess
G=100;
G_STOR = [];
while (G>delta);
    t=0;
    X=[]; X=[X,x0];
    Lambda=[]; Lambda=[Lambda,lambda0];
    
    % Run the test shooting
    for i=1:3;
        ei=zeros(3,1); ei(i)=1;
        X=[X,x0];
        Lambda=[Lambda,lambda0+eps.*ei];
    end;

    % Solve for lambda0 as well as for the perturbed
    % initial conditions
    while(t<=T); 
        for i=1:4; % Position 1 is the nominal system 
            x=X(:,i);
            lambda=Lambda(:,i);
            %u=atan2(lambda(4),lambda(3)); 
            u1 = -lambda(1)+lambda(3)*x(2);
            u2 = -lambda(2) - lambda(3)*x(1);
            dx=[u1;u2;x(1)*u2-x(2)*u1];
            dlambda=[-lambda(3)*u2;lambda(3)*u1;0];
            x=x+dt.*dx;
            lambda=lambda+dt.*dlambda;
            X(:,i)=x;
            Lambda(:,i)=lambda;
        end;
        t=t+dt;
    end;
    %x=X(:,1);
    %Compute the gradient dg/dlambda0
    G=(lambda(1)^2+lambda(2)^2+(lambda(3)+ p*x(3))^2);
    dG=[];
    for i=2:4;
        lambda=Lambda(:,i); x=X(:,i);
        Gi=(lambda(1)^2+lambda(2)^2+(lambda(3)+ p*x(3))^2);
        dG=[dG,(Gi-G)/eps];
    end;

    % gradient descent
    lambda0=lambda0-gamma.*dG';
    G_STOR = [G_STOR G];
end;
plot(G_STOR);
%figure;
%%
% Compute the optimal solution

    
    
        