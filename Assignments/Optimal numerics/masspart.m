%% MASSPARTICLE                  %%
clc;
clear all;
close all;
h=0.04;

% Numerical parameters
eps=0.01;     % test shooting
delta=0.1;   % cutoff for optimization
gamma=0.1;    % descent step size

T=1; dt=0.01; 

% let x=[x,y,vx,vy]'
x0=zeros(4,1);   % initial condition

lambda0=[0;0.4;1;0.2]; % initial lambda0 guess
G=100;
G_STOR = [];
while (G>delta);
    t=0;
    X=[]; X=[X,x0];
    Lambda=[]; Lambda=[Lambda,lambda0];
    
    % Run the test shooting
    for i=1:4;
        ei=zeros(4,1); ei(i)=1;
        X=[X,x0];
        Lambda=[Lambda,lambda0+eps.*ei];
    end;

    % Solve for lambda0 as well as for the perturbed
    % initial conditions
    while(t<=T); 
        for i=1:5; % Position 1 is the nominal system 
            x=X(:,i);
            lambda=Lambda(:,i);
            u=atan2(lambda(4),lambda(3));
            dx=[x(3);x(4);cos(u);sin(u)];
            dlambda=[0;0;-lambda(1);-lambda(2)];
            x=x+dt.*dx;
            lambda=lambda+dt.*dlambda;
            X(:,i)=x;
            Lambda(:,i)=lambda;
        end;
        t=t+dt;
    end;

    %Compute the gradient dg/dlambda0
    G=(x(2)-h)^2+x(4)^2+lambda(1)^2+(lambda(3)+1)^2;
    dG=[];
    for i=2:5;
        lambda=Lambda(:,i); x=X(:,i);
        Gi=(x(2)-h)^2+x(4)^2+lambda(1)^2+(lambda(3)+1)^2;
        dG=[dG,(Gi-G)/eps];
    end;

    % gradient descent
    lambda0=lambda0-gamma.*dG';
    G_STOR = [G_STOR G];
end;
plot(G_STOR);
figure;
%%
% Compute the optimal solution
t=0; TT=[]; X=[]; 
x=x0; lambda=lambda0;
while (t<=T);
    TT=[TT,t];
    X=[X,x];
    u=atan2(lambda(4),lambda(3));
    dx=[x(3);x(4);cos(u);sin(u)];
    dlambda=[0;0;-lambda(1);-lambda(2)];
    x=x+dt.*dx;
    lambda=lambda+dt.*dlambda;
    t=t+dt;
end;
X=[X,x];
% Plot solution
hold off;
plot(X(1,:),X(2,:));
    
    
        