%% Questions 3a - Testshooting                  %%
clc;
clear all;
close all;

% Numerical parameters
eps=0.01;     % test shooting
delta=0.01;   % cutoff for optimization
gamma=0.1;    % descent step size

n=3; T=1; dt=0.01; k=0;
x0=zeros(n,1);   % initial condition
lambda0=[0.1;0.1;0.1]; % initial lambda0 guess

G=100;
G_STOR = [];

while (G>delta);
    t=0;
    X=[]; X=[X,x0];
    Lambda=[]; Lambda=[Lambda,lambda0];
    
    % Run the test shooting
    for i=1:n;
        ei=zeros(n,1); ei(i)=1;
        X=[X,x0];
        Lambda=[Lambda,lambda0+eps.*ei];
    end;

    % Solve for lambda0 as well as for the perturbed
    % initial conditions
    while(t<=T); 
        for i=1:n+1; % Position 1 is the nominal system 
            x=X(:,i);
            lambda=Lambda(:,i);
            
            u1=lambda(3)*x(2) - lambda(1);
            u2=-lambda(3)*x(1) - lambda(2);
                        
            dx=[u1;u2;(x(1)*u2-x(2)*u1)];            
            dlambda=[(x(1)*lambda(3)^2 + lambda(3)*lambda(2));
                     (x(2)*lambda(3)^2 - lambda(3)*lambda(1));
                     0];
            
            x=x+dt.*dx;
            lambda=lambda+dt.*dlambda;
            X(:,i)=x;
            Lambda(:,i)=lambda;
        end;
        t=t+dt;
    end;

    %Compute the gradient dg/dlambda0
    G=0.5*(lambda(1)^2+lambda(2)^2+(lambda(3)+2*x(3))^2);
    dG=[];
    for i=2:n+1;
        lambda=Lambda(:,i); x=X(:,i);
        Gi=0.5*(lambda(1)^2+lambda(2)^2+(lambda(3)+2*x(3))^2);
        dG=[dG,(Gi-G)/eps];
    end;

    % gradient descent: variable step-size
    k=k+1;
    stepdir=gamma/k;
    lambda0=lambda0-stepdir.*dG';
    G_STOR = [G_STOR G];
end;
plot(G_STOR);
    