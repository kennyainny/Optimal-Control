%% NOSESHAPE                  %%
%% Newton's optimal nose cone %%
%% shape using testshooting   %%

% And, replace x with t and r with x as compared to notes
clc;
clear all;
close all;
global x lambda;

% Numerical parameters
eps=0.01;     % test shooting
delta=0.05;   % cutoff for optimization
gamma=0.1;    % descent step size

T=1; % the maximum Length (denoted in notes by l)
dt=0.01; 

a=1;  
x0=a;   % initial condition

lambda0=1; % initial lambda0 guess
G=100;
G_STOR = [];
while (G>delta);
    t=0;
    X=[]; X=[X;x0];
    Lambda=[]; Lambda=[Lambda;lambda0];

    X=[X,x0];
    Lambda=[Lambda;lambda0+eps]; % Test shooting

    while(t<=T);
        U=[]; 
        for i=1:2;
            x=X(i);
            lambda=Lambda(i);
            % can't find closed loop solution for u
            u=fsolve('optunose',-1,optimset('Display','off'));
            dx=-u;
            dlambda=-u^3/(1+u^2);
            x=x+dt.*dx;
            lambda=lambda+dt.*dlambda;
            X(i)=x;
            Lambda(i)=lambda;
        end;
        t=t+dt;
    end;

    % Compute the derivative
    G=norm(Lambda(1)-(X(1)))^2;
    dG=[];
    Gi=norm(Lambda(2)-(X(2)))^2;
    dG=[(Gi-G)/eps];

    % gradient descent
    lambda0=lambda0-gamma.*dG;
    
    G_STOR = [G_STOR G];
end;
plot(G_STOR)
%%
% Optimal solution
t=0; TT=[]; X=[]; 
x=x0; lambda=lambda0;

while (t<=T);
    TT=[TT,t];
    X=[X,x];
    u=fsolve('optunose',-1,optimset('Display','off'));
    dx=-u;
    dlambda=-u^3/(1+u^2);
    x=x+dt.*dx;
    lambda=lambda+dt.*dlambda;
    t=t+dt;
end;
figure;
% Plot solution
hold off;
plot(TT,1.5-0.5.*X); % We used a transformation on x to get the nice-looking equations
hold on;
plot(TT,-(1.5-0.5.*X));
ylabel('x');
xlabel('t');
    
    
    
        