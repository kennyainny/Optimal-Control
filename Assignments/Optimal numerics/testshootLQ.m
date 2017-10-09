%% TESTSHOOTLQ %%
close all;
clear all;
clc;
format short; format compact;

% Cost matrices and system dynamics
Q=eye(2); S=0.5.*eye(2); R=1;
A=[0,1;0,0]; B=[0;1]; n=2;

% Variou s constants
eps=0.01;   % step in testshooting
delta=0.01; % cutoff for gradient descent
gamma=0.1;  % gradient descent step size

T=1; dt=0.01; 
x0=[1;1];
lambda0=[0;0];  % initial guess
G=100;  % lambda0 cost
G_STOR = [];
while (G>delta);
    t=0;
    X=[]; X=[X,x0];
    Lambda=[]; Lambda=[Lambda,lambda0];

    % Test shooting- create variations of lambda0 along N directions
    for i=1:n;
        ei=zeros(n,1); ei(i)=1;
        X=[X,x0];
        Lambda=[Lambda,lambda0+eps.*ei];
    end;

    % Solve forward in time for each perturbed case
    while(t<=T); 
        for i=1:n+1; % Position 1 is the nominal system 
            x=X(:,i);
            lambda=Lambda(:,i);
            u=-inv(R)*B'*lambda;
            dx=A*x+B*u;
            dlambda=-Q*x-A'*lambda;
            x=x+dt.*dx;
            lambda=lambda+dt.*dlambda;
            X(:,i)=x;
            Lambda(:,i)=lambda;
        end;
        t=t+dt;
    end;

    %Compute the gradient dg/dlambda0
    G=norm(Lambda(:,1)-S*X(:,1))^2;
    dG=[];
    for i=2:n+1;
        Gi=norm(Lambda(:,i)-S*X(:,i))^2;
        dG=[dG,(Gi-G)/eps];
    end;

    % gradient descent
    lambda0=lambda0-gamma.*dG';

    G_STOR = [G_STOR G];
end;
plot(G_STOR);
disp(lambda0)
%% Compute the optimal solution
t=0; TT=[]; X=[]; 
x=x0; lambda=lambda0;

while (t<=T);
    TT=[TT,t];
    X=[X,x];
    u=-inv(R)*B'*lambda;
    dx=A*x+B*u;
    dlambda=-Q*x-A'*lambda;
    x=x+dt.*dx;
    lambda=lambda+dt.*dlambda;
    t=t+dt;
end;

% Plot the solution
figure;
subplot(2,1,1)
plot(TT,X(1,:));
ylabel('x_1');
subplot(2,1,2)
plot(TT,X(2,:));
ylabel('x_2');
xlabel('t');
    
    
    
        