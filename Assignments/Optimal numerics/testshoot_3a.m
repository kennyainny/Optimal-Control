%% TESTSHOOT Problem 3a %%
close all;
clear all;
clc;
format short; format compact;

% Cost matrices and system dynamics
%Q=eye(2); S=0.5.*eye(2); R=1;
%A=[0,1;0,0]; B=[0;1];

% Various constants
eps=0.01;   % step in testshooting
delta=0.001; % cutoff for gradient descent
gamma=0.1;  % gradient descent step size

n=3; T=1; dt=0.01; 
x0=[0;0;0];
lambda0=[0.1;0.1;0.1];  % initial guess
G = 100;  % lambda0 cost
G_STOR = [];

%x_31=-lambda0(3)/2; x_1=x0(1); x_2=x0(2);
l3=lambda0(3); x_1=x0(1); x_2=x0(2);

A = [0, l3, 0;
     -l3, 0, 0; 
     -l3*x_1, -l3*x_2, 0];

B = [-1, 0, 0;
     0, -1, 0;
     x_2, -x_1, 0];
 
C = [l3^2, 0, 0;
     0, l3^2, 0;
     0, 0, 0];

D = [0, l3, 0;
     -l3, 0, 0;
     0, 0, 0];
 
 
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

    % Solve forward in time for each perturbed case in X and lambda
    while(t<=T); 
        for i=1:n+1; % Position 1 is the nominal system 
            x=X(:,i);
            lambda=Lambda(:,i);   
            
            %x_31=-lambda(3)/2; x_1=x(1); x_2=x(2);
            l3=lambda(3); x_1=x(1); x_2=x(2);
            dx = A*x+B*lambda;
            dlambda = C*x+D*lambda;
            
            x=x+dt.*dx;
            lambda=lambda+dt.*dlambda;
            
            X(:,i)=x;
            Lambda(:,i)=lambda;
        end;
        t=t+dt;
    end;

    %Compute the gradient dg/dlambda0
    %G=norm(Lambda(:,1)-2*X(:,1))^2;
    G=(lambda(1)^2+lambda(2)^2+(lambda(3)+2*x(3))^2);
    dG=[];
    for i=2:n+1;
        lambda=Lambda(:,i); x=X(:,i);
        Gi=(lambda(1)^2+lambda(2)^2+(lambda(3)-2*x(3))^2);
        %Gi=norm(Lambda(:,i)-2*X(:,i))^2;
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
    
    
    
        