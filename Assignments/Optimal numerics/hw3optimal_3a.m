clc;
clear all;
close all;

% Numerical parameters
eps=0.01;     % test shooting
delta=0.001;   % cutoff for optimization
gamma=0.1;    % descent step size

T=1; dt=0.01;
S=[0 0 0;0 0 0;0 0 -2];
x0=[0;0;0];   % initial condition

lambda0=[0.1;0.1;0.1]; % initial lambda0 guess
G=100; % lambda0 cost
G_STOR = [];
while (G>delta);
    t=0;
    X=[]; X=[X,x0];
    Lambda=[]; Lambda=[Lambda,lambda0];
    
    % Run the test shooting - create variations of lambda0 
    % along N directions
    for i=1:3;
        ei=zeros(3,1); ei(i)= 1;
        X=[X,x0];
        Lambda=[Lambda,lambda0+eps.*ei];
    end;
    
    % Solve for lambda0 as well as for the perturbed
    % initial conditions
    while(t<=T);
        for i=1:4; 
            x=X(:,i);
            lambda=Lambda(:,i);
            
            u = [lambda(3)*x(2)-lambda(1); -lambda(3)*x(1)-lambda(2)];
            dx = [u(1); u(2); x(1)*u(2)-x(2)*u(1)];
            dlambda = [-lambda(3)*u(2); lambda(3)*u(1); 0];
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
    for i=2:4;
        Gi=norm(Lambda(:,i)-S*X(:,i))^2;
        dG=[dG,(Gi-G)/eps];
    end;
    % gradient descent
    lambda0=lambda0-gamma.*dG';
    G_STOR = [G_STOR G];
end;
plot(G_STOR);
% figure;
%%
% Compute the optimal solution
t=0; TT=[]; X=[];
x=x0; lambda=lambda0;
while (t<=T);
    TT=[TT,t];
    X=[X,x];
    u = [lambda(3)*x(2)-lambda(1); -lambda(3)*x(1)-lambda(2)];
    
    dx = [u(1); u(2); x(1)*u(2)-x(2)*u(1)];
    
    dlambda = [-lambda(3)*u(2); lambda(3)*u(1); 0];
    
    x=x+dt.*dx;
    lambda=lambda+dt.*dlambda;
    t=t+dt;
end;
X=[X,x];
% Plot solution
figure;
subplot(3,1,1)
plot(X(1,:));
ylabel('x_1');
subplot(3,1,2)
plot(X(2,:));
ylabel('x_2');
xlabel('t')
subplot(3,1,3)
plot(X(3,:));
ylabel('x_3');
xlabel('t')


