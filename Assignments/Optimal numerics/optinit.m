% System and cost parameters
A=[1,2;2,-1];
Q=eye(2); q=1;

% Gradient descent step size
gamma=0.01;

% Numerics
t=0; dt=0.01; tf=1;
NL=[]; k=10; JJ=[];

% Initial guess
x0=[1;1];

for (j=1:k);
    X=[]; J=0; i=0; x=x0;
    while (t<=tf);  % solve for x forward in time
        X=[X,x]; 
        J=J+dt*(x'*Q*x-q);
        x=x+dt.*(A*x); % forward Euler
        t=t+dt; i=i+1;
    end;
    t=tf; lambda=[0;0];
    while (t>=0);   % solve for lambda backward in time
        lambda=lambda-dt.*(-2.*Q*X(:,i)-A'*lambda); % backwards Euler
        t=t-dt; i=i-1;
    end;
    x0=x0-gamma.*lambda;
    NL=[NL;norm(lambda)];
    JJ=[JJ;J];
end;

hold off;
figure(1);
subplot(2,1,1);
plot(JJ);
subplot(2,1,2);
plot(NL);
