n=5; m=1; 

tf=1; dt=0.001;

% Generate a random system and 
% boundary conditions
A=randn(n,n);
B=randn(n,m);

x0=randn(n,1); 
xT=randn(n,1);

% Generate initial costate value

M=[A,-B*B';zeros(n,n),-A'];
N=expm(tf.*M);
N1=N(1:n,1:n);
N2=N(1:n,n+1:2*n);
lambda0=inv(N2)*(xT-N1*x0);

t=0; 
% time and output vectors
T=[]; Y=[];

x=x0; lambda=lambda0;
while (t<=tf);
    T=[T;t]; Y=[Y;x(1)];
    
    u=-B'*lambda;
    x=x+dt.*(A*x+B*u);
    lambda=lambda+dt.*(-A'*lambda);
    t=t+dt;
    
end;

hold off;
plot(T,Y);
hold on;
plot(0,x0(1),'o'); plot(tf,xT(1),'o');
    
    