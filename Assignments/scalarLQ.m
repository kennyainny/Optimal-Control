a=1; b=2; q=1; r=0.5; s=2;
tf=1; dt=0.01;

x0=1.5;

P=[];
p=s; t=tf;

while (t>=0);
    P=[p;P];
    p=p-dt*(-2*a*p-q+b^2/r*p^2);
    t=t-dt;
end;
    
x=x0; t=0;
X=[]; U=[]; T=[]; i=1;
while (t<=tf);
    X=[X;x];
    T=[T;t];
    p=P(i);
    u=-b*p/r*x;
    U=[U;u];
    x=x+dt*(a*x+b*u);
    t=t+dt;
    i=i+1;
end;
    
figure(1);
hold off;
plot(T,X);
figure(2);
hold off;
plot(T,U);

 % Analytic solution
 M=[a,-b^2/r;-q,-a];
 x2=x0; t=0;
 X2=[]; U2=[];
 while (t<=tf);
     X2=[X2;x2];
     z=expm(M.*(t-tf))*[1;s];
     p=z(2)/z(1);
     u=-b*p/r*x2;
     U2=[U2;u];
     x2=x2+dt*(a*x2+b*u);
     t=t+dt;
 end;
 
 figure(1);
 hold on;
 plot(T,X2,'.-');
 figure(2);
 hold on;
 plot(T,U2,'.-');
