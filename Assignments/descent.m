clear all
% constant step length
gamma=0.02;
gamma=0.5;

% Armijo parameters
alpha=0.8; beta=0.6;

% number of iterations
N=40;

% initial conditions
%u=0.4;
u=-0.4;

% k=1 - fixed step length
% k=2 - gamma(k)=gamma/k
% k=3 - Armijo
k=3;

U=[]; G=[];
for i=1:N;
    % get the function values
    [g,dg]=fundata(u);
    U=[U,u]; G=[G,g];

    if (k==1); %% Fixed step-size gradient descent
        stepdir=-gamma*dg';

    elseif (k==2); %% gamma(i)=gamma/i;
        stepdir=-gamma/i*dg';
        
    else;  %% Amrijo step-size gradient descent
        j=1;
        [ga,dga]=fundata(u-beta*dg');    
        while (ga-g>=-alpha*beta^j*norm(dg)^2);
            j=j+1;
            [ga,dga]=fundata(u-beta^j*dg');
        end;
        stepdir=-beta^j*dg';
    end;

    u=u+stepdir;

end;

hold off;
figure(1);
subplot(2,1,1);
plot(1:N,G);
subplot(2,1,2);
plot(1:N,U);

hold off;
figure(2);
T=[]; Gg=[];
for t=-2:0.1:2;
    T=[T;t];
    [g,dg]=fundata(t);
    Gg=[Gg,g];
end;
plot(T,Gg);
hold on;
plot(U,G,'o');