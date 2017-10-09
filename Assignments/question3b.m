
clear all;
clc

% Fixed value numerics
dt=0.01; tf=1; T=[];

%% Optimal computation with alpha = 0.25

X_1 = []; X_2 = [];
J_1 = 0; J_2 = 0;

alpha = 0.25;
tau = (2-alpha)/(3*alpha);

for j=0:100
    t=dt*j;
    T=[T;t];
    if(t<tau)
        x_t = -alpha*t +1;
        X_1 = [X_1; x_t];
        J_1 = J_1 + dt*0.5*(x_t - alpha)^2;
    else
        x_t = alpha*(t-2*tau)+1;
        X_2 = [X_2; x_t];
        J_2 = J_2 + dt*0.5*(x_t - alpha)^2;
    end
    
end

disp(strcat('Minimum cost for alpha = 0.25 is : ', {' '}, num2str(J_1 + J_2)));

X_alpha1 = [X_1;X_2];

%% Optimal computation with alpha = 1.0

X_1 = []; X_2 = [];
J_1 = 0; J_2 = 0;

alpha = 1.0;
tau = (2-alpha)/(3*alpha);

for j=0:100
    t=dt*j;
    %T=[T;t];
    if(t<tau)
        x_t = -alpha*t +1;
        X_1 = [X_1; x_t];
        J_1 = J_1 + dt*0.5*(x_t - alpha)^2;
    else
        x_t = alpha*(t-2*tau)+1;
        X_2 = [X_2; x_t];
        J_2 = J_2 + dt*0.5*(x_t - alpha)^2;
    end
    
end

disp(strcat('Minimum cost for alpha = 1 is : ', {' '}, num2str(J_1 + J_2)));

X_alpha2 = [X_1;X_2];

%% Optimal computation with alpha = 1.75

X_1 = []; X_2 = [];
J_1 = 0; J_2 = 0;

alpha = 1.75;
tau = (2-alpha)/(3*alpha);

for j=0:100
    t=dt*j;
    %T=[T;t];
    if(t<tau)
        x_t = -alpha*t +1;
        X_1 = [X_1; x_t];
        J_1 = J_1 + dt*0.5*(x_t - alpha)^2;
    else
        x_t = alpha*(t-2*tau)+1;
        X_2 = [X_2; x_t];
        J_2 = J_2 + dt*0.5*(x_t - alpha)^2;
    end
    
end

disp(strcat('Minimum cost for alpha = 1.75 is : ', {' '}, num2str(J_1 + J_2)));

X_alpha3 = [X_1;X_2];

%% Graph solution

figure(1);
clf;
plot(T,X_alpha1,'r',T,X_alpha2,'g.',T,X_alpha3,'b--'); hold on;
title('Optimal Switching Time: X(t) vs t'); xlabel('Time, t'); ylabel('State Variable, X(t)'); hold on;
legend('alpha = 0.25','alpha = 1.0','alpha = 1.75'); 
hold off;
