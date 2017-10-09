clc;
clear all;

alpha =1;

%% u = +1, k = even
k=0;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

figure(1);
clf;
y = alpha*sin(T-t) - cos(T-t);
plot(t,y); hold on;

k=k+1;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

y = alpha*sin(T-t) - cos(T-t);
plot(t,y); hold on;

k=k+1;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

y = alpha*sin(T-t) - cos(T-t);
plot(t,y); hold on;

k=k+1;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

y = alpha*sin(T-t) - cos(T-t);
plot(t,y); hold on;

line([0 T],[0 0], 'LineStyle','- -'); hold on;
title('(lambda2 vs t [k*pi, (k+1)*pi])  |  (u = +1, alpha = 1)');
legend('t = [0, pi]', 't = [pi, 2pi]','t = [2pi, 3pi]', 't = [3pi, 4pi]'); 
hold off;

%% u = -1
k=0;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

figure(2);
clf;
y = alpha*sin(T-t) + cos(T-t);
plot(t,y); hold on;

k=k+1;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

y = alpha*sin(T-t) + cos(T-t);
plot(t,y); hold on;

k=k+1;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

y = alpha*sin(T-t) + cos(T-t);
plot(t,y); hold on;

k=k+1;
t0 = k*pi; T=(k+1)*pi;
t = t0:0.1:T;

y = alpha*sin(T-t) + cos(T-t);
plot(t,y); hold on;

line([0 T],[0 0], 'LineStyle','- -'); hold on;
title('(lambda2 vs t [k*pi, (k+1)*pi])  |  (u = -1, alpha = 1)');
legend('t = [0, pi]', 't = [pi, 2pi]','t = [2pi, 3pi]', 't = [3pi, 4pi]'); 
hold off;



