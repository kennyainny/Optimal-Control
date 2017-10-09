clear all;
clc
clf

x=-1; y=0; r=1;
th = 0:0.1:4*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit)
hold on


x=1; y=0; r=1;
th = 0:0.1:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit)
hold on

title('(X1 vs X2)'); xlabel('X1(t)'); ylabel('X2(t)');
legend('u=-1 => (x1+1)^2 + x2^2 = 1', 'u=+1 => (1-x1)^2 + x2^2 = 1'); 
hold on

ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
hold off

