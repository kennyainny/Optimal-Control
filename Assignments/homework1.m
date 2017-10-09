tic
clear all

% Armijo parameters
alpha=0.5; beta=0.5;
 
% number of iterations
N=0;
 
% initial conditions
u_armijo= [1,1,1,1,1].';
u_line= [1,1,1,1,1].';

%calculate the min value for convergence
b=[-2;1;-1;3;1];
Q=[5 0 8 -1 -3;0 10 9 7 11;8 9 25 0 6;-1 7 0 19 5;-3 11 6 5 18];
Umin = -inv(Q)*b;
[Gmin,dg_min]=costFunc(Umin);
g=0;

%Armijo and Line Search Vectors
Armijo=[]; Linesearch=[];


while(abs(Gmin-g)>= 0.0001 && N<100);
    
     %Line Search Algorithm
    [g,dg]=costFunc(u_line);
    Linesearch=[Linesearch,g];
    u_line = u_line - gammaMin(u_line) * dg';
    
    
    %Armijo Algorithm
    [g,dg]=costFunc(u_armijo);
    Armijo=[Armijo,g];
    j=1;
    [ga,dga]=costFunc(u_armijo-beta*dg');    
    while (ga-g>=-alpha*beta^j*norm(dg)^2);
           j=j+1;
           [ga,dga]=costFunc(u_armijo-beta^j*dg');
    end;
    
    stepdir=-beta^j*dg'; 
    u_armijo=u_armijo+stepdir;
    N=N+1;
end;

 
hold off;
figure(1);
subplot(2,1,1);
plot(1:N,Linesearch,'b', [1,N+2], [Gmin Gmin],'r-.'); title('Line Search Algorithm'); xlabel('Iterations, n'); ylabel('Cost, g(u)');
%PLOT(X,Y,'y-',X,Y,'go') plot(x,y,'r',x1,y1,'r');axis equal;title('XY position'); xlabel('x'); ylabel('y');
hold on;
subplot(2,1,2);
plot(1:N,Armijo,'b', [1,N+2], [Gmin Gmin],'r-.'); title('Armijo Algorithm'); xlabel('Iterations, n'); ylabel('Cost, g(u)');

toc

