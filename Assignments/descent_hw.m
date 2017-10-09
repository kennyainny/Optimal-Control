clear all;

% Armijo parameters
alpha=0.5; beta=0.5;

% number of iterations
N = 1;

% initial conditions
u_armijo = [1; 1; 1; 1; 1];
u_ls = [1; 1; 1; 1; 1];

G_amijo=[]; G_ls=[]; 
difference=[]; 
diff_amijo = -1; diff_ls = -1;
convergence = 0.001;

while( (diff_amijo<0) || (diff_amijo>convergence) || (diff_ls>convergence) && N<20);

    % get the function values
    [g_ls, gamma_ls, grad_ls] = linesearch(u_ls);    
    [g_armijo,dg_armijo]=armijo(u_armijo);

    G_ls=[G_ls,g_ls];
    G_amijo=[G_amijo,g_armijo];

    u_ls = u_ls - gamma_ls * transpose(grad_ls);
    
    j=1;
        [ga,dga]=armijo(u_armijo-beta*dg_armijo');    
        while ((ga-g_armijo>=-alpha*beta^j*norm(dg_armijo)^2) && j <200)
            j=j+1;            
            [ga,dga]=armijo(u_armijo-beta^j*dg_armijo');
        end;
        
    stepdir=-beta^j*dg_armijo;
    u_armijo=u_armijo+stepdir';
    
    len = size(G_ls);

    if(len>1)
    diff_amijo = abs(G_amijo(len) - G_amijo(len-1));
    diff_ls = abs(G_ls(len) - G_ls(len-1));
    end
    
    N = N+1;
end;

hold off;
figure(1);
subplot(2,1,1);
plot(1:N,G_ls);
subplot(2,1,2);
plot(1:N,G_amijo);

%hold off;
%figure(2);
%T=[]; Gg=[];
%for t=-2:0.1:2;
%    T=[T;t];
%    [g,dg]=fundata(t);
%    Gg=[Gg,g];
%end;
%plot(T,Gg);
%hold on;
%plot(U,G,'o');

