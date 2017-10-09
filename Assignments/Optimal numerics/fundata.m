function[g,dg]=fundata(lambda,Lambda,x_3,X)
eps=0.01; n=3;

 g=(lambda(1)^2+lambda(2)^2+(lambda(3)+2*x_3)^2);
 dg=[];
 for i=2:n+1;
     lambda=Lambda(:,i); x=X(:,i);
     Gi=(lambda(1)^2+lambda(2)^2+(lambda(3)+2*x_3)^2);
     dg=[dg,(Gi-g)/eps];
 end;

