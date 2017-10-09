function[g_uk, gamma, grad]=linesearch(Uk);
%initialize variables
Q = [ 5 0 8 -1 -3; 0 10 9 7 11; 8 9 25 0 6; -1 7 0 19 5; -3 11 6 5 18];
B = [-2; 1; -1; 3; 1];

g_uk = 0.5 * transpose(Uk) * Q * Uk + transpose(B)*Uk;

%Umin = - inv(Q)* B;
%compute the grad vector
grad = transpose(Uk)*Q + transpose(B);
grad_transpose = Q*Uk + B;
%compute gamma min
gamma = (transpose(grad)*grad) * (inv(grad*Q*grad_transpose));


