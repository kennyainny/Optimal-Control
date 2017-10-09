function[c] = optunose(u)
global x lambda;
c=(x*u^2*(3+u^2)/((1+u^2)^2)-lambda)^2;
end

