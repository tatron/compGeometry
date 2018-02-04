function [V, D, C, U] = problem2(num_pts, v, numEigs, timeEnd)
%problem2 Summary of this function goes here
%   Detailed explanation goes here

a = 0;
b = 2*pi;
x = linspace(a,b,num_pts+1)';
x = x(1:end-1);
C = [cos(x),sin(x)];
seg = zeros(num_pts,2);
seg(:,1) = (1:num_pts);
seg(1:end-1,2) = (2:num_pts);
seg(end,2) = 1;

vX = v(x);

[Stiff, Mass] = getFEMmats(C, seg);

%Eigenvalue Problem  
[V,D] = eigs(Stiff,Mass,numEigs,'sm');

%Heat Equation
dt = max(x(2:end,1)-x(1:end-1,1))^2;
t = (0:dt:timeEnd);
U = zeros(num_pts,length(t));
U(:,1) = vX(:);

for i = 2:length(t)
    U(:,i) = (Mass + 0.5*(t(i)-t(i-1))*Stiff)\(Mass-0.5*...
        (t(i)-t(i-1))*Stiff)*U(:,i-1);
end    

end