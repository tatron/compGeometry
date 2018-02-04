function [V, D, c, U] = problem1(num_pts, v, numEigs, timeEnd)
%problem1 Summary of this function goes here
%   Detailed explanation goes here

a = 0;
b = 1;
c = linspace(a,b,num_pts)';
c = [c,zeros(num_pts,1)];
seg = zeros(num_pts-1,2);
seg(:,1) = (1:num_pts-1);
seg(:,2) = (2:num_pts);
vC = v(c(:,1));

[Stiff, Mass] = getFEMmats(c, seg);

%Eigenvalue Problem  
[V,D] = eigs(Stiff,Mass,numEigs,'sm');

%Heat Equation
dt = max(c(2:end,1)-c(1:end-1,1))^2;
t = (0:dt:timeEnd);
U = zeros(num_pts-2,length(t));
U(:,1) = vC(2:end-1);

MassInt = Mass(2:end-1,2:end-1);
StiffInt = Stiff(2:end-1,2:end-1);

for i = 2:length(t)
    U(:,i) = (MassInt + 0.5*(t(i)-t(i-1))*StiffInt)\(MassInt-0.5*...
        (t(i)-t(i-1))*StiffInt)*U(:,i-1);
end    
U = [zeros(1,length(t));U;zeros(1,length(t))];

end