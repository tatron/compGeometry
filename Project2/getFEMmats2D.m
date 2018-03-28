function [Mass,Stiff] = getFEMmats2D(pt,trg)
%getFEMmats2D Summary of this function goes here
%   Detailed explanation goes here

num_pt  = length(pt);
num_trg = length(trg);

Stiff = sparse(num_pt,num_pt);
Mass = sparse(num_pt,num_pt);

for i=1:num_trg
    p = trg(i,:);
    v = pt(p,:);
    
    Area = 1/2*norm(cross(v(1,:)-v(2,:),v(1,:)-v(3,:)));
    
    theta(1) = atan2(2*Area,dot(v(2,:)-v(1,:),v(3,:)-v(1,:)));
    theta(2) = atan2(2*Area,dot(v(1,:)-v(2,:),v(3,:)-v(2,:)));
    theta(3) = atan2(2*Area,dot(v(1,:)-v(3,:),v(2,:)-v(3,:)));
    
    Stiff(p(1),p(1)) = Stiff(p(1),p(1)) + 1/2*(cot(theta(2)) ...
        + cot(theta(3)));
    Stiff(p(2),p(2)) = Stiff(p(2),p(2)) + 1/2*(cot(theta(1)) ...
        + cot(theta(3)));
    Stiff(p(3),p(3)) = Stiff(p(3),p(3)) + 1/2*(cot(theta(1)) ...
        + cot(theta(2))); 
    
    Stiff(p(1),p(2)) = -1/2*cot(theta(3)) + Stiff(p(1),p(2));
    Stiff(p(2),p(1)) = Stiff(p(1),p(2));
    Stiff(p(1),p(3)) = -1/2*cot(theta(2)) + Stiff(p(1),p(3));
    Stiff(p(3),p(1)) = Stiff(p(1),p(3));
    Stiff(p(2),p(3)) = -1/2*cot(theta(1)) + Stiff(p(2),p(3));
    Stiff(p(3),p(2)) = Stiff(p(2),p(3));
    
    Mass(p(1),p(1)) = Mass(p(1),p(1)) + 1/6*Area;
    Mass(p(2),p(2)) = Mass(p(2),p(2)) + 1/6*Area;
    Mass(p(3),p(3)) = Mass(p(3),p(3)) + 1/6*Area;
    
    Mass(p(1),p(2)) = Mass(p(1),p(2)) + 1/12*Area;
    Mass(p(2),p(1)) = Mass(p(1),p(2));
    Mass(p(1),p(3)) = Mass(p(1),p(3)) + 1/12*Area;
    Mass(p(3),p(1)) = Mass(p(1),p(3));
    Mass(p(2),p(3)) = Mass(p(2),p(3)) + 1/12*Area;
    Mass(p(3),p(2)) = Mass(p(2),p(3));
end    

end