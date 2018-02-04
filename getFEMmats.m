function [StiffM, MassM] = getFEMmats(C, seg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

num_pts = length(C);
StiffM = sparse(num_pts,num_pts);
MassM = sparse(num_pts,num_pts);
for i = 1:length(seg)
    pt_1 = seg(i,1);
    pt_2 = seg(i,2);
    v_1 = C(pt_1,:);
    v_2 = C(pt_2,:);
    segArea = norm(v_2-v_1);
    g = segArea^2;
    StiffM(pt_1,pt_2) = -1/segArea;
    StiffM(pt_2,pt_1) = StiffM(pt_1,pt_2);
    StiffM(pt_1,pt_1) = 1/segArea + StiffM(pt_1,pt_1);
    StiffM(pt_2,pt_2) = 1/segArea + StiffM(pt_2,pt_2);
    
    MassM(pt_1,pt_2) = 1/6*segArea;
    MassM(pt_2,pt_1) = MassM(pt_1,pt_2);
    MassM(pt_1,pt_1) = 1/3*segArea + MassM(pt_1,pt_1);
    MassM(pt_2,pt_2) = 1/3*segArea + MassM(pt_2,pt_2);
end    

end