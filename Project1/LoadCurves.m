function [pt1 Lseg1 v1 pt2 Lseg2 v2] = LoadCurves



%%%%%%% Load Trefoil Knot
figure;
[pt1 Lseg1] = ReadLObj('TrefoilKnot.Lobj');
ViewLObj(pt1,Lseg1); 
title('TrefoilKnot');
num_pt1 = size(pt1,1);
v1 = zeros(num_pt1,1);
v1([num_pt1/2-3:num_pt1/2 + 3]) = 10;
%%%%%% Load a curve from Radon Transform

M=256;% The size of phantom
P = phantom('Modified Shepp-Logan',M);
h = 1/2;
Proj_angle = [0:h:360-h];
pt2 = radon(P,Proj_angle)';
figure;
imshow(P);
title('Modified Shepp-Logan Phantom');
num_pt2 = size(pt2,1);
Lseg2 = zeros(num_pt2,2);
Lseg2(:,1) = [1:num_pt2];
Lseg2(1:end-1,2) = [2:num_pt2];
Lseg2(end,2) = 1;
v2 = zeros(num_pt2,1);
v2([num_pt2/2-3:num_pt2/2 + 3]) = 10;