%alpha = @(s) [sin(s)+2*sin(2*s);cos(s)-2*cos(2*s); sin(3*s)];
% Unit Circle
alpha = @(s) [cos(s);sin(s);0*s];
%gradAlpha = @(s) [-sin(s);cos(s);0*s];


a = -pi;
b = pi;
h = 0.1;
%x = 0:.05:2*pi;
x = a:h:b;
%x = [x,b];

numEigs = 8;
num_pt = length(x);
seg = zeros(num_pt,2);
seg(:,1) = 1:num_pt;
seg(1:end-1,2) = 2:num_pt;
seg(end,2) = 1;
num_seg = size(seg,1);

MassMat = sparse(num_pt,num_pt);
StiffMat = sparse(num_pt,num_pt);
for i = 1:num_seg
    p_1 = seg(i,1);
    p_2 = seg(i,2);
    v_1 = alpha(x(p_1));
    v_2 = alpha(x(p_2));
    segArea = norm(v_2-v_1);
    g = segArea^2;
    StiffMat(p_1,p_2) = -1/segArea;
    StiffMat(p_2,p_1) = StiffMat(p_1,p_2);
    StiffMat(p_1,p_1) = 1/segArea + StiffMat(p_1,p_1);
    StiffMat(p_2,p_2) = 1/segArea + StiffMat(p_2,p_2);
    
    MassMat(p_1,p_2) = 1/6*segArea;
    MassMat(p_2,p_1) = MassMat(p_1,p_2);
    MassMat(p_1,p_1) = 1/3*segArea + MassMat(p_1,p_1);
    MassMat(p_2,p_2) = 1/3*segArea + MassMat(p_2,p_2);
end    

[V,D] = eigs(StiffMat,MassMat,numEigs,'sm');
cLim = [min(min(V)),max(max(V))];

trace = alpha(x);
str = 'Eigenfunction %d';
hPlot = gobjects(1,8);
for j = 1:numEigs
    hPlot(j) = subplot(2,4,j);
    set(hPlot(j),'clim',cLim);
    title(sprintf(str,j));
    patch(trace(1,:),trace(2,:),trace(3,:),V(:,j),'FaceColor','none',...
        'EdgeColor','interp');
end    
colormap jet;