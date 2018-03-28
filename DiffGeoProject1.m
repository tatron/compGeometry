%DiffGeoProject1 
%Date: 01/31/2018 
%This script runs the first project fo my computational differential 
%geomtery reading course. 

Derror = zeros(numEigs,length(numPtsTest));
%%  Problem 1
num_pts = 10;
a = 0;
b = 1;
C = linspace(a,b,num_pts)';
C = [C,zeros(num_pts,1)];
seg = zeros(num_pts-1,2);
seg(:,1) = (1:num_pts-1);
seg(:,2) = (2:num_pts);
v = @(x) exp(-(x-0.5).^2./0.1^2);
vC = v(C(:,1));

[StiffMat, MassMat] = getFEMmats(C, seg);

%%  P1; Eigenvalue Problem  
numEigs = 8;
[V,D] = eigs(StiffMat,MassMat,numEigs,'sm');

str = 'Eigenfunction %d';
hPlot = gobjects(1,numEigs);
for j = 1:numEigs
    hPlot(j) = subplot(2,numEigs/2,j);
    plot(C(:,1),V(:,j));
    title(sprintf(str,j));
end    

Dtrue = (pi*(0:numEigs-1)).^2;
Dtrue = diag(Dtrue);
Derror(:,idx) = sqrt(sum((D-Dtrue).^2))';

%%  P1; Heat Equation
h = 0.0125;
t = (0:h:0.4);
cMat = zeros(num_pts-2,length(t));
cMat(:,1) = vC(2:end-1);

MassM = MassMat(2:end-1,2:end-1);
StiffM = StiffMat(2:end-1,2:end-1);

for i = 2:length(t)
    cMat(:,i) = (MassM + 0.5*(t(i)-t(i-1))*StiffM)\(MassM-0.5*...
        (t(i)-t(i-1))*StiffM)*cMat(:,i-1);
end    
cMat = [zeros(1,length(t));cMat;zeros(1,length(t))];

figure
for i= 1:length(t)
    plot(C,cMat(:,i));
    axis([0 1 -0.1 1.1]);
    drawnow;
end    

%%  Problem 2
num_pts = 250;
a = 0;
b = 2*pi;
x = linspace(a,b,num_pts+1);
x = x(1:end-1);
C = [cos(x);sin(x)]';
seg = zeros(num_pts,2);
seg(:,1) = (1:num_pts);
seg(1:end-1,2) = (2:num_pts);
seg(end,2) = 1;

v = @(theta) exp(-(theta-pi).^2./0.1^2);
vX = v(x);

[StiffMat, MassMat] = getFEMmats(C, seg);

%%  P2; Eigenvalue Problem
numEigs = 8;

[V,D] = eigs(StiffMat,MassMat,numEigs,'sm');

str = 'Eigenfunction %d';
hPlot = gobjects(1,numEigs);
for j = 1:numEigs
    hPlot(j) = subplot(2,numEigs/2,j);
    plot3([C(:,1);C(1,1)],[C(:,2);C(1,2)],[V(:,j);V(1,j)]);
    title(sprintf(str,j));
end    

%%  P2; Heat Equation
h = 0.005;
t = (0:h:1);
cMat = zeros(num_pts,length(t));
cMat(:,1) = vX(:);

for i = 2:length(t)
    cMat(:,i) = (MassMat + 0.5*(t(i)-t(i-1))*StiffMat)\(MassMat-0.5*...
        (t(i)-t(i-1))*StiffMat)*cMat(:,i-1);
end    

figure
for i= 1:length(t)
    plot3([C(:,1);C(1,1)],[C(:,2);C(1,2)],[cMat(:,i);cMat(1,i)]);
    axis([-1.1 1.1 -1.1 1.1 -0.1 1.1]);
    drawnow;
end    

%%  Problem 3  
[pt1, Lseg1, v1, pt2, Lseg2, v2] = LoadCurves;
[StiffMatTK, MassMatTK] = getFEMmats(pt1, Lseg1);
[StiffMatSL, MassMatSL] = getFEMmats(pt2, Lseg2);

%%  P3; Eigenvalue Problem  Trefoil Knot
numEigs = 8;

[V,D] = eigs(StiffMatTK,MassMatTK,numEigs,'sm');

cLim = [min(min(V)),max(max(V))];
str = 'Eigenfunction %d';
hPlot = gobjects(1,numEigs);
for j = 1:numEigs
    hPlot(j) = subplot(2,numEigs/2,j);
    set(hPlot(j),'clim',cLim);
    title(sprintf(str,j));
    patch([pt1(:,1);pt1(1,1)],[pt1(:,2);pt1(1,2)],[pt1(:,3);pt1(1,3)],...
        [V(:,j);V(1,j)],'FaceColor','none','EdgeColor','interp');
end    
colormap jet;

%%  P3; Heat Equation   Trefoil Knot
h = 0.005;
t = (0:h:1);
cMat = zeros(length(pt1),length(t));
cMat(:,1) = v1(:);

for i = 2:length(t)
    cMat(:,i) = (MassMatTK + 0.5*(t(i)-t(i-1))*StiffMatTK)...
        \(MassMatTK-0.5*(t(i)-t(i-1))*StiffMatTK)*cMat(:,i-1);
end    

cLim = [min(min(cMat)),max(max(cMat))];
for i= 1:length(t)
    caxis(cLim);
    colormap jet;
    patch([pt1(:,1);pt1(1,1)],[pt1(:,2);pt1(1,2)],[pt1(:,3);pt1(1,3)],...
        [cMat(:,i);cMat(1,i)],'FaceColor','none','EdgeColor','interp');
    axis([-2.6 2.6 -3.1 3.1 -1.1 1.1]);
    drawnow;
end 

%%  P3; Eigenvalue Problem  Shepp Logan
numEigs = 8;

[V,D] = eigs(StiffMatSL,MassMatSL,numEigs,'sm');

%%  P3; Heat Equation   Shepp Logan
h = 0.005;
t = (0:h:1);
cMat = zeros(length(pt2),length(t));
cMat(:,1) = v2(:);

for i = 2:length(t)
    cMat(:,i) = (MassMatSL + 0.5*(t(i)-t(i-1))*StiffMatSL)...
        \(MassMatSL-0.5*(t(i)-t(i-1))*StiffMatSL)*cMat(:,i-1);
end 