%Project1 Script

numPts = 2.^(3:8);
numPts_sz = length(numPts);

%%  Problem 1
v1 = @(x) exp(-(x-0.5).^2./0.1^2);
v2 = @(theta) exp(-(theta-pi).^2./0.1^2);

numEigs = 8;
timeEnd = 0.15;

V = cell(numPts_sz,2);
D = cell(numPts_sz,2);
c = cell(numPts_sz,2);
U = cell(numPts_sz,2);

D1_true = (pi*(0:numEigs-1)).^2;
D1_true = diag(D1_true);
D1_error = zeros(numEigs,numPts_sz);

D2_true = zeros(numEigs,1);
D2_true(1) = 0;
D2_true(2:2:end) = (1:length(D2_true(2:2:end))).^2;
D2_true(3:2:end) = (1:length(D2_true(3:2:end))).^2;
D2_true = diag(D2_true);
D2_error = zeros(numEigs,numPts_sz);

for i = 1:numPts_sz
    [V{i,1},D{i,1},c{i,1},U{i,1}] = problem1(numPts(i), v1, numEigs, ...
        timeEnd);
    D1_error(:,i) = sqrt(sum((D{i,1}-D1_true).^2))';
    [V{i,2},D{i,2},c{i,2},U{i,2}] = problem2(numPts(i), v2, numEigs, ...
        timeEnd);
    D2_error(:,i) = sqrt(sum((D{i,2}-D2_true).^2))';
end 

%%  Problem 1 Plots

xL = min(min(c{end,1}(:,1)));
xU = max(max(c{end,1}(:,1)));
yL = min(min(V{end,1}));
yU = max(max(V{end,1}));
str = 'Eigenfunction %d';
hPlot1 = gobjects(1,numEigs);
for j = 1:numEigs
    hPlot1(j) = subplot(2,numEigs/2,j);
    plot(c{end,1}(:,1),V{end,1}(:,j));
    axis([xL xU yL yU]);
    title(sprintf(str,j));
end    

figure
loglog(numPts,D1_error);
title('Problem 1; Log-Log Plot of Eigenvalue Error vs. Number of Points');
axis tight;
xlabel('Number of Points in Discretization');
ylabel('Eigenvalue Error');

logX = log(numPts(end)) - log(numPts(end-2));
logY = log(D1_error(end,end)) - log(D1_error(end,end-2));
ratio1 = -logY/logX;

yL = min(min(U{end,1}));
yU = max(max(U{end,1}));
figure
for i= 1:size(U{end,1},2)
    plot(c{end,1},U{end,1}(:,i));
    axis([xL xU yL yU]);
    drawnow;
end    

%%  Problem 2 Plots

xL = min(min(c{end,2}(:,1)));
xU = max(max(c{end,2}(:,1)));
yL = min(min(c{end,2}(:,2)));
yU = max(max(c{end,2}(:,2)));
zL = min(min(V{end,2}));
zU = max(max(V{end,2}));

str = 'Eigenfunction %d';
hPlot2 = gobjects(1,numEigs);
for j = 1:numEigs
    hPlot2(j) = subplot(2,numEigs/2,j);
    plot3([c{end,2}(:,1);c{end,2}(1,1)],[c{end,2}(:,2);c{end,2}(1,2)],...
        [V{end,2}(:,j);V{end,2}(1,j)]);
    title(sprintf(str,j));
    axis([xL xU yL yU zL zU]);
    title(sprintf(str,j));
end    

figure
loglog(numPts,D2_error);
title('Problem 2; Log-Log Plot of Eigenvalue Error vs. Number of Points');
axis tight;
xlabel('Number of Points in Discretization');
ylabel('Eigenvalue Error');

logX = log(numPts(end)) - log(numPts(end-2));
logY = log(D2_error(end,end)) - log(D2_error(end,end-2));
ratio2 = -logY/logX;

zL = min(min(U{end,2}));
zU = max(max(U{end,2}));
figure
for i= 1:size(U{end,2},2)
    plot3([c{end,2}(:,1);c{end,2}(1,1)],[c{end,2}(:,2);c{end,2}(1,2)], ...
        [U{end,2}(:,i);U{end,2}(1,i)]);
    axis([xL xU yL yU zL zU]);
    drawnow;
end

%%  Problem 3

[ptTK, LsegTK, vTK, ptSL, LsegSL, vSL] = LoadCurves;

%%  Problem 3 Continued
[StiffTK, MassTK] = getFEMmats(ptTK, LsegTK);
[StiffSL, MassSL] = getFEMmats(ptSL, LsegSL);

[Vtk, Dtk] = eigs(StiffTK,MassTK,numEigs,'sm');
[Vsl, Dsl] = eigs(StiffSL,MassSL,numEigs,'sm');

hTK = (1/(min(full(diag(StiffTK)))/2))^2;
tTK = (0:hTK:0.4);
Utk = zeros(length(ptTK),length(tTK));
Utk(:,1) = vTK(:);

for i = 2:length(tTK)
    Utk(:,i) = (MassTK + 0.5*(tTK(i)-tTK(i-1))*StiffTK)...
        \(MassTK-0.5*(tTK(i)-tTK(i-1))*StiffTK)*Utk(:,i-1);
end    

hSL = (1/(min(full(diag(StiffSL)))/2))^2;
tSL = (0:0.04:0.4);
Usl = zeros(length(ptSL),length(tSL));
Usl(:,1) = vSL(:);

for i = 2:length(tSL)
    Usl(:,i) = (MassSL + 0.5*(tSL(i)-tSL(i-1))*StiffSL)...
        \(MassSL-0.5*(tSL(i)-tSL(i-1))*StiffSL)*Usl(:,i-1);
end    

%%  Problem 3; Trefoil Knot Plots

cLim = [min(min(Vtk)),max(max(Vtk))];
str = 'Eigenfunction %d';
hPlot3 = gobjects(1,numEigs);
for j = 1:numEigs
    hPlot3(j) = subplot(2,numEigs/2,j);
    set(hPlot3(j),'clim',cLim);
    title(sprintf(str,j));
    patch([ptTK(:,1);ptTK(1,1)],[ptTK(:,2);ptTK(1,2)],...
        [ptTK(:,3);ptTK(1,3)],[Vtk(:,j);Vtk(1,j)],'FaceColor','none',...
        'EdgeColor','interp');
end    
colormap jet;

cLim = [min(min(Utk)),max(max(Utk))];
figure
for i= 1:length(tTK)
    clf
    patch([ptTK(:,1);ptTK(1,1)],[ptTK(:,2);ptTK(1,2)],...
        [ptTK(:,3);ptTK(1,3)],[Utk(:,i);Utk(1,i)],'FaceColor','none',...
        'EdgeColor','interp');
    caxis(cLim);
    colormap jet;
    pause(0.1)
    hold on
end 