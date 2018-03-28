%%Project 2

%%  Problem 1

el = 3;
numEigs = 16;
trueEigs = zeros(numEigs,1);
idx = 1;
for i=0:el
    trueEigs(idx:(idx+2*i)) = i*(i+1);
    idx = idx + 2*i + 1;
end    

pt   = cell(5,1);
trg  = cell(5,1);

Mass  = cell(5,1);
Stiff = cell(5,1);

EigFunc = cell(5,1);
EigValD = cell(5,1);
EigVal  = zeros(5,numEigs);

[pt{1},trg{1}] = ReadOFF('sphere1k.off');
[pt{2},trg{2}] = ReadOFF('sphere2k.off');
[pt{3},trg{3}] = ReadOFF('sphere4k.off');
[pt{4},trg{4}] = ReadOFF('sphere8k.off');
[pt{5},trg{5}] = ReadOFF('sphere16k.off');

for i=1:5
    [Mass{i},Stiff{i}] = getFEMmats2D(pt{i},trg{i});
    [EigFunc{i},EigValD{i}] = eigs(Stiff{i},Mass{i},numEigs,'sm');
    EigVal(i,:) = diag(EigValD{i});
end    

EigValError = EigVal - trueEigs';

plot(log2(abs(EigValError(:,2:end))))
xlabel('Eigenvalue')
ylabel('log_2(Error)')
axis tight
hold on
plot(log2(1./(1000*2.^[0:4])),'r+-')
title('Eigenvalue Error for 1000*2^k Point Spheres')

%%  Problem 2

[pt,trg] = ReadOFF('kitten.off');
[Mass,Stiff] = getFEMmats2D(pt,trg);
[EigFunc,EigVal] = eigs(Stiff,Mass,numEigs,'sm');

%Heat Equation
timeEnd = 0.5;
dt = 0.01;
t = (0:dt:timeEnd);
U = zeros(length(pt),length(t));
for i=1:length(pt)
    U(i,1) = v(pt(i,:));
end    

for i = 2:length(t)
    U(:,i) = (Mass + 0.5*(t(i)-t(i-1))*Stiff)\(Mass-0.5*...
        (t(i)-t(i-1))*Stiff)*U(:,i-1);
end   

%%  Problem 3

[pt,trg] = ReadOFF('face.off');
TR = triangulation(trg,pt);
fe = freeBoundary(TR);
[~,Stiff_M] = getFEMmats2D(pt,trg);
[~,Mass_bM] = getFEMmats(pt,fe);

[V,D] = eigs(Stiff_M,Mass_bM,16,'sm');

%Ask how to make movies
%Ask how to change the patch color