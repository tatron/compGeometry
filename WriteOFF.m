function WriteOFF(fname,pt,trg)

fid = fopen(fname,'w');
fprintf(fid,'OFF\n');
trg = trg - 1;
num_pt  = size(pt,1);
num_trg = size(trg,1);
fprintf(fid,'%d %d %d \n',[num_pt num_trg 0]);

for i=1:num_pt
    fprintf(fid,'%f ',pt(i,:));
    fprintf(fid,'\n');
end;

for i=1:num_trg
    fprintf(fid,'%d ',3);
    fprintf(fid,'%d ',trg(i,:));
    fprintf(fid,'\n');
end;

fclose(fid);
