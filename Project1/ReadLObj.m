function [Lpt,Lseg] = ReadLObj(fname)
% [PtCoord,lineInd] = ReadObjLines(fname)
%
% Read in multiple lines in an obj file


fid = fopen(fname,'r');
shape_type = fscanf(fid,'%s',1);

if shape_type ~= 'L'
    disp('not a valide obj line file');
    return;
end;

a = fscanf(fid,'%f',2);
num_pts = a(2);
Lpt = fscanf(fid,'%f',num_pts*3);
Lpt = reshape(Lpt,3,num_pts);
Lpt =Lpt';

num_seg = fscanf(fid,'%d',1);

Lseg = fscanf(fid,'%f',num_seg*2);
Lseg = reshape(Lseg,2,num_seg);
Lseg =Lseg';


fclose(fid);