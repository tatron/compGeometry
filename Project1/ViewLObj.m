function ViewLObj(varargin)

color ='r';
if ischar(varargin{1})
    %fpath=varargin{1};
    fname=varargin{1};
    [Lpt,Lseg]= ReadLObj(fname);
    if nargin==2
        color = varargin{2};
    end
    
else Lpt=varargin{1};
    Lseg=varargin{2};
    if nargin==3
        color = varargin{3};
    end
end


num_seg = size(Lseg,1);

for i=1:num_seg
    plot3(Lpt(Lseg(i,:),1),Lpt(Lseg(i,:),2),Lpt(Lseg(i,:),3),color);
    hold on;
end

axis equal

axis off

end