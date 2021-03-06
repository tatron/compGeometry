function ViewPC(pt,varargin)

if nargin < 1
    disp(' Parameters: ViewPC(pt,f,varargin) ' );
    return
end
s = 1;

%figure;
switch nargin
    case 1
        plot3(pt(:,1),pt(:,2),pt(:,3),'b.','MarkerSize',3);
    case 2
        f = varargin{1};
        scatter3(pt(:,1),pt(:,2),pt(:,3),s,f,'MarkerSize',5);
    case 3
        f = varargin{1};
        s = varargin{2};
        scatter3(pt(:,1),pt(:,2),pt(:,3),s,f);
end


axis off;
axis equal;