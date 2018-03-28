function [v_pt] = v(pt)
%v Summary of this function goes here
%   Detailed explanation goes here

v_pt = sin(acos(pt(3)/norm(pt))); 

end

