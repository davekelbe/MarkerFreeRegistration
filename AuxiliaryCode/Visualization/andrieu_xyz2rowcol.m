function [ row,col ] = andrieu_xyz2rowcol( axis_a,axis_e, axis_c, axis_r, x,y,z )
%ANDRIEU_XYZ2ROWCOL Converts an x,y,z point to row, col of the Andrieu Img.
%   
%
%
%
%   (C) David Kelbe, Rochester Institute of Technology 

[seg_a1, seg_e1, ~] = cart2sph(x, y, z);
seg_a1 = rad2deg(seg_a1);
seg_e1 = rad2deg(seg_e1);
is_lt0 = (seg_a1<0);
seg_a1(is_lt0) = seg_a1(is_lt0) + 360;

col = interp1(axis_a,axis_c,seg_a1);
row = interp1(axis_e,axis_r,seg_e1);


end

