function [ a_blocked, e_blocked, xy_blocked ] = lidar_stage_shadow( data_a, data_e, data_xy )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rows = round(4*data_e);
columns = ceil((1/.255)*data_a);
rows = rows+360;
columns = columns +1;

nr = 180*4;
nc = max(columns);
np = numel(data_xy);

I = zeros(nr,nc);
for p = 1:np;
        I(rows(p),columns(p))= data_xy(p);
end
figure; imagesc(I); axis image; colorbar

zero = (I<0.3);
%figure; imagesc(zero); axis image;

cumzero = cumsum(~zero);
%figure; imagesc(cumzero);
[row,col] = find(cumzero==1);
hold on
plot(col,row,'-r');
e_blocked = (row/4)-90;
a_blocked = .255*(col-1);
xy_blocked = I(row,col);
for i = 1:numel(e_blocked);
    xy_blocked(i) = I(row(i),col(i));
end

end

