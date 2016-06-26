function [ I12ieq_all ] = andrieu_taper( data_x,data_y,data_z,I12ieq,seg_z,seg_y,seg_index,axis_a,axis_e, path_png,seg_fill )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

I12ieq_all = I12ieq;

[seg_a1, seg_e1, ~] = cart2sph(seg_z(:,1), seg_z(:,2),seg_z(:,3));
[seg_a2, seg_e2, ~] = cart2sph(seg_y(:,1), seg_y(:,2),seg_y(:,3));
seg_a1 = rad2deg(seg_a1);
seg_e1 = rad2deg(seg_e1);
seg_a2 = rad2deg(seg_a2);
seg_e2 = rad2deg(seg_e2);
is_lt0 = (seg_a1<0);
seg_a1(is_lt0) = seg_a1(is_lt0) + 360;
is_lt0 = (seg_a2<0);
seg_a2(is_lt0) = seg_a2(is_lt0) + 360;

axis_c = 1:numel(axis_a);
axis_r = 1:numel(axis_e);
seg_col1 = interp1(axis_a,axis_c,seg_a1);
seg_row1 = interp1(axis_e,axis_r,seg_e1);
seg_col2 = interp1(axis_a,axis_c,seg_a2);
seg_row2 = interp1(axis_e,axis_r,seg_e2);

n_sub = size(seg_z,1);
cmap = jet(100)*255;
color = zeros(n_sub,3);
index = seg_fill;
index = index + abs(min(index));
index = index/ max(index);
index = round(index*100);
index(index==0) = 1;
figure('visible','off');
for s = 1:n_sub
    color(s,:) = cmap(index(s),:);
end

for s = 1:n_sub
    sub_x = data_x(seg_index{s});
    sub_y = data_y(seg_index{s});
    sub_z = data_z(seg_index{s});
    
    [sub_a, sub_e,~] = cart2sph(sub_x,sub_y,sub_z);
    sub_a = rad2deg(sub_a);
    sub_e = rad2deg(sub_e);
    is_a_neg = sub_a<0;
    sub_a(is_a_neg)=sub_a(is_a_neg)+360;
    sub_row = round(interp1(axis_e,axis_r,sub_e));
    sub_col = round(interp1(axis_a,axis_c,sub_a));
    is_remove = isnan(sub_row) | isnan(sub_col);
    sub_row = sub_row(~is_remove);
    sub_col = sub_col(~is_remove);
    %is_in_row = (sub_row > 1 & sub_row < n_row);
    %is_in_col = (sub_col > 1 & sub_col < n_col);
    %sub_row = sub_row(is_in_row);
    %sub_col = sub_col(is_in_col);
    for p = 1:numel(sub_row);
        %fprintf('pixel %g\n',p);
        val = I12ieq(sub_row(p),sub_col(p),3);
        I12ieq_all(sub_row(p),sub_col(p),1) = color(s,1);
        I12ieq_all(sub_row(p),sub_col(p),2) = color(s,2);
        I12ieq_all(sub_row(p),sub_col(p),3) = color(s,3);

    end
end

I12ieq_all(:,:,1) = fliplr(I12ieq_all(:,:,1));
I12ieq_all(:,:,2) = fliplr(I12ieq_all(:,:,2));
I12ieq_all(:,:,3) = fliplr(I12ieq_all(:,:,3));
%
figure;
imshow(I12ieq_all)
colormap('gray');
axis off; axis image
hold on
set(gca,'Units','normalized','Position',[0 0 1 1]);
[n_row, n_col,~] = size(I12ieq);
set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
for s = 1:n_sub
    plot([n_col-seg_col1(s),n_col-seg_col2(s)],[seg_row1(s) seg_row2(s)],'-g','linewidth', .75);
end
f = getframe(gcf);
filename_png_andrieu_all = sprintf('%s%s',path_png,'andrieu_all.png');
imwrite(f.cdata,filename_png_andrieu_all);
%}
end

