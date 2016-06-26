function [ seg_percent_behind,seg_percent_inlier, I12_norm, is_valid ] = andrieu_check_showthrough( data_x,data_y,data_z,I12r,colorLUT,...
    seg_index,seg_iter,axis_a,axis_e,cbar_minmax,...
    seg_z_center,seg_y_center,seg_z_left, seg_y_left, seg_z_right, seg_y_right)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

I12r_all = zeros([size(I12r),3]);
I12r_all(:,:,1)= I12r;
I12r_all(:,:,2)= I12r;
I12r_all(:,:,3)= I12r;
I12_norm = I12r_all;
I12_norm = I12_norm./max(I12_norm(:));
I12_norm1 = I12_norm(:,:,1);
I12_norm2 = I12_norm(:,:,2);
I12_norm3 = I12_norm(:,:,3);

mask_inliers = zeros(size(I12r));
mask_seg = zeros(size(I12r));

axis_c = 1:numel(axis_a);
axis_r = 1:numel(axis_e);
% Center 
% Left side 
[seg_row3, seg_col3] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_z_left(:,1), seg_z_left(:,2),seg_z_left(:,3));
[seg_row4, seg_col4] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_y_left(:,1), seg_y_left(:,2),seg_y_left(:,3));
is_valid2 = (abs(seg_col3-seg_col4)<size(I12r,2)/3);
% Right side 
[seg_row5, seg_col5] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_z_right(:,1), seg_z_right(:,2),seg_z_right(:,3));
[seg_row6, seg_col6] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_y_right(:,1), seg_y_right(:,2),seg_y_right(:,3));
is_valid3 = (abs(seg_col5-seg_col6)<size(I12r,2)/3);

is_valid = ( is_valid2 & is_valid3);
seg_row3 = seg_row3(is_valid);
seg_col3 = seg_col3(is_valid);
seg_row4 = seg_row4(is_valid);
seg_col4 = seg_col4(is_valid);
seg_row5 = seg_row5(is_valid);
seg_col5 = seg_col5(is_valid);
seg_row6 = seg_row6(is_valid);
seg_col6 = seg_col6(is_valid);
colorLUT  = colorLUT(is_valid);
seg_index = seg_index(is_valid);

if  nargin>13;
    n_sub = numel(seg_row3);
else
    n_sub = size(seg_z_center,1);
end
cmap = jet(100);
if ~isempty(cbar_minmax)
    colorLUT(colorLUT<cbar_minmax(1)) = cbar_minmax(1);
    colorLUT(colorLUT>cbar_minmax(2)) = cbar_minmax(2);
end

colorLUT_norm = colorLUT - sign(min(colorLUT))*abs((min(colorLUT)));
colorLUT_norm = colorLUT_norm/ max(colorLUT_norm);
colorLUT_norm = round(colorLUT_norm*100);
colorLUT_norm(colorLUT_norm==0) = 1;
color = uint8(255*cmap(colorLUT_norm,:));

% For each segment 
for s = 1:n_sub
    % Get x,y,z of points in current segment 
    sub_x = data_x(seg_index{s});
    sub_y = data_y(seg_index{s});
    sub_z = data_z(seg_index{s});
    % Convert to a,e
    [sub_a, sub_e,~] = cart2sph(sub_x,sub_y,sub_z);
    sub_a = rad2deg(sub_a);
    sub_e = rad2deg(sub_e);
    is_a_neg = sub_a<0;
    sub_a(is_a_neg)=sub_a(is_a_neg)+360;
    % Find row, col
    sub_row = round(interp1(axis_e,axis_r,sub_e));
    sub_col = round(interp1(axis_a,axis_c,sub_a));
    % Remove NaN values  
    is_remove = isnan(sub_row) | isnan(sub_col);
    sub_row = sub_row(~is_remove);
    sub_col = sub_col(~is_remove);
    % Color each pixel satisfying criterion
    for p = 1:numel(sub_row);
        %fprintf('pixel %g\n',p);
       % val = I12ieq(sub_row(p),sub_col(p),3);
        %I12ieq_all(sub_row(p),sub_col(p),1) = color(s,1);
        %I12ieq_all(sub_row(p),sub_col(p),2) = color(s,2);
        %I12ieq_all(sub_row(p),sub_col(p),3) = color(s,3);
        mask_inliers(sub_row(p),sub_col(p)) = seg_iter(s);
    end
end
%{
cmap = jet(n_sub);
cmap(1,:) = [0 0 0];
figure; imagesc(mask_inliers);
colormap(cmap);
axis image
%}

[n_row, n_col,~] = size(I12r);
for s = 1:n_sub
    x = [seg_col3(s) seg_col4(s) seg_col6(s) seg_col5(s)];
    y = [seg_row3(s) seg_row4(s),seg_row6(s) seg_row5(s)];
    if (max(x) - min(x))> n_col/2;
       x(x<n_col/2)= x(x<n_col/2) + n_col;
       e2 = [n_col n_col; 0 n_row];
       [xi,yi,ii] = polyxpoly([x x(1)],[y y(1)],e2(1,:),e2(2,:));
       %{
       figure;
       hold on
       plot(e2(1,:),e2(2,:),'-r');
       plot([x x(1)],[y y(1)],'-g')
       scatter(xi,yi,'b');
       %}
       ii = sort(ii(:,1));
       if ii(1) == 1 && ii(2) == 2;
           xn1 = [xi(2) xi(1) x(4)];
           yn1 = [yi(2) yi(1) y(4)];
           xn2 = [xi(2) x(1) x(2) x(3) xi(1)]-n_col;
           yn2 = [yi(2) y(1) y(2) y(3) yi(1)];
       end
       if ii(1) == 2 && ii(2) == 4;
           xn1 = [xi(2) xi(1) x(3) x(4)];
           yn1 = [yi(2) yi(1) y(3) y(4)];
           xn2 = [xi(2) x(1) x(2) xi(1)]-n_col;
           yn2 = [yi(2) y(1) y(2) yi(1)];
       end
       if ii(1) == 2 && ii(2) == 3;
           xn1 = [xi(2) x(1) xi(1)];
           yn1 = [yi(2) y(1) yi(1)];
           xn2 = [xi(2) xi(1) x(2) x(3) x(4)];
           yn2 = [yi(2) yi(1) y(2) y(3) y(4)];
       end
       if ii(1) == 3 && ii(2) == 4;
           xn1 = [xi(2) x(2) xi(1)]-n_col;
           yn1 = [yi(2) y(2) yi(1)];
           xn2 = [xi(2) xi(1) x(3) x(4) x(1)];
           yn2 = [yi(2) yi(1) y(3) y(4) y(1)];
       end
       if ii(1) == 1 && ii(2) == 4;
           xn1 = [xi(2) x(2) xi(3)];
           yn1 = [yi(2) y(2) yi(3)];
           xn2 = [xi(2) x(4) x(1) x(2) xi(1)]-n_col;
           yn2 = [yi(2) y(4) y(1) y(2) yi(1)];
       end
       %{
       figure; imagesc(mask_inliers);
       hold on
       plot(x,y,'-r','linewidth',1)
       plot(xn1,yn1,'-r','linewidth',1)  
       plot(xn2,yn2,'-r','linewidth',1)
       %} 
       mask_current_seg = roipoly(mask_seg,xn1,yn1) | roipoly(mask_seg,xn2,yn2);
    else
       mask_current_seg = roipoly(mask_seg,x,y);
    end
    mask_seg(mask_current_seg) = seg_iter(s);
end
%{
cmap = jet(n_sub);
cmap(1,:) = [0 0 0];
figure; imagesc(mask_seg);
colormap(cmap);
axis image
%}

seg_z_range = sqrt(sum((seg_z_center.^2),2));
seg_y_range = sqrt(sum((seg_y_center.^2),2));
seg_min_range = min([seg_y_range seg_z_range],[],2);
seg_percent_behind = zeros(n_sub,1);
seg_percent_inlier = zeros(n_sub,1);

% Initialize masks
is_inlier = false(n_row,n_col);
is_farther = false(n_row,n_col);
is_seg = false(n_row,n_col);

% For each segment 
for s = 1:n_sub;
    % Plot mask images 
    %{
    figure; 
    subplot(2,1,1); imagesc(mask_inliers==s); axis image;
    subplot(2,1,2); imagesc(mask_seg==s); axis image;
    %}
    % Show image for each iteration 
    %{
    I12_temp1 = I12_norm(:,:,1);
    I12_temp2 = I12_norm(:,:,2);
    I12_temp3 = I12_norm(:,:,3);
    is_inlier = (mask_inliers == s);
    is_farther = (I12ieq > seg_min_range(s));
    is_seg = (mask_seg == s);
    I12_temp1(is_seg&is_farther&~is_inlier)= 1;
    I12_temp2(is_seg&is_farther&~is_inlier)= 0;
    I12_temp3(is_seg&is_farther&~is_inlier)= 0;
    I12_temp1(is_seg&~is_farther&~is_inlier)= 0;
    I12_temp2(is_seg&~is_farther&~is_inlier)= 0;
    I12_temp3(is_seg&~is_farther&~is_inlier)= 1;
    I12_temp1(is_inlier)= 0;
    I12_temp2(is_inlier)= 1;
    I12_temp3(is_inlier)= 0;
    I12_temp = zeros(size(I12_norm));
    I12_temp(:,:,1) = I12_temp1;
    I12_temp(:,:,2) = I12_temp2;
    I12_temp(:,:,3) = I12_temp3;
    figure; 
    imshow(I12_temp); 
    axis image;
    %}
    
    % Update mask with latest segment 
    is_inlier_current = (mask_inliers == s);
    is_farther_current = (I12r > seg_min_range(s)) | (I12r==0);
    is_seg_current = (mask_seg == s);
    
    n_inlier = sum(sum(is_inlier_current));
    n_farther = sum(sum(is_seg_current&is_farther_current&~is_inlier_current)); 
    n_seg = sum(sum(is_seg_current));
    seg_percent_behind(s) = n_farther./n_seg;
    seg_percent_inlier(s) = n_inlier./n_seg;
    
    is_inlier = is_inlier_current | is_inlier;
    is_farther = is_farther_current | is_farther;
    is_seg = is_seg_current | is_seg;

end

    I12_norm1(is_seg&is_farther&~is_inlier)= 1;
    I12_norm2(is_seg&is_farther&~is_inlier)= 0;
    I12_norm3(is_seg&is_farther&~is_inlier)= 0;
    I12_norm1(is_seg&~is_farther&~is_inlier)= 0;
    I12_norm2(is_seg&~is_farther&~is_inlier)= 0;
    I12_norm3(is_seg&~is_farther&~is_inlier)= 1;
    I12_norm1(is_inlier)= 0;
    I12_norm2(is_inlier)= 1;
    I12_norm3(is_inlier)= 0;
    
    I12_norm = zeros(size(I12_norm));
    I12_norm(:,:,1) = I12_norm1;
    I12_norm(:,:,2) = I12_norm2;
    I12_norm(:,:,3) = I12_norm3;
 
%}
end

