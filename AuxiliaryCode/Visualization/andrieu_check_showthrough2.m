function [ seg_percent_behind,seg_percent_inlier, I12_norm ] = andrieu_check_showthrough2( data_x,data_y,data_z,I12r,colorLUT,...
    seg_index,seg_iter,axis_a,axis_e,cbar_minmax,...
    seg_col, seg_row,seg_z_center, seg_y_center)
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

n_seg = numel(seg_index);
% For each segment
for s = 1:numel(unique(seg_iter));
    ix_seg = find(seg_iter==s);
    if isempty(ix_seg);
        continue
    end
    for i = 1:numel(ix_seg);
        % Get x,y,z of points in current segment
        sub_x = data_x(seg_index{ix_seg(i)});
        sub_y = data_y(seg_index{ix_seg(i)});
        sub_z = data_z(seg_index{ix_seg(i)});
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
end
%{
cmap = jet(n_seg);
cmap(1,:) = [0 0 0];
figure; imagesc(mask_inliers);
colormap(cmap);
axis image
%}

[n_row, n_col,~] = size(I12r);
for s = 1:numel(unique(seg_iter));
    mask_current_seg = false(n_row,n_col);
    ix_seg = find(seg_iter==s);
    if isempty(ix_seg);
        continue
    end
    if numel(ix_seg) > 1;
        %{
        figure
        hold on
        plot(seg_col{ix_seg(1)},seg_row{ix_seg(1)},'r')
        plot([seg_col{ix_seg(2)}+n_col],seg_row{ix_seg(2)},'g')
        foo = 1;
        %}
    end
    for i = 1:numel(ix_seg);
        x = seg_col{ix_seg(i)};
        y = seg_row{ix_seg(i)};
        if any(~isfinite(x)) || any(~isfinite(y));
            foo = 1;
        end
        %
        mask_current_seg = mask_current_seg | ...
            roipoly(mask_current_seg,x,y);
        foo = 1;
    end
    mask_seg(mask_current_seg) = seg_iter(s);
end
%{
cmap = jet(n_seg);
cmap(1,:) = [0 0 0];
figure; imagesc(mask_seg);
colormap(cmap);
axis image
%}

seg_z_range = sqrt(sum((seg_z_center.^2),2));
seg_y_range = sqrt(sum((seg_y_center.^2),2));
seg_min_range = max([seg_y_range seg_z_range],[],2);
seg_percent_behind = zeros(n_seg,1);
seg_percent_inlier = zeros(n_seg,1);

% Initialize masks
is_inlier = false(n_row,n_col);
is_farther = false(n_row,n_col);
is_seg = false(n_row,n_col);

% For each segment
for s = 1:numel(unique(seg_iter));
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
    is_farther = (I12r > seg_min_range(s));
    is_seg = (mask_seg == s);
    I12_temp1(is_seg&is_farther&~is_inlier)= 1; % Red - farther
    I12_temp2(is_seg&is_farther&~is_inlier)= 0;
    I12_temp3(is_seg&is_farther&~is_inlier)= 0;
    I12_temp1(is_seg&~is_farther&~is_inlier)= 0; % Blue - closer
    I12_temp2(is_seg&~is_farther&~is_inlier)= 0;
    I12_temp3(is_seg&~is_farther&~is_inlier)= 1;
    I12_temp1(is_inlier)= 0; % green - inlier
    I12_temp2(is_inlier)= 1;
    I12_temp3(is_inlier)= 0;
    I12_temp = zeros(size(I12_norm));
    I12_temp(:,:,1) = I12_temp1;
    I12_temp(:,:,2) = I12_temp2;
    I12_temp(:,:,3) = I12_temp3;
    %
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

