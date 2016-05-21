function [ I_ROI_id,I_ROI_range, I_ROI_is_v, I_SEG_id, I_results_v, I_mask_vo2, pixel_results] = ROI_analysis( ROI,I12r,I12xy,I12z0,I12e,seg_row,seg_col,seg_iter,seg_inlier,seg_r, is_valid, t_angular_resolution )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Preliminary calculations 
[n_col, n_row] = size(I12r);
n_ROI = ROI.NumROI;

%% Truth ROI segment 

%Initialize images to store ROI segment number and range 
I_ROI_id = zeros(n_col,n_row);
I_ROI_range = zeros(n_col,n_row);
I_ROI_is_o = false(n_col,n_row);
I_ROI_is_v = false(n_col,n_row);
% Initialize array to store ROI range 
ROI_range = zeros(n_ROI,1);
for r = 1:n_ROI;
    n_ROI_pts = size(ROI.Points{r},1);
    % Color points in image 
    for p = 1:n_ROI_pts
    I_ROI_id( ROI.Points{r}(p,3), ROI.Points{r}(p,2)) = r;
    end
    % Determine range of ROI 
    ROI_range(r) = median(I12r(I_ROI_id==r));
    I_ROI_range((I_ROI_id==r)) = ROI_range(r);
    if strfind(char(ROI.Name(r)),'o') 
        I_ROI_is_o((I_ROI_id==r)) = true;
    else
        I_ROI_is_v((I_ROI_id==r)) = true;
    end
end
clear p
% Make union of both occluded and visible
I_ROI_is_vo = I_ROI_is_o | I_ROI_is_v;

% Display ROI segment image 
%{
cmap = jet(n_ROI);
cmap(1,:) = [0 0 0];
figure; imagesc(I_ROI_id);
colormap(cmap);
axis image
% Display ROI range image
cmap = jet(n_ROI);
cmap(1,:) = [0 0 0];
figure; imagesc(I_ROI_range);
axis image
% Display ROI is either occluded or visible image
figure; imagesc(I_ROI_is_vo);
axis image
% Display ROI is_visible image
figure; imagesc(I_ROI_is_v);
axis image
% Display ROI is_occluded image
figure; imagesc(I_ROI_is_o);
axis image
%}



%% Detected SEG image 
% For each segment
I_SEG_id = zeros(size(I_ROI_id));
I_SEG_inlier = zeros(size(I_ROI_id));
for s = 1:numel(unique(seg_iter));
    if ~is_valid(s);
        continue
    end
    mask_current_seg = false(n_col,n_row);
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
    end
    I_SEG_id(mask_current_seg) = seg_iter(s);
    I_SEG_inlier(mask_current_seg) = seg_inlier(s);
end

%{
cmap = jet(n_ROI);
cmap(1,:) = [0 0 0];
figure; imagesc(I_SEG_id);
colormap(cmap);
axis image
%}

I_SEG_is = (I_SEG_id>0);

%% ROC Curve

%[x,y] = perfcurve(I_ROI_is_vo(:),I_SEG_inlier(:),true);

%% FA/Hit/CR/Miss

% Compute false alarms, misses, hits, and correct rejections 
I_fa_v = ~I_ROI_is_v & I_SEG_is;
I_fa_vo = ~I_ROI_is_vo & I_SEG_is;
I_miss_v = I_ROI_is_v & ~I_SEG_is;
I_miss_vo = I_ROI_is_vo & ~I_SEG_is;
I_cr_v = ~I_ROI_is_v & ~I_SEG_is;
I_cr_vo = ~I_ROI_is_vo & ~I_SEG_is;
I_hit_v = I_ROI_is_v & I_SEG_is;
I_hit_vo = I_ROI_is_vo & I_SEG_is;

% Results image given both visible and occluded 
I_results_vo = zeros(n_col,n_row);
I_results_vo(I_fa_vo) = 1;
I_results_vo(I_miss_vo) = 2;
I_results_vo(I_cr_vo) = 3;
I_results_vo(I_hit_vo) = 4;

% Results image given only visible 
I_results_v = zeros(n_col,n_row);
I_results_v(I_fa_v) = 1;
I_results_v(I_miss_v) = 2;
I_results_v(I_cr_v) = 3;
I_results_v(I_hit_v) = 4;

% Results image with differentiation between visible and occluded 
I_mask_vo2 = zeros(n_col, n_row);
I_mask_vo2(I_fa_v&~I_ROI_is_vo) = 1; % FA 
I_mask_vo2(I_miss_vo&I_ROI_is_v) = 2; % Visible Miss 
I_mask_vo2(I_miss_vo&I_ROI_is_o) = 3; % Occluded Miss
I_mask_vo2(I_cr_v&~I_ROI_is_vo) = 4; % CR
I_mask_vo2(I_hit_vo&I_ROI_is_v) = 5; % Visible Hit 
I_mask_vo2(I_hit_vo&I_ROI_is_o) = 6; % Occluded Hit 

%{
% Display results 
cmap = zeros(4,3);
cmap(1,:) = [ 1 0 0];
cmap(2,:) = [ 0 0 1];
cmap(3,:) = [ 0 0 0];
cmap(4,:) = [ 0 1 0];
figure; imagesc(I_results_vo); axis image; colormap(cmap);
figure; imagesc(I_results_v); axis image; colormap(cmap);
% Differentiation 
cmap = zeros(6,3);
cmap(1,:) = [1 0 0]; %red 
cmap(2,:) = [0 0 1]; %blue
cmap(3,:) = [0 1 1]; % cyan 
cmap(4,:) = [0 0 0]; %black
cmap(5,:) = [0 1 0]; %green
cmap(6,:) = [1 1 0]; % yellow
figure; imagesc(I_mask_vo2); axis image; colormap(cmap);
%}

%% Classification of laser returns (pixels)
n_v = sum(sum(I_ROI_is_v));
n_o = sum(sum(I_ROI_is_o));
n_vo = sum(sum(I_ROI_is_vo));
n_n = sum(sum(~I_ROI_is_vo));

n.FA = sum(sum(I_mask_vo2==1));
n.miss_v = sum(sum(I_mask_vo2==2));
n.miss_o = sum(sum(I_mask_vo2==3));
n.CR = sum(sum(I_mask_vo2==4));
n.hit_v = sum(sum(I_mask_vo2==5));
n.hit_o = sum(sum(I_mask_vo2==6));
n.hit_vo = n.hit_v + n.hit_o;
n.miss_vo = n.miss_v + n.miss_o;

p.FA = n.FA./n_n; % Probability of FA 
p.miss_o = n.miss_o./n_o; % Probability of miss given occluded 
p.miss_v = n.miss_v./n_v; % Probability of miss given visible
p.miss_vo = n.miss_vo./n_vo; % Probability of miss given visible and occluded
p.CR  = n.CR./n_n; % Probability of Correct Rejection
p.hit_o = n.hit_o./n_o; % Probability of hits given occluded 
p.hit_v = n.hit_v./n_v; % Probability of hits given visible 
p.hit_vo = n.hit_vo./n_vo; % Probability of hits given visible and occluded

%% Pixel Classification as a function of other variables 

[ n_of_r, p_of_r ] = Pixel_Classification( I12r,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'range' );
[ n_of_xy, p_of_xy ] = Pixel_Classification( I12xy,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'xy range' );
[ n_of_z0, p_of_z0 ] = Pixel_Classification( I12z0,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'height above ground' );
[ n_of_e, p_of_e ] = Pixel_Classification( I12e,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'elevation angle' );

%[ p_of_dia ] = Pixel_Classification_seg( I_SEG_id,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,seg_r,'diameter [m]' );

%% Area-based classification 
I12area = (2*I12r*tand(t_angular_resolution/2)).^2;
na_v = sum(sum(I12area(I_ROI_is_v)));
na_o = sum(sum(I12area(I_ROI_is_o)));
na_vo = sum(sum(I12area(I_ROI_is_vo)));
na_n = sum(sum(I12area(~I_ROI_is_vo)));

a.FA = sum(sum(I12area(I_mask_vo2==1)));
a.miss_v = sum(sum(I12area(I_mask_vo2==2)));
a.miss_o = sum(sum(I12area(I_mask_vo2==3)));
a.CR = sum(sum(I12area(I_mask_vo2==4)));
a.hit_v = sum(sum(I12area(I_mask_vo2==5)));
a.hit_o = sum(sum(I12area(I_mask_vo2==6)));
a.hit_vo = a.hit_v + a.hit_o;
a.miss_vo = a.miss_v + a.miss_o;

pa.FA = a.FA./na_n; % Probability of FA 
pa.miss_o = a.miss_o./na_o; % Probability of miss given occluded 
pa.miss_v = a.miss_v./na_v; % Probability of miss given visible
pa.miss_vo = a.miss_vo./na_vo; % Probability of miss given visible and occluded
pa.CR  = a.CR./na_n; % Probability of Correct Rejection
pa.hit_o = a.hit_o./na_o; % Probability of hits given occluded 
pa.hit_v = a.hit_v./na_v; % Probability of hits given visible 
pa.hit_vo = a.hit_vo./na_vo; % Probability of hits given visible and occluded


[a_of_r, pa_of_r]  = Pixel_Classification_Area( I12r,I12area,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'range' );
[a_of_xy, pa_of_xy] = Pixel_Classification_Area( I12xy,I12area,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'xy range' );
[a_of_z0, pa_of_z0] = Pixel_Classification_Area( I12z0,I12area,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'height above ground' );
[a_of_e, pa_of_e] = Pixel_Classification_Area( I12e,I12area,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,'elevation angle' );


pixel_results.n = n;
pixel_results.n_of_r = n_of_r;
pixel_results.n_of_xy = n_of_xy;
pixel_results.n_of_z0 = n_of_z0;
pixel_results.n_of_e = n_of_e;
pixel_results.p = p;
pixel_results.p_of_r = p_of_r;
pixel_results.p_of_xy = p_of_xy;
pixel_results.p_of_z0 = p_of_z0;
pixel_results.p_of_e = p_of_e;
pixel_results.a = a;
pixel_results.a_of_r = a_of_r;
pixel_results.a_of_xy = a_of_xy;
pixel_results.a_of_z0 = a_of_z0;
pixel_results.a_of_e = a_of_e;
pixel_results.pa = pa;
pixel_results.pa_of_r = pa_of_r;
pixel_results.pa_of_xy = pa_of_xy;
pixel_results.pa_of_z0 = pa_of_z0;
pixel_results.pa_of_e = pa_of_e;

end

