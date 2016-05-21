function [ n_stems, ba, atdbh, match_SEG, match_TREE ] ...
    = forestry_evaluation_tree( seg_z, seg_y, seg_r, tree, I12z0,I12x,I12y,I12z,I12r,...
    I_SEG_id, I_ROI_id,I_ROI_range, ...
    NEON_DBH_10over, NEON_DBH_10under, path_mat, info_site, info_plot )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Initial computations
[n_row, n_col] = size(I12z0);
n_seg = size(seg_z,1);
n_tree = size(tree,1);

% Make SEG_xy_range image
seg_range = sqrt(sum(((seg_z(:,1:2) + seg_y(:,1:2))/2).^2,2));
I_SEG_range = zeros(n_row,n_col);
for s = 1:n_seg;
    I_SEG_range(I_SEG_id==s) = seg_range(s);
end
% Make is breast height image
%figure; imagesc(I_SEG_id); axis image;
%figure; imagesc(I_ROI_id); axis image;
%figure; imagesc(I12z0); axis image
I_is_bh = I12z0>1.2&I12z0<1.4;
%figure; imagesc(I_is_bh); axis image
I_SEG_bh = I_SEG_id;
I_SEG_bh(~I_is_bh) = 0;
I_ROI_bh = I_ROI_id;
I_ROI_bh(~I_is_bh) = 0;
%figure; imagesc(I_SEG_bh); axis image
%figure; imagesc(I_ROI_bh); axis image
%{
% Number of ROI stems within circular radius range
unique_val = unique(I_ROI_bh);
n_stems_ROI = numel(unique_val);
if unique_val(1)==0;
    n_stems_ROI = n_stems_ROI - 1;
end
% Number of SEG stems within circular radius range
unique_val = unique(I_SEG_bh);
n_stems_SEG = numel(unique_val);
if unique_val(1)==0;
    n_stems_SEG = n_stems_SEG - 1;
end
% Output
% n_stems_ROI
% n_stems_SEG
%}
% Outputs
% n_row
% n_col
% n_seg
% I_SEG_bh
% I_ROI_bh
% I_is_bh
%% Breast Height SEG,ROI,NEON arrays

%Find SEG which are at breast height
unique_SEG_id = unique(I_SEG_id(I_is_bh));
if unique_SEG_id(1) == 0;
    unique_SEG_id = unique_SEG_id(2:end);
end
% BUG??
seg_r = seg_r(:,1);
unique_SEG_id = unique_SEG_id(unique_SEG_id<numel(seg_r));
unique_SEG_x = (seg_z(unique_SEG_id,1) + seg_y(unique_SEG_id,1))/2;
unique_SEG_y = (seg_z(unique_SEG_id,2) + seg_y(unique_SEG_id,2))/2;
unique_SEG_y = -unique_SEG_y;
unique_SEG_r = seg_r(unique_SEG_id);
unique_SEG_xy = sqrt(unique_SEG_x.^2 + unique_SEG_y.^2);
n_unique_SEG = numel(unique_SEG_r);

% Find ROI which are at breast height
unique_ROI_id = unique(I_ROI_id(I_is_bh));
if unique_ROI_id(1) == 0;
    unique_ROI_id = unique_ROI_id(2:end);
end
n_unique_ROI = numel(unique_ROI_id);
unique_ROI_x = zeros(n_unique_ROI,1);
unique_ROI_y = zeros(n_unique_ROI,1);
unique_ROI_z = zeros(n_unique_ROI,1);
unique_ROI_r = zeros(n_unique_ROI,1);
for i = 1:n_unique_ROI;
    is_in = (I_ROI_id==unique_ROI_id(i)&I_is_bh);
    unique_ROI_x(i) = nanmedian(I12x(is_in));
    unique_ROI_y(i) = -nanmedian(I12y(is_in));
    unique_ROI_z(i) = nanmedian(I12z(is_in));
    I12xdist = I12r*.0044; %2*tand(.25/2)
    I12xdist(~is_in) = 0;
    xdist = sum(I12xdist,2);
    xdist(xdist==0) = nan;
    unique_ROI_r(i) = nanmedian(xdist)/2;
end
[unique_ROI_theta, unique_ROI_phi,unique_ROI_range] = cart2sph(unique_ROI_x, unique_ROI_y,unique_ROI_z);
unique_ROI_range = unique_ROI_range + unique_ROI_r/2;
[unique_ROI_x, unique_ROI_y, ~] = sph2cart(unique_ROI_theta,unique_ROI_phi,unique_ROI_range);
unique_ROI_xy = sqrt(unique_ROI_x.^2 + unique_ROI_y.^2);

unique_NEON_x = NEON_DBH_10over.x;
unique_NEON_y = NEON_DBH_10over.y;
unique_NEON_xy = sqrt(unique_NEON_x.^2 + unique_NEON_y.^2);
unique_NEON_r = NEON_DBH_10over.dbh/(2*100);
n_unique_NEON = numel(unique_NEON_x);

n_tree = numel(tree);
for t = 1:n_tree;
    zbh = tree(t).loc(3,1)+ 1.3;
    tree(t).dbh = 2*interp1(tree(t).loc(3,:),tree(t).r,zbh);
end

unique_TREE_x = zeros(n_tree,1);
unique_TREE_y = zeros(n_tree,1);
unique_TREE_r = zeros(n_tree,1);
for t = 1:n_tree;
    unique_TREE_x(t) = tree(t).loc(1,1);
    unique_TREE_y(t) = -tree(t).loc(2,1); %>>>>>>>> Note y flip
    unique_TREE_r(t) = tree(t).dbh/2;
end
unique_TREE_xy = sqrt(unique_TREE_x.^2 + unique_TREE_y.^2);


% Outputs
% unique_roi_x,y,z,r
% unique_SEG_x,y,z,r
% unique_NEON_x,y,z,r
% unique_TREE_x,y,z,r

%% Adjust registration

filepath_homography = sprintf('%shomography_%03.0f-%02.0f.mat',path_mat,info_site,info_plot);
if exist(filepath_homography,'file');
    load(filepath_homography);
else
    % Plot DBH maps
    %
    figure('position', [587 95 1026 832]);
    hold on
    n_unique_TREE = numel(unique_TREE_x);
    %for s = 1:n_unique_ROI;
    %        h = filledCircle([unique_ROI_x(s) unique_ROI_y(s)],unique_ROI_r(s),1000,'b');
    %        set(h,'FaceAlpha',.5);
    %end
    %for s = 1:n_unique_SEG;
    %        h = filledCircle([unique_SEG_x(s) unique_SEG_y(s)],unique_SEG_r(s),1000,'r');
    %        set(h,'FaceAlpha',.5);
    %end
    for s = 1:n_unique_NEON
        h = filledCircle([unique_NEON_x(s) unique_NEON_y(s)],unique_NEON_r(s),1000,'g');
        set(h,'FaceAlpha',.5);
    end
    for s = 1:n_unique_TREE;
        h = filledCircle([unique_TREE_x(s) unique_TREE_y(s)],unique_TREE_r(s),1000,'m');
        set(h,'FaceAlpha',.5);
    end
    h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
    axis equal;
    axis([-10 10 -10 10]);
    xlabel('Easting');
    ylabel('Northing');
    grid on
    
    [x,y] = ginput();
    p_neon = [x(1:2:end-1)';y(1:2:end-1)'];%'; ones(1,n_p)]
    p_lidar = [x(2:2:end)';y(2:2:end)'];%; ones(1,n_p)]
    homography = homography_solve(p_lidar, p_neon);
    save(filepath_homography,'homography');
end
y = homography_transform([unique_TREE_x';unique_TREE_y'],homography);
unique_TREE_x = y(1,:)';
unique_TREE_y = y(2,:)';
y = homography_transform([unique_SEG_x';unique_SEG_y'],homography);
unique_SEG_x = y(1,:)';
unique_SEG_y = y(2,:)';
y = homography_transform([unique_ROI_x';unique_ROI_y'],homography);
unique_ROI_x = y(1,:)';
unique_ROI_y = y(2,:)';

%{
figure('position', [587 95 1026 832]);
hold on
n_unique_TREE = numel(unique_TREE_x);
for s = 1:n_unique_ROI;
    h = filledCircle([unique_ROI_x(s) unique_ROI_y(s)],unique_ROI_r(s),1000,'b');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_SEG;
    h = filledCircle([unique_SEG_x(s) unique_SEG_y(s)],unique_SEG_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_NEON
    h = filledCircle([unique_NEON_x(s) unique_NEON_y(s)],unique_NEON_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_TREE;
    h = filledCircle([unique_TREE_x(s) unique_TREE_y(s)],unique_TREE_r(s),1000,'m');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-10 10 -10 10]);
xlabel('Easting');
ylabel('Northing');
grid on
%}

%% Stem Count
r_step = 2;
rmax = max([unique_ROI_xy; unique_SEG_xy; 10]);
r_array = 0:r_step:rmax+r_step;
n_r = numel(r_array)-1;
%{
n_stems.SEG_of_xy = zeros(n_r,1);
n_stems.ROI_of_xy = zeros(n_r,1);
I_is_used_SEG = false(n_row,n_col);
I_is_used_ROI = false(n_row,n_col);
for r = 1:n_r;
    I_is_in_SEG = (I_SEG_range>r_array(r)&I_SEG_range<=r_array(r+1))&(~I_is_used_SEG)&(I_is_bh);
    I_is_in_ROI = (I_ROI_range>r_array(r)&I_ROI_range<=r_array(r+1))&(~I_is_used_ROI)&(I_is_bh);
    SEG_id = unique(I_SEG_id(I_is_in_SEG));
    n_SEG_id = numel(SEG_id);
    ROI_id = unique(I_ROI_id(I_is_in_ROI));
    n_ROI_id = numel(ROI_id);
    % Ensure no duplicate counting: If a part of a segment is used, it
    % is disqualified from future counts
    for i = 1:n_SEG_id;
        I_is_used_SEG(I_SEG_id==SEG_id(i))=true;
    end
    for i = 1:n_ROI_id;
        I_is_used_ROI(I_ROI_id==ROI_id(i))=true;
    end
    n_stems.SEG_of_xy(r) = n_SEG_id;
    n_stems.ROI_of_xy(r) = n_ROI_id;
end
%}

n_stems.NEON10o_of_xy = zeros(n_r,1);
n_stems.TREE_of_xy = zeros(n_r,1);
n_stems.SEG_of_xy = zeros(n_r,1);
n_stems.ROI_of_xy = zeros(n_r,1);
for r = 1:n_r;
    is_in_NEON = unique_NEON_xy>=r_array(r)&unique_NEON_xy<r_array(r+1);
    is_in_TREE = unique_TREE_xy>=r_array(r)&unique_TREE_xy<r_array(r+1);
    is_in_SEG = unique_SEG_xy>=r_array(r)&unique_SEG_xy<r_array(r+1);
    is_in_ROI = unique_ROI_xy>=r_array(r)&unique_ROI_xy<r_array(r+1);
    is_10o_NEON = unique_NEON_r >= .10/2;
    is_10o_TREE = unique_TREE_r >= .10/2;
    is_10o_SEG = unique_SEG_r >= .10/2;
    is_10o_ROI = unique_ROI_r >= .10/2;
    is_site_NEON = (unique_NEON_x >=-10 & unique_NEON_x <=10 & unique_NEON_y >=-10 & unique_NEON_y<=10);
    is_site_TREE = (unique_TREE_x >=-10 & unique_TREE_x <=10 & unique_TREE_y >=-10 & unique_TREE_y<=10);
    is_site_SEG = (unique_SEG_x >=-10 & unique_SEG_x <=10 & unique_SEG_y >=-10 & unique_SEG_y<=10);
    is_site_ROI = (unique_ROI_x >=-10 & unique_ROI_x <=10 & unique_ROI_y >=-10 & unique_ROI_y<=10);
    n_stems.NEON10o_of_xy(r) = sum(is_in_NEON&is_10o_NEON&is_site_NEON);
    n_stems.TREE_of_xy(r) = sum(is_in_TREE&is_10o_TREE&is_site_TREE);
    n_stems.SEG_of_xy(r) = sum(is_in_SEG&is_10o_SEG&is_site_SEG);
    n_stems.ROI_of_xy(r) = sum(is_in_ROI&is_10o_ROI&is_site_ROI);
end

% Output:
% n_stems_SEG_of_xy
% n_stems_ROI_of_xy
% n_stems_over10_NEON_of_xy
% n_stems_NEON_over
% n_stems_NEON_under
% n_stems_NEON_all
%% Stem count within square site
% Number of NEON stems within square site

%{
n_stems.NEON_over = numel(NEON_DBH_10over.x);
n_stems.NEON_under = numel(NEON_DBH_10under.module);
n_stems.NEON_all = n_stems.NEON_over + n_stems.NEON_under;

unique_is_in = unique_SEG_x >=-10 & unique_SEG_x <=10 & unique_SEG_y >=-10 & unique_SEG_y<=10;
n_stems.SEG_over = sum(unique_is_in&(2*unique_SEG_r>.1));
n_stems.SEG_under = sum(unique_is_in&(unique_SEG_r<=.1)&(unique_SEG_r>=.01));
n_stems.SEG_all = n_stems.SEG_over + n_stems.SEG_under;

unique_is_in = unique_ROI_x >=-10 & unique_ROI_x <=10 & unique_ROI_y >=-10 & unique_ROI_y<=10;
n_stems.ROI_over = sum(unique_is_in&(2*unique_ROI_r>.1));
n_stems.ROI_under = sum(unique_is_in&(unique_ROI_r<=.1)&(unique_ROI_r>=.01));
n_stems.ROI_all = n_stems.ROI_over + n_stems.ROI_under;
%}
n_stems.NEON_all = sum(n_stems.NEON10o_of_xy);
n_stems.TREE_all = sum(n_stems.TREE_of_xy);
n_stems.SEG_all = sum(n_stems.SEG_of_xy);
n_stems.ROI_all = sum(n_stems.ROI_of_xy);
% Ouput:
% n_stems.NEON_all
% n_stems.SEG_all
% n_stems.ROI_all

%% Basal Area
%{
r_step = 2;
rmax = max([unique_ROI_xy; unique_SEG_xy; 10]);
r_array = 0:r_step:rmax+r_step;
n_r = numel(r_array)-1;
ba.SEG_of_xy = zeros(n_r,1);
ba.ROI_of_xy = zeros(n_r,1);
I_is_used_SEG = false(n_row,n_col);
I_is_used_ROI = false(n_row,n_col);
for r = 1:n_r;
    I_is_in_SEG = (I_SEG_range>r_array(r)&I_SEG_range<=r_array(r+1))&(~I_is_used_SEG)&(I_is_bh);
    I_is_in_ROI = (I_ROI_range>r_array(r)&I_ROI_range<=r_array(r+1))&(~I_is_used_ROI)&(I_is_bh);
    SEG_id = unique(I_SEG_id(I_is_in_SEG));
    n_SEG_id = numel(SEG_id);
    ROI_id = unique(I_ROI_id(I_is_in_ROI));
    n_ROI_id = numel(ROI_id);
    % Ensure no duplicate counting: If a part of a segment is used, it
    % is disqualified from future counts
    for i = 1:n_SEG_id;
        I_is_used_SEG(I_SEG_id==SEG_id(i))=true;
    end
    for i = 1:n_ROI_id;
        I_is_used_ROI(I_ROI_id==ROI_id(i))=true;
    end
    radii_SEG = unique_SEG_r(ismember(unique_SEG_id,SEG_id));
    radii_ROI = unique_ROI_r(ismember(unique_ROI_id,ROI_id));
    ba.SEG_of_xy(r) = sum(pi*radii_SEG.^2);
    ba.ROI_of_xy(r) = sum(pi*radii_ROI.^2);
end
ba.AREA_of_xy = pi*r_array.^2;

ba.NEON10o_of_xy = zeros(n_r,1);
for r = 1:n_r;
    is_in = unique_NEON_xy>=r_array(r)&unique_NEON_xy>r_array(r+1);
    ba.NEON10o_of_xy(r) = sum(is_in);
end
%}

ba.NEON10o_of_xy = zeros(n_r,1);
ba.TREE_of_xy = zeros(n_r,1);
ba.SEG_of_xy = zeros(n_r,1);
ba.ROI_of_xy = zeros(n_r,1);
for r = 1:n_r;
    is_in_NEON = unique_NEON_xy>=r_array(r)&unique_NEON_xy<r_array(r+1);
    is_in_TREE = unique_TREE_xy>=r_array(r)&unique_TREE_xy<r_array(r+1);
    is_in_SEG = unique_SEG_xy>=r_array(r)&unique_SEG_xy<r_array(r+1);
    is_in_ROI = unique_ROI_xy>=r_array(r)&unique_ROI_xy<r_array(r+1);
    is_10o_NEON = unique_NEON_r >= .10/2;
    is_10o_TREE = unique_TREE_r >= .10/2;
    is_10o_SEG = unique_SEG_r >= .10/2;
    is_10o_ROI = unique_ROI_r >= .10/2;
    is_site_NEON = unique_NEON_x >=-10 & unique_NEON_x <=10 & unique_NEON_y >=-10 & unique_NEON_y<=10;
    is_site_TREE = unique_TREE_x >=-10 & unique_TREE_x <=10 & unique_TREE_y >=-10 & unique_TREE_y<=10;
    is_site_SEG = unique_SEG_x >=-10 & unique_SEG_x <=10 & unique_SEG_y >=-10 & unique_SEG_y<=10;
    is_site_ROI = unique_ROI_x >=-10 & unique_ROI_x <=10 & unique_ROI_y >=-10 & unique_ROI_y<=10;
    ba.NEON10o_of_xy(r) = sum(pi*unique_NEON_r(is_in_NEON&is_10o_NEON&is_site_NEON).^2);
    ba.TREE_of_xy(r) = sum(pi*unique_TREE_r(is_in_TREE&is_10o_TREE&is_site_TREE).^2);
    ba.SEG_of_xy(r) = sum(pi*unique_SEG_r(is_in_SEG&is_10o_SEG&is_site_SEG).^2);
    ba.ROI_of_xy(r) = sum(pi*unique_ROI_r(is_in_ROI&is_10o_ROI&is_site_ROI).^2);
end

ba.a_of_xy = pi* (r_array(2:end).^2 - (r_array(1:end-1)).^2);
ba.a_of_xy(6:end) = 0;

% Output:
% ba_SEG_of_xy
% ba_ROI_of_xy
% ba_AREA_of_xy
% ba_NEON10o_of_xy
%% Basal Area Square Site
% basal area within square site

%{
ba.NEON_over = sum(pi*(NEON_DBH_10over.dbh/(2*100)).^2);
ba.NEON_under = sum(pi*(NEON_DBH_10under.dbh/(2*100)).^2);
ba.NEON_all = ba.NEON_over + ba.NEON_under; % m^2
ba.AREA = 100; % 10 x 10 m

unique_is_over = unique_SEG_x >=-10 & unique_SEG_x <=10 ...
    & unique_SEG_y >=-10 & unique_SEG_y<=10 & 2*unique_SEG_r>.1;
unique_is_under = unique_SEG_x >=-10 & unique_SEG_x <=10 ...
    & unique_SEG_y >=-10 & unique_SEG_y<=10 & 2*unique_SEG_r<=.1;
ba.SEG_over = sum(pi*(unique_SEG_r(unique_is_over).^2));
ba.SEG_under = sum(pi*(unique_SEG_r(unique_is_under).^2));
ba.SEG_all = ba.SEG_over + ba.SEG_under;

unique_is_over = unique_ROI_x >=-10 & unique_ROI_x <=10 ...
    & unique_ROI_y >=-10 & unique_ROI_y<=10 & 2*unique_ROI_r>.1;
unique_is_under = unique_ROI_x >=-10 & unique_ROI_x <=10 ...
    & unique_ROI_y >=-10 & unique_ROI_y<=10 & 2*unique_ROI_r<=.1;
ba.ROI_over = sum(pi*(unique_ROI_r(unique_is_over).^2));
ba.ROI_under = sum(pi*(unique_ROI_r(unique_is_under).^2));
ba.ROI_all = ba.ROI_over + ba.ROI_under;
%}

ba.NEON_all = sum(ba.NEON10o_of_xy);
ba.TREE_all = sum(ba.TREE_of_xy);
ba.SEG_all = sum(ba.SEG_of_xy);
ba.ROI_all = sum(ba.ROI_of_xy);
ba.a_all = 10.^2;

% Ouput:
% ba.NEON_over/under/all
% ba.SEG_over/under/all
% ba.ROI_over/under/all
% ba.AREA
%% DBH Maps
%{
figure;
hold on
%titlestr = sprintf('DBH map - A%02.0f',info_site);
title('NEON DBH Map');
for s = 1:n_stems_NEON_over
    rectangle('Position',[NEON_DBH_10over.x(s)-NEON_DBH_10over.dbhm(s)/2 ...
        NEON_DBH_10over.y(s)-NEON_DBH_10over.dbhm(s)/2 ...
        NEON_DBH_10over.dbhm(s) NEON_DBH_10over.dbhm(s)],...
        'Curvature',[1 1]);
end
axis equal;
axis([-12 12 -12 12]);
xlabel('Easting');
ylabel('Northing');
grid on

figure;
hold on
title('ROI DBH Map');
for s = 1:n_unique_ROI;
        rectangle('Position',[unique_roi_x(s)-unique_roi_r(s)/2 ...
        unique_roi_y(s)-unique_roi_r(s)/2 ...
        unique_roi_r(s) unique_roi_r(s)],...
        'Curvature',[1 1]);
end
    axis equal;
axis([-12 12 -12 12]);
xlabel('Easting');
ylabel('Northing');
grid on

figure;
hold on
title('SEG DBH Map');
for s = 1:n_unique_SEG;
        rectangle('Position',[unique_SEG_x(s)-unique_SEG_r(s) ...
        unique_SEG_y(s)-unique_SEG_r(s) ...
        unique_SEG_r(s)*2 unique_SEG_r(s)*2],...
        'Curvature',[1 1]);
end
    axis equal;
axis([-12 12 -12 12]);
xlabel('Easting');
ylabel('Northing');
grid on
%}

atdbh.ROI_x = unique_ROI_x;
atdbh.ROI_y = unique_ROI_y;
atdbh.ROI_r = unique_ROI_r;
atdbh.SEG_x = unique_SEG_x;
atdbh.SEG_y = unique_SEG_y;
atdbh.SEG_r = unique_SEG_r;
atdbh.NEON_x = unique_NEON_x;
atdbh.NEON_y = unique_NEON_y;
atdbh.NEON_r = unique_NEON_r;
atdbh.TREE_x = unique_TREE_x;
atdbh.TREE_y = unique_TREE_y;
atdbh.TREE_r = unique_TREE_r;

% ouptut
% unique.
foo = 1;
%% Diameter comparison
delta_xy = .75;
match_SEG  = match_stem( unique_SEG_x, unique_SEG_y, unique_SEG_r, unique_SEG_id,...
    unique_ROI_x, unique_ROI_y, unique_ROI_r, unique_ROI_id,...
    unique_NEON_x, unique_NEON_y, unique_NEON_r, ...
    delta_xy );
%%
%{
figure('position', [587 95 1026 832]);
hold on
for s = 1:n_unique_ROI;
    h = filledCircle([unique_ROI_x(s) unique_ROI_y(s)],unique_ROI_r(s),1000,'b');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_SEG;
    h = filledCircle([unique_SEG_x(s) unique_SEG_y(s)],unique_SEG_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on

hold on
for m = 1:numel(match_SEG.ROI_roi_r);
    plot([match_SEG.ROI_roi_x(m) match_SEG.ROI_seg_x(m)],...
        [match_SEG.ROI_roi_y(m) match_SEG.ROI_seg_y(m)],'-k');
end
%}
%%
%
figure('position', [587 95 1026 832]);
hold on
for s = 1:n_unique_SEG;
    h = filledCircle([unique_SEG_x(s) unique_SEG_y(s)],unique_SEG_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_NEON
    h = filledCircle([unique_NEON_x(s) unique_NEON_y(s)],unique_NEON_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on

hold on
for m = 1:numel(match_SEG.NEON_neon_r);
    plot([match_SEG.NEON_neon_x(m) match_SEG.NEON_seg_x(m)],...
        [match_SEG.NEON_neon_y(m) match_SEG.NEON_seg_y(m)],'-k');
end
%}
%%

unique_TREE_id = (1:n_tree)';
match_TREE  = match_stem( unique_TREE_x, unique_TREE_y, unique_TREE_r, unique_TREE_id,...
    unique_ROI_x, unique_ROI_y, unique_ROI_r, unique_ROI_id,...
    unique_NEON_x, unique_NEON_y, unique_NEON_r, ...
    delta_xy );
match_TREE.ROI_tree_id = match_TREE.ROI_seg_id;
match_TREE.ROI_tree_r = match_TREE.ROI_seg_r;
match_TREE.ROI_tree_xy = match_TREE.ROI_seg_xy;
match_TREE.ROI_tree_x = match_TREE.ROI_seg_x;
match_TREE.ROI_tree_y = match_TREE.ROI_seg_y;
match_TREE.NEON_tree_id = match_TREE.NEON_seg_id;
match_TREE.NEON_tree_r = match_TREE.NEON_seg_r;
match_TREE.NEON_tree_xy = match_TREE.NEON_seg_xy;
match_TREE.NEON_tree_x = match_TREE.NEON_seg_x;
match_TREE.NEON_tree_y = match_TREE.NEON_seg_y;
match_TREE = rmfield(match_TREE,'ROI_seg_id');
match_TREE = rmfield(match_TREE,'ROI_seg_r');
match_TREE = rmfield(match_TREE,'ROI_seg_xy');
match_TREE = rmfield(match_TREE,'ROI_seg_x');
match_TREE = rmfield(match_TREE,'ROI_seg_y');
match_TREE = rmfield(match_TREE,'NEON_seg_id');
match_TREE = rmfield(match_TREE,'NEON_seg_r');
match_TREE = rmfield(match_TREE,'NEON_seg_xy');
match_TREE = rmfield(match_TREE,'NEON_seg_x');
match_TREE = rmfield(match_TREE,'NEON_seg_y');

%%
%{
figure('position', [587 95 1026 832]);
hold on
for s = 1:n_unique_ROI;
    h = filledCircle([unique_ROI_x(s) unique_ROI_y(s)],unique_ROI_r(s),1000,'b');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_tree;
    h = filledCircle([unique_TREE_x(s) unique_TREE_y(s)],unique_TREE_r(s),1000,'m');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on

hold on
for m = 1:numel(match_TREE.ROI_roi_r);
    plot([match_TREE.ROI_roi_x(m) match_TREE.ROI_tree_x(m)],...
        [match_TREE.ROI_roi_y(m) match_TREE.ROI_tree_y(m)],'-k');
end
%}
%%
%{
figure('position', [587 95 1026 832]);
hold on
for s = 1:n_tree;
    h = filledCircle([unique_TREE_x(s) unique_TREE_y(s)],unique_TREE_r(s),1000,'m');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_NEON
    h = filledCircle([unique_NEON_x(s) unique_NEON_y(s)],unique_NEON_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on

hold on
for m = 1:numel(match_TREE.NEON_neon_r);
    plot([match_TREE.NEON_neon_x(m) match_TREE.NEON_tree_x(m)],...
        [match_TREE.NEON_neon_y(m) match_TREE.NEON_tree_y(m)],'-k');
end
%}


% Plot DBH maps
%
figure('position', [587 95 1026 832]);
hold on
n_unique_TREE = numel(unique_TREE_x);
for s = 1:n_unique_ROI;
    h = filledCircle([unique_ROI_x(s) unique_ROI_y(s)],unique_ROI_r(s),1000,'b');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_SEG;
    h = filledCircle([unique_SEG_x(s) unique_SEG_y(s)],unique_SEG_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_NEON
    h = filledCircle([unique_NEON_x(s) unique_NEON_y(s)],unique_NEON_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_unique_TREE;
    h = filledCircle([unique_TREE_x(s) unique_TREE_y(s)],unique_TREE_r(s),1000,'m');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on
for m = 1:numel(match_SEG.ROI_roi_r);
    plot([match_SEG.ROI_roi_x(m) match_SEG.ROI_seg_x(m)],...
        [match_SEG.ROI_roi_y(m) match_SEG.ROI_seg_y(m)],'-k');
end
for m = 1:numel(match_SEG.NEON_neon_r);
    plot([match_SEG.NEON_neon_x(m) match_SEG.NEON_seg_x(m)],...
        [match_SEG.NEON_neon_y(m) match_SEG.NEON_seg_y(m)],'-k');
end
for m = 1:numel(match_TREE.ROI_roi_r);
    plot([match_TREE.ROI_roi_x(m) match_TREE.ROI_tree_x(m)],...
        [match_TREE.ROI_roi_y(m) match_TREE.ROI_tree_y(m)],'-k');
end
for m = 1:numel(match_TREE.NEON_neon_r);
    plot([match_TREE.NEON_neon_x(m) match_TREE.NEON_tree_x(m)],...
        [match_TREE.NEON_neon_y(m) match_TREE.NEON_tree_y(m)],'-k');
end

% Outputs
% match.ROI_roi_id
% match.ROI_roi_r
% match.ROI_seg_id
% match.ROI_seg_r
% match.NEON_neon_id
% match.NEON_neon_r
% match.NEON_neon_id
% match.NEON_neon_r


end

