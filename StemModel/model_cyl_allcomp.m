function [xyz_alpha, r_alpha, is_inlier_alpha,E_alpha,...
    M, r_kasa, r_pratt,xyz_left, xyz_right,return_str,coverage,t_inlier,...
    status] = model_cyl_allcomp( ...
    all_x,all_y,all_z,...
    t_r_adj,t_theta,t_r1_max,t_r1_min,t_taper_min, ...
    t_taper_max, t_sampling,t_min_density,t_min_obj_size,...
    t_min_pts, t_coverage, block_axes, options_verbose,line_a1norm,...
    axis_a, axis_e, I12ieq)
%MODEL_CYL produces a cylinder model for a point subset
%
% ********** Inputs ***********
%
% all_x         - all points within voxel (x-values)
% all_y         - all points within voxel (y-values)
% all_z         - all points within voxel (z-values)
% t_theta       - maximum user-defined lean angle [degrees]
% t_r1_max      - maximum user_defined radius [m]
% t_r1_min      - minimum user-defined radius [m]
% t_taper_min   - minimum user-defined taper [degrees]
% t_taper_max   - maximum user-defined taper [degrees]
% t_sampling    - Sensor sampling parameter = 2tan(divergence/2) [deg]
% t_min_density - minimum density is sub-voxel
% t_min_obj_size - Sensor minimum detectable object size at nominal range [m]
% t_min_pts     - Minimum number of points at range based on system specs
% options_verbose - Run-time option, if true, images are displayed
% line_a1norm   - Initial vector of tree growth
%
% ********* Outputs ***********
%
% xyz_alpha     -
% r_alpha       -
% is_inlier_alpha -
% E_alpha       -
% M             -
% r_kasa        -
% r_pratt       -
% xyz_left      -
% xyz_right     -
% return_str    -
%
% ********* Example ************
%
%
% Copyright (C) 2013, David Kelbe, Rochester Institute of Technology

%% Initialization

% Set empty output in case of early return
r_alpha = [];
r_kasa = [];
r_pratt = [];
xyz_alpha = [];
is_inlier_alpha = [];
E_alpha = [];
M = [];
xyz_left = [];
xyz_right = [];
coverage = [];
t_inlier = [];
status = [];
% Initialize return_str as false
return_str.empty_histogram = false;
return_str.npts_dense = false;
return_str.empty_alpha_shape = false;
return_str.side_points = false;
return_str.angle = false;
return_str.taper = false;
return_str.radius = false;
return_str.coverage = false;
return_str.success = false;

% Two cases of visualization output for efficient debugging
%{
options_point_subset = false;
options_principal_direction = false;
options_orthogonalization_1 = false;
options_insufficient_density = false;
options_density = false;
options_remove_insufficient_pts = false;
options_best_alpha_shape =false;
options_alpha_shapes = false;
options_side_points = false;
options_RANSAC_axis_fit = false;
options_orthogonalization_2 = false;
options_taper = false;
options_homography = false;
options_circ = false;
options_circzoom = false;
options_inliers = false;
options_coverage = false;
options_theta = false;
options_return = true;
%}

%
options_point_subset = false;
options_principal_direction = false;
options_orthogonalization_1 = false;
options_insufficient_density = false;
options_density = false;
options_remove_insufficient_pts = false;
options_best_alpha_shape =false;
options_alpha_shapes = false;
options_side_points = false;
options_RANSAC_axis_fit = false;
options_orthogonalization_2 = false;
options_homography = false;
options_circzoom = false;
options_inliers = false;
options_return = true;
%}

% Other circle fitting methods
options_kasa = false;
options_pratt = false;
options_hough = false;
options_hough_matlab = false;

% Initial computations
n_all = numel(all_z);
info_point_spacing = norm(mean([all_x,all_y,all_z]))*t_sampling;
all_xyz = [all_x all_y all_z];

% Plot point subset
if options_verbose && options_point_subset;
    h1 = figure('position',[890   651   720   270]);
    scatter3(all_x, all_y, all_z, 10, sqrt(all_x.^2 + all_y.^2 + all_z.^2), 'filled');
    hold on
    %axis(block_axes);
    xlabel('x'); ylabel('y'); zlabel('z');
    daspect([1 1 1 ]);
    hc = colorbar; ylabel(hc, 'Range from Sensor [m]');
    campos([0 0 0]);
    title('Initial Point Subset');
    h2 = figure;
    I12ieq_c = color_andrieu( I12ieq,axis_a, axis_e, all_x, all_y, all_z, true(1,n_all) );
    imshow(I12ieq_c)
    hold on
    colormap('gray');
    axis off; axis image
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(h2,'Units','pixels','Position',[900 305 n_col/2 n_row/2]);
    pause
end
%}

if options_verbose && options_point_subset; delete(h1); delete(h2); end
% Outputs:      options             - user-defined run-type options
%               n_all               - number of points (all)
%               info_point_spacing  - nominal spacing between laser scans
%               all_xyz             - concatenated for computation
%% Determine principle direction, GSM, and Density

% z1 is the vector from the sensor to the mean of the data points
line_z1 = [mean(all_x),mean(all_y),mean(all_z)]';
%line_z1 = [median(all_x),median(all_y),median(all_z)]';

% Wiggle nominal direction to find projection with maximum density
if isempty(line_a1norm);
    % Find grid of angle locations within t_theta
    o = tand(t_theta);
    x_centers = linspace(-o,o,5);
    y_centers = linspace(-o,o,5);
    [x_2d,y_2d] = meshgrid(x_centers,y_centers);
    n_angles = numel(x_2d);
    x_1d = reshape(x_2d, [n_angles,1]);
    y_1d = reshape(y_2d, [n_angles,1]);
    % Update to give circular cone
    is_lt_t_theta  = sqrt(x_1d.^2 + y_1d.^2)<=o;
    x_1d = x_1d(is_lt_t_theta);
    y_1d = y_1d(is_lt_t_theta);
    n_angles = numel(x_1d);
    % Initiate arrays for loop
    new_a1 = [x_1d y_1d ones(n_angles,1)];
    new_a1 = new_a1./repmat(sqrt(sum(new_a1.^2,2)),[1,3]);
    angle_max_count = zeros(n_angles,1);
    angle_M1 = zeros(3,3,n_angles);
    %figure;
    for a = 1:n_angles;
        %GSO
        % First vector is principal growth direction
        q1 = new_a1(a,:)';
        line_y1 = line_z1 + q1;
        % Second vector is from sensor to object
        p2 = (line_z1+line_y1)/2;
        w2 = p2-dot(p2,q1)*q1;
        q2 = w2/norm(w2);
        % Third vector orthogonal
        p3 = [1 1 1]';
        w3 = p3 - dot(p3,q1)*q1 - dot(p3,q2)*q2;
        q3 = w3/norm(w3);
        M1 = [q1'; q2'; q3'];
        all_wuv = M1*all_xyz';
        all_w = all_wuv(1,:)';
        all_u = all_wuv(2,:)';
        all_v = all_wuv(3,:)';
        ctrs{1} = min(all_u)-t_min_obj_size:t_min_obj_size:max(all_u)+t_min_obj_size;
        ctrs{2} = min(all_v)-t_min_obj_size:t_min_obj_size:max(all_v)+t_min_obj_size;
        ctrs_w = min(all_w)-t_min_obj_size:t_min_obj_size:max(all_w)+t_min_obj_size;
        % Histogram
        [X2d_count,~,~] = hist3_kelbe([all_u all_v],ctrs);
        % Expected density
        [X2d_exp] = local_point_density(M1, ctrs{1}, ctrs{2}, ctrs_w);%, all_u, all_v, all_w, sqrt(all_x.^2 + all_y.^2) );
        X2d_exp = fliplr(rot90(X2d_exp,-1));
        X2d_ratio = X2d_count./ X2d_exp;
        % Update
        angle_M1(:,:,a) = M1;
        angle_max_count(a) = max(X2d_ratio(:));
        %{
    %subplot(1,2,1); imagesc(X2d_ratio); axis image
    %subplot(1,2,2); hold on; scatter(x_2d_reshape(r),y_2d_reshape(r),50,max_count(r),'filled');
    %daspect([1 1 1]);
    %pause
        %}
    end
    
    % Get best iteration
    [~,a_best] = max(angle_max_count);
    a1 = new_a1(a_best,:)';
    M1 = angle_M1(:,:,a_best);
    line_y1 = line_z1 + a1;
    
else
    % First vector is principal growth direction
    q1 = line_a1norm;
    line_y1 = line_z1 + q1;
    % Second vector is from sensor to object
    p2 = (line_z1+line_y1)/2;
    w2 = p2-dot(p2,q1)*q1;
    q2 = w2/norm(w2);
    % Third vector orthogonal
    p3 = [1 1 1]';
    w3 = p3 - dot(p3,q1)*q1 - dot(p3,q2)*q2;
    q3 = w3/norm(w3);
    M1 = [q1'; q2'; q3'];
end

% Repeat to get uvw points
all_wuv = M1*all_xyz';
all_w = all_wuv(1,:)';
all_u = all_wuv(2,:)';
all_v = all_wuv(3,:)';
% Histogram
ctrs{1} = min(all_u)-t_min_obj_size:t_min_obj_size:max(all_u)+t_min_obj_size;
ctrs{2} = min(all_v)-t_min_obj_size:t_min_obj_size:max(all_v)+t_min_obj_size;
ctrs_w = min(all_w)-t_min_obj_size:t_min_obj_size:max(all_w)+t_min_obj_size;
[X2d_count,~,bins] = hist3_kelbe([all_u all_v],ctrs);
[ X2d_exp] = local_point_density(M1, ctrs{1}, ctrs{2}, ctrs_w);%, all_u, all_v, all_w, sqrt(all_x.^2 + all_y.^2) );
if isempty(X2d_exp);
    return_str.empty_histogram = true;
    return
end
X2d_exp = fliplr(rot90(X2d_exp,-1));
X2d_ratio = X2d_count./ X2d_exp;

if options_verbose && options_principal_direction
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1); imagesc(X2d_ratio); axis image
    colorbar; caxis([0 1]);
    if isempty(line_a1norm);
        subplot(1,3,2);
        hold on;
        scatter(x_1d,y_1d,50,angle_max_count,'filled');
        scatter(x_1d(a_best),y_1d(a_best),100,angle_max_count(a_best));
        daspect([1 1 1]);
    end
    subplot(1,3,3);
    scatter(all_v,all_w,10,all_u,'filled');
    daspect([1 1 1]);
end

%{
% Check if stem is vertical enough
line_a1 = line_z1 - line_y1;
line_a1norm = line_a1/norm(line_a1);
if line_a1norm(3)<0;
    line_a1norm = -line_a1norm;
end
line_theta = acosd(line_a1norm(3));
if line_theta >t_theta
    return_str = 'angle';
    return
end
%}

clear o x_centers y_centers x_2d y_2d n_angles a x_1d y_1d
clear is_lt_t_theta new_a1 angle_max_count angle_M1
clear q1 p2 w2 q2 p3 w3 q3 all_xyz line_z1 line_y1 line_a1norm
clear X2d_count X2d_exp a_best a1 all_wuv ctrs_w
% Outputs:  all_w                   - Projection of all points onto a1norm
%           all_u                   - Projection of points onto sensor line
%           all_v                   - Projection of points onto 3rd basis
%           M1                      - Transformation from inital G.S.O.
%% Show orthogonalization
if options_verbose && options_orthogonalization_1
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1);
    scatter(all_v,all_w,20,all_u,'filled');
    hold on
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    scatter(all_u,all_w,20,all_v,'filled');
    hold on;
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(all_u,all_v,20,all_w,'filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Orthogonalization based on RANSAC');
    pause
end

if options_verbose && options_orthogonalization_1; delete(h); end

%% Density Estimation
X2d_iscount = X2d_ratio > t_min_density;
X2d_dilate = imdilate(X2d_iscount,strel('square',5));
bin_ind = sub2ind([numel(ctrs{1}), numel(ctrs{2})],  bins(:,1),bins(:,2));
is_dense = X2d_dilate(bin_ind);
all_ratio = X2d_ratio(bin_ind);
%}
%is_dense = true(size(all_u));
dense_u = all_u(is_dense);
dense_v = all_v(is_dense);
dense_w = all_w(is_dense);
%dense_ratio = all_ratio(is_dense);
%all_ix = (1:numel(all_u))';
%dense_ix = all_ix(is_dense);

if numel(dense_u) < t_coverage*t_min_pts;
    return_str.npts_dense = true;
end

if options_verbose && options_insufficient_density && return_str.npts_dense;
    h = figure('position', [872   568   794   408],'color','red');
    subplot(1,2,1);
    scatter3(all_v(is_dense),all_u(is_dense),all_w(is_dense),20,all_ratio(is_dense),'filled');
    hold on
    scatter3(all_v(~is_dense),all_u(~is_dense),all_w(~is_dense),20,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    grid on; xlabel('v'); ylabel('u'); zlabel('w');
    colorbar; caxis([0 1])
    campos([0 0 0]);
    title('Density inliers');
    subplot(1,2,2);
    imagesc(X2d_ratio);
    colorbar
    caxis([0 1])
    axis image
    pause
end
if options_verbose && options_insufficient_density && return_str.npts_dense;
    delete(h);
end
if return_str.npts_dense
    return
end

if options_verbose && options_density
    h = figure('position', [872   568   794   408]);
    subplot(1,2,1);
    scatter3(all_v(is_dense),all_u(is_dense),all_w(is_dense),20,all_ratio(is_dense),'filled');
    hold on
    scatter3(all_v(~is_dense),all_u(~is_dense),all_w(~is_dense),20,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    grid on; xlabel('v'); ylabel('u'); zlabel('w');
    colorbar; caxis([0 1])
    campos([0 0 0]);
    title('Density inliers');
    subplot(1,2,2);
    imagesc(X2d_ratio);
    colorbar
    caxis([0 1])
    axis image
    pause
end

if options_verbose && options_density; delete(h); end
clear ctrs X2d_count bins bin_ind X2d_dilate X2d_iscount
% Outputs:  dense_w                 - Points of sufficient density
%           dense_u                 - Points of sufficient density
%           dense_v                 - Points of sufficient density


clear all_is_line_inlier line_a1 line_theta
%{
% Outputs:  line_z1                 - lower (in z) point on a 3D line
%           line_y1                 - upper (in z) point on a 3D line
%           line_a1norm             - normal vector in growth direction
% %
% %% Gram Schmidt Orthogonalization 1
% %{
% % Initial Gram Schmidt Orthogonalization
% all_xyz = [all_x all_y all_z];
% % First vector is principal growth direction
% q1 = line_a1norm;
% % Second vector is from sensor to object
% p2 = (line_z1+line_y1)/2;
% w2 = p2-dot(p2,q1)*q1;
% q2 = w2/norm(w2);
% % Third vector orthogonal
% p3 = [1 1 1]';
% w3 = p3 - dot(p3,q1)*q1 - dot(p3,q2)*q2;
% q3 = w3/norm(w3);
% all_w = dot(all_xyz,repmat(q1',[n_all,1]),2); % Projection onto a1norm
% all_u = dot(all_xyz,repmat(q2',[n_all,1]),2); % Projection onto second basis vector
% all_v = dot(all_xyz,repmat(q3',[n_all,1]),2); % Projection onto third basis vector
% M1 = [q1'; q2'; q3'];
% %}
% % Show initial orthogonalization
% if options_verbose && options_orthogonalization_1
%     h = figure('position', [872   568   794   408]);
%     subplot(1,3,1);
%     scatter(all_v,all_w,20,all_u,'filled');
%     hold on
%     daspect([1 1 1]);
%     xlabel('v'); ylabel('w');
%     subplot(1,3,2)
%     scatter(all_u,all_w,20,all_v,'filled');
%     hold on;
%     daspect([1 1 1]);
%     xlabel('u'); ylabel('w');
%     subplot(1,3,3)
%     hold on
%     scatter(all_u,all_v,20,all_w,'filled');
%     daspect([1 1 1]);
%     xlabel('u'); ylabel('v');
%     suptitle('Orthogonalization based on RANSAC');
%     pause
% end
%
% clear q1 p2 w2 q2 p3 w3 q3 all_xyz line_z1 line_y1 line_a1norm
% if options_verbose && options_orthogonalization_1; delete(h); end
% % Outputs:  all_w                   - Projection of all points onto a1norm
% %           all_u                   - Projection of points onto sensor line
% %           all_v                   - Projection of points onto 3rd basis
% %           M1                      - Transformation from inital G.S.O.
% %% Density Estimation
% % (could just use alpha shapes)
%
% corners_wuv = M1*corners_xyz;
% info_density_step = t_min_obj_size;%6*info_point_spacing;
% ctrs_u = min(corners_wuv(2,:)):info_density_step:max(corners_wuv(2,:));
% ctrs_v = min(corners_wuv(3,:)):info_density_step:max(corners_wuv(3,:));
% ctrs_w = min(corners_wuv(1,:)):info_density_step:max(corners_wuv(1,:));
%
% %{
% figure;
% scatter3(all_u, all_v, all_w,10,'r','filled');
% hold on
% scatter3(corners_wuv(2,:),corners_wuv(3,:),corners_wuv(1,:),50,'b','filled')
% %}
%
% % Compute histogram and dilate
% [ X2d_exp] = local_point_density(M1, ctrs_u, ctrs_v, ctrs_w );
% if isempty(X2d_exp);
%     return_str = 'empty histogram';
%     return
% end
% X2d_exp = fliplr(rot90(X2d_exp,-1));
%
% % Find maximum azimuth and elevation differences
%
% ctrs_uv{1} = ctrs_u;
% ctrs_uv{2} = ctrs_v;
% [X2d_count,~,bins] = hist3_kelbe([all_u all_v],ctrs_uv);
% %X2d_count = rot90(X2d_count);
%
% %{
% % Compute histogram and dilate
% [ X2d_exp, spacing_u, spacing_v ] = local_point_density( all_u,all_v,all_w, M1, info_point_spacing );
% if isempty(X2d_exp);
%     return
% end
% X2d_exp = rot90(X2d_exp,-1);
%
% % Find maximum azimuth and elevation differences
%
% ctrs{1} = spacing_u;
% ctrs{2} = spacing_v;
% [X2d_count,~,bins] = hist3_kelbe([all_u all_v],ctrs);
% %X2d_count = rot90(X2d_count);
% %}
%
% %n_points_exp_max = 9*t_r/info_point_spacing;
% X2d_ratio = X2d_count./X2d_exp;
% %is_dense= X2d_ratio>t_density_cutoff;
% X2d_iscount = X2d_ratio > .2;
% X2d_dilate = imdilate(X2d_iscount,strel('square',3));
% bin_ind = sub2ind([numel(ctrs_u), numel(ctrs_v)],  bins(:,1),bins(:,2));
% is_dense = X2d_dilate(bin_ind);
% all_ratio = X2d_ratio(bin_ind);
% %}
% %is_dense = true(size(all_u));
% dense_u = all_u(is_dense);
% dense_v = all_v(is_dense);
% dense_w = all_w(is_dense);
% %dense_ratio = all_ratio(is_dense);
% %all_ix = (1:numel(all_u))';
% %dense_ix = all_ix(is_dense);
%
% if numel(dense_u) < 0.9*t_min_pts;
%     return_str = 'number of points - dense';
% end
%
% if options_verbose && options_density
%     if strcmp(return_str,'number of points - dense')
%         h = figure('position', [872   568   794   408],'color','red');
%     else
%         h = figure('position', [872   568   794   408]);
%     end
%     subplot(1,2,1);
%     scatter3(all_u(is_dense),all_v(is_dense),all_w(is_dense),20,all_ratio(is_dense),'filled');
%     hold on
%     scatter3(all_u(~is_dense),all_v(~is_dense),all_w(~is_dense),20,[.5 .5 .5],'filled');
%     %scatter3(all_u(is_dense),all_v(is_dense),all_w(is_dense),20,'b','filled');
%     %hold on
%     %scatter3(all_u(~is_dense),all_v(~is_dense),all_w(~is_dense),20,'r','filled');
%     daspect([1 1 1]);
%     grid on; xlabel('v'); ylabel('w');
%     colorbar; caxis([0 1])
%     campos([0 0 0]);
%     title('Density inliers');
%     subplot(1,2,2);
%     imagesc(X2d_ratio);
%     colorbar
%     caxis([0 1])
%     axis image
%     pause
% end
%
% if options_verbose && options_density; delete(h); end
% if strcmp(return_str,'number of points - dense')
%     return
% end
% clear ctrs X2d_count bins bin_ind X2d_dilate X2d_iscount is_dense
% % Outputs:  dense_w                 - Points of sufficient density
% %           dense_u                 - Points of sufficient density
% %           dense_v                 - Points of sufficient density
% %}
%}
%% Best alpha shape for uv projection
% Alpha shape of dense points in uv projection ('nadir' view)
%if isempty(dense_u)
%    foo = 1;
%    return
%end
info_point_spacing = norm(mean([all_x(is_dense),all_y(is_dense),all_z(is_dense)]))*t_sampling;

[~,Suv_dense] = alphavol([dense_u dense_v],info_point_spacing);

%
if isempty(Suv_dense.bnd)
    return_str.empty_alpha_shape = true;
    return
end
%}
% Find the largest component
[ bnd_uv_x_alpha, bnd_uv_y_alpha, is_alpha ] = alpha_maxcomp( Suv_dense.bnd, dense_u,dense_v );

% Plot results
if options_verbose && options_best_alpha_shape
    h = figure('position',[1031 542 443 412]);
    plot(bnd_uv_x_alpha,bnd_uv_y_alpha,'-r','LineWidth',2)
    hold on
    scatter(dense_u(is_alpha),dense_v(is_alpha),10,'b','filled');
    scatter(dense_u(~is_alpha),dense_v(~is_alpha),10,'r','filled');
    scatter(all_u(~is_dense),all_v(~is_dense),10,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    title('Best alpha shape');
    %pause
end

% Subset points based on alpha shape
alpha_u = dense_u(is_alpha);
alpha_v = dense_v(is_alpha);
alpha_w = dense_w(is_alpha);

if options_verbose && options_best_alpha_shape; delete(h); end
% Outputs:  alpha_w                 - Projection of best alpha points onto a1norm
%           alpha_u                 - Projection of best alpha points onto sensor line
%           alpha_v                 - Projection of best alpha points onto 3rd basis
%           is_alpha                - Logical index to best alpha points
%% Alpha shapes for remaining projections
%{
 % Test example
 %  [va,ua,wa] = sphere;
 %  va = va(:); ua = ua(:); wa = wa(:);
 %  point_spacing = 1;
%}
% Adjust point spacing in case of variablity within block
%info_point_spacing_min = min(sqrt(sum(([all_x,all_y,all_z]).^2,2)))*t_sampling;
info_point_spacing_max = max(sqrt(sum(([all_x,all_y,all_z]).^2,2)))*t_sampling;
%info_point_spacing_mean = mean(sqrt(sum(([all_x,all_y,all_z]).^2,2)))*t_sampling;


% Alpha shapes with "best" points for vw and uw projections
[~,Svw_alpha] = alphavol([alpha_v alpha_w],info_point_spacing_max);
[~,Suw_alpha] = alphavol([alpha_u alpha_w],info_point_spacing_max);

if isempty(Svw_alpha.bnd) || isempty(Suw_alpha.bnd);
    return_str.empty_alpha_shape = true;
    return
end
% Find largest alpha shape for remaining projections
[ ~, ~, ~, bnd_vw_x, bnd_vw_y, ~ ] = alpha_all_comp( Svw_alpha.bnd, alpha_v,alpha_w );
[ ~, ~, ~, bnd_uw_x, bnd_uw_y, ~ ] = alpha_all_comp( Suw_alpha.bnd, alpha_u,alpha_w );
n_gv = numel(bnd_vw_x);
n_gu = numel(bnd_uw_x);
%{
% Plot results
if options_verbose && options_alpha_shapes
    h = figure('position',[865 549 806 405]);
    subplot(1,3,1);
    % plot(Svw.x(Svw.b)',Svw.y(Svw.b)','-r','LineWidth',2)
    plot(bnd_vw_x,bnd_vw_y,'-r','LineWidth',2)
    hold on
    scatter(all_v(is_alpha),all_w(is_alpha),10,'b','filled');
    scatter(all_v(~is_alpha),all_w(~is_alpha),10,'k','filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    plot(bnd_uw_x,bnd_uw_y,'-r','LineWidth',2)
    hold on
    scatter(all_u(is_alpha),all_w(is_alpha),10,'b','filled');
    scatter(all_u(~is_alpha),all_w(~is_alpha),10,'k','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    plot(bnd_uv_x_alpha,bnd_uv_y_alpha,'-r','LineWidth',2)
    hold on
    scatter(all_u(is_alpha),all_v(is_alpha),10,'b','filled');
    scatter(all_u(~is_alpha),all_v(~is_alpha),10,'k','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Alpha shapes');
    %pause
end
%}

% Plot results
if options_verbose && options_alpha_shapes
    h = figure('position',[865 549 806 405]);
    subplot(1,3,1);
    hold on
    % plot(Svw.x(Svw.b)',Svw.y(Svw.b)','-r','LineWidth',2)
    for g = 1:n_gv
        plot(bnd_vw_x{g},bnd_vw_y{g},'-r','LineWidth',2)
    end
    scatter(dense_v(is_alpha),dense_w(is_alpha),10,'b','filled');
    scatter(dense_v(~is_alpha),dense_w(~is_alpha),10,'g','filled');
    scatter(all_v(~is_dense),all_w(~is_dense),10,[.85 .85 .85],'filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    hold on
    for g = 1:n_gu;
        plot(bnd_uw_x{g},bnd_uw_y{g},'-r','LineWidth',2)
    end
    scatter(dense_u(is_alpha),dense_w(is_alpha),10,'b','filled');
    scatter(dense_u(~is_alpha),dense_w(~is_alpha),10,'g','filled');
    scatter(all_u(~is_dense),all_w(~is_dense),10,[.85 .85 .85],'filled');
    daspect([1 1 1])
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    plot(bnd_uv_x_alpha,bnd_uv_y_alpha,'-r','LineWidth',2)
    hold on
    scatter(dense_u(is_alpha),dense_v(is_alpha),10,'b','filled');
    scatter(dense_u(~is_alpha),dense_v(~is_alpha),10,'g','filled');
    scatter(all_u(~is_dense),all_v(~is_dense),10,[.85 .85 .85],'filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Alpha shapes');
    pause
end

clear Suv_alpha Svw_alpha Suw_alpha bnd_uv_x_alpha bnd_uv_y_alpha
if options_verbose && options_alpha_shapes; delete(h); end
% Outputs:  bnd_vw_x                - x Boundary points for vw projection
%           bnd_vw_y                - y Boundary points for vw projection
%           bnd_uw_x                - x Boundary points for uw projection
%           bnd_uw_y                - y Boundary points for uw projection
%% Find side points

f_base = 45; % Fixed value 360/8

side1_vw_x_cell = cell(n_gv,1);
side1_vw_y_cell = cell(n_gv,1);
side2_vw_x_cell = cell(n_gv,1);
side2_vw_y_cell = cell(n_gv,1);
for g = 1:n_gv;
    % Array of adjacent edge points
    edge_vw_x = [bnd_vw_x{g}(1:end-1) circshift(bnd_vw_x{g}(1:end-1),[-1 0])];
    edge_vw_y = [bnd_vw_y{g}(1:end-1) circshift(bnd_vw_y{g}(1:end-1),[-1 0])];
    % Difference of adjacent edge points
    delta_vw_x = edge_vw_x(:,1) - edge_vw_x(:,2);
    delta_vw_y = edge_vw_y(:,1) - edge_vw_y(:,2);
    delta_vw = [delta_vw_x delta_vw_y];
    delta_vw = delta_vw./repmat(sqrt(delta_vw(:,1).^2 + delta_vw(:,2).^2),[1,2]);
    % Angle of edge segments
    theta_vw =  atan2d(delta_vw(:,1),delta_vw(:,2));
    % Logical defining side for edge points
    is_sidevw1 = (theta_vw<f_base&theta_vw>-f_base);
    is_sidevw2 = (theta_vw>180-f_base)|(theta_vw<-180+f_base);
    % Array of side points
    side1_vw_x_cell{g} = edge_vw_x(is_sidevw1)';
    side1_vw_y_cell{g} = edge_vw_y(is_sidevw1)';
    side2_vw_x_cell{g} = edge_vw_x(is_sidevw2)';
    side2_vw_y_cell{g} = edge_vw_y(is_sidevw2)';
end
side1_vw_x = [side1_vw_x_cell{:}]';
side1_vw_y = [side1_vw_y_cell{:}]';
side2_vw_x = [side2_vw_x_cell{:}]';
side2_vw_y = [side2_vw_y_cell{:}]';


side1_uw_x_cell = cell(n_gu,1);
side1_uw_y_cell = cell(n_gu,1);
for g = 1:n_gu;
    % Array of adjacent edge points
    edge_uw_x = [bnd_uw_x{g}(1:end-1) circshift(bnd_uw_x{g}(1:end-1),[-1 0])];
    edge_uw_y = [bnd_uw_y{g}(1:end-1) circshift(bnd_uw_y{g}(1:end-1),[-1 0])];
    % Difference of adjacent edge points
    delta_uw_x = edge_uw_x(:,1) - edge_uw_x(:,2);
    delta_uw_y = edge_uw_y(:,1) - edge_uw_y(:,2);
    delta_uw = [delta_uw_x delta_uw_y];
    delta_uw = delta_uw./repmat(sqrt(delta_uw(:,1).^2 + delta_uw(:,2).^2),[1,2]);
    % Angle of edge segments
    theta_uw = atan2d(delta_uw(:,1),delta_uw(:,2));
    % Logical defining side for edge points
    is_sideuw1 = (theta_uw<f_base&theta_uw>-f_base);
    % Array of side points
    side1_uw_x_cell{g} = edge_uw_x(is_sideuw1)';
    side1_uw_y_cell{g} = edge_uw_y(is_sideuw1)';
end
side1_uw_x = [side1_uw_x_cell{:}]';
side1_uw_y = [side1_uw_y_cell{:}]';

if options_verbose && options_side_points
    h = figure('position',[865 549 806 405]);
    subplot(1,2,1);
    hold on
    scatter(all_v(~is_dense), all_w(~is_dense),10,[.85 .85 .85],'filled')
    scatter(alpha_v, alpha_w,10,[.5 .5 .5],'filled')
    scatter(side1_vw_x, side1_vw_y,20,'r','filled')
    %is_noside = (~is_sidevw1 & ~is_sidevw2);
    %scatter(edge_vw_x(is_noside,1), edge_vw_y(is_noside,1),50,[.5 .5 .5],'filled')
    scatter(side2_vw_x, side2_vw_y,20,'g','filled')
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,2,2);
    hold on
    scatter(all_u(~is_dense), all_w(~is_dense),10,[.85 .85 .85],'filled')
    scatter(alpha_u, alpha_w,10,[.5 .5 .5],'filled')
    scatter(side1_uw_x, side1_uw_y,20,'r','filled')
    %scatter(edge_uw_x(~is_sideuw1,1), edge_uw_y(~is_sideuw1,1),50,[.5 .5 .5],'filled')
    xlabel('u'); ylabel('w');
    daspect([1 1 1]);
    suptitle('Side points');
    pause
end

clear delta_vw_x delta_vw_y delta_vw theta_vw
clear delta_uw_x delta_uw_y delta_uw theta_uw
clear is_sidevw1 is_sidevw2 is_sideuw1
clear edge_vw_x edge_vw_y edge_uw_x edge_uw_y
if options_verbose && options_side_points; delete(h); end
% Outputs:  f_base                  - Fixed value of 45 degrees (360/8)
%           side1_vw_x              - Side points (x) of vw projection
%           side1_vw_y              - Side points (y) of vw projection
%           side2_vw_x              - Side points (x) of vw projection
%           side2_vw_y              - Side points (y) of vw projection
%           side1_uw_x              - Side points (x) of vw projection
%           side1_uw_y              - Side points (y) of vw projection
%% RANSAC axes fit

% Base (a) and top (b) points for side 1 (left) and 2 (right) of vw proj.
if numel(side1_vw_x)<3 || numel(side2_vw_x)<3 || numel(side1_uw_x)<3;
    return_str.side_points = true;
    return;
end
% RANSAC
if options_verbose && options_RANSAC_axis_fit
    [S1_a_vw,S1_b_vw, is_inlier_vw1,E_vw1,statusvw1] = ransac_2line([side1_vw_x side1_vw_y]',info_point_spacing/2);
    [S2_a_vw,S2_b_vw, is_inlier_vw2,E_vw2,statusvw2] = ransac_2line([side2_vw_x side2_vw_y]',info_point_spacing/2);
    [S1_a_uw,S1_b_uw, is_inlier_uw1,E_uw1,statusuw1] = ransac_2line([side1_uw_x side1_uw_y]',info_point_spacing/2);
else
    [S1_a_vw,S1_b_vw, ~,~,statusvw1] = ransac_2line([side1_vw_x side1_vw_y]',info_point_spacing/2);
    [S2_a_vw,S2_b_vw, ~,~,statusvw2] = ransac_2line([side2_vw_x side2_vw_y]',info_point_spacing/2);
    [S1_a_uw,S1_b_uw, ~,~,statusuw1] = ransac_2line([side1_uw_x side1_uw_y]',info_point_spacing/2);
end

status.vw1 = statusvw1;
status.vw2 = statusvw2;
status.uw1 = statusuw1;

% Reformat points into x's and ys'
s1x_vw = [S1_a_vw(1) S1_b_vw(1)]; % x values for side 1
s1y_vw = [S1_a_vw(2) S1_b_vw(2)]; % y values for side 1
s2x_vw = [S2_a_vw(1) S2_b_vw(1)]; % x values for side 2
s2y_vw = [S2_a_vw(2) S2_b_vw(2)]; % y values for side 2
s1x_uw = [S1_a_uw(1) S1_b_uw(1)];
s1y_uw = [S1_a_uw(2) S1_b_uw(2)];

% Average v and w values for edges to get center line
s12x_vw = (s1x_vw + s2x_vw)/2;
s12y_vw = (s1y_vw + s2y_vw)/2;

% Show results of axis fitting
if options_verbose && options_RANSAC_axis_fit
    h = figure('position', [872   568   794   408]);
    subplot(1,2,1);
    hold on
    plot(s12x_vw, s12y_vw,'-y','linewidth',2)
    plot(s1x_vw, s1y_vw,'-b','linewidth',2)
    plot(s2x_vw, s2y_vw,'-g','linewidth',2)
    scatter(all_v(~is_dense), all_w(~is_dense),10,[.85 .85 .85],'filled')
    scatter(alpha_v,alpha_w,20,[.5 .5 .5],'filled');
    %scatter(side1_vw_x(is_inlier_vw1),side1_vw_y(is_inlier_vw1),50,E_vw1(is_inlier_vw1),'filled');
    %scatter(side1_vw_x(~is_inlier_vw1),side1_vw_y(~is_inlier_vw1),50,E_vw1(~is_inlier_vw1));
    scatter(side1_vw_x(is_inlier_vw1),side1_vw_y(is_inlier_vw1),20,'b','filled');
    scatter(side1_vw_x(~is_inlier_vw1),side1_vw_y(~is_inlier_vw1),20,'r','filled');
    scatter(side2_vw_x(is_inlier_vw2),side2_vw_y(is_inlier_vw2),20,'g','filled');
    scatter(side2_vw_x(~is_inlier_vw2),side2_vw_y(~is_inlier_vw2),20,'m','filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,2,2)
    hold on
    plot(s1x_uw, s1y_uw,'-b','linewidth',2);
    scatter(all_u(~is_dense), all_w(~is_dense),10,[.85 .85 .85],'filled')
    scatter(alpha_u,alpha_w,20,[.5 .5 .5],'filled');
    scatter(side1_uw_x(is_inlier_uw1),side1_uw_y(is_inlier_uw1),20,'b','filled');
    scatter(side1_uw_x(~is_inlier_uw1),side1_uw_y(~is_inlier_uw1),20,'r','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    pause
end
%}
clear edge_vw_x edge_vw_y edge_uw_x edge_uw_y
clear S1_a_vw S1_b_vw S2_a_vw S2_b_vw
clear S1_a_uw S1_b_uw
if options_verbose && options_RANSAC_axis_fit; delete(h); end
%           s1x_vw                  - Left side x-growth direction in vw proj
%           s1y_vw                  - Left side y-growth direction in vw proj
%           s2x_vw                  - Right side x-growth direction in vw proj
%           s2y_vw                  - Right side y-growth direction in vw proj
%           s12x_vw                 - Average x-growth direction in vw proj
%           s12y_vw                 - Average y-growth direction in vw proj
%           s1x_uw                  - Front side x-growth direction in uw proj
%           s1y_uw                  - Front side y-growth direction in uw proj
%% Grahm Schmidt Orthogonalization 2

% Preparation for slope computations
dv_vw = s12x_vw(2) - s12x_vw(1);
dw_vw = s12y_vw(2) - s12y_vw(1);
du_uw = s1x_uw(2) - s1x_uw(1);
dw_uw = s1y_uw(2) - s1y_uw(1);

% G.S.M Re-orthogonalizaation
w1_adj = [dv_vw/dw_vw du_uw/dw_uw 1]';
w1_adj = w1_adj*sign(w1_adj(3));
q1 = w1_adj/norm(w1_adj);
p2 = [mean(alpha_v); mean(alpha_u); mean(alpha_w)]; % vector to center of tree points
w2 = p2-dot(p2,q1)*q1;
q2 = w2/norm(w2);
p3 = [1 1 1]';
w3 = p3 - dot(p3,q1)*q1 - dot(p3,q2)*q2;
q3 = w3/norm(w3);
X_adj = [all_w all_u all_v]';
M2 = [0 0 1; 0 1 0; 1 0 0];
M3 = [q1'; q2'; q3'];
orthog_wuv = (M3*M2*X_adj)';
orthog_w = orthog_wuv(:,1);
orthog_u = orthog_wuv(:,2);
orthog_v = orthog_wuv(:,3);

% Alternative, simpler form
M32 = [q3'; q2'; q1'];
X_adj2 = [all_v all_u all_w]';
orthog_vuw = (M32*X_adj2)';
orthog_v2 = orthog_vuw(:,1);
orthog_u2 = orthog_vuw(:,2);
orthog_w2 = orthog_vuw(:,3);

% Show results of reorthogonalization
if options_verbose && options_orthogonalization_2
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1);
    hold on
    scatter(orthog_v,orthog_w,20,'b','filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    hold on
    scatter(orthog_u,orthog_w,20,'b','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(orthog_u,orthog_v,20,'b','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Reorthogonalization');
    pause
end
%}

clear w1_adj q1 p2 w2 q2 p3 w3 q3 X_adj
clear dv_vw dw _vw du_uw dw_uw
if options_verbose && options_orthogonalization_2; delete(h); end
% Outputs:  orthog_w                - Orthogonalized points in w basis
%           orthog_u                - Orthogonalized points in u basis
%           orthog_v                - Orthogonalized points in v basis
%           M2                      - Transformation matrix 2
%           M3                      - Transformation matrix 3
%% Adjust side vectors

wmin = s12y_vw(2);
wmax = s12y_vw(1);
% Find u value at wmin and wmax points
dudw = (s1x_uw(1) - s1x_uw(2))/(s1y_uw(1) - s1y_uw(2));
dvdw12 = (s12x_vw(1) - s12x_vw(2))/(s12y_vw(1) - s12y_vw(2));
dvdw1 = (s1x_vw(1) - s1x_vw(2))/(s1y_vw(1) - s1y_vw(2));
dvdw2 = (s2x_vw(1) - s2x_vw(2))/(s2y_vw(1) - s2y_vw(2));

u0 = -dudw*s1y_uw(1) + s1x_uw(1);
v0_12 = -dvdw12*s12y_vw(1) + s12x_vw(1);
v0_1 = -dvdw1*s1y_vw(1) + s1x_vw(1);
v0_2 = -dvdw2*s2y_vw(1) + s2x_vw(1);

u = @(w) dudw*w + u0;
v12 = @(w) dvdw12*w + v0_12;
v1 = @(w) dvdw1*w + v0_1;
v2 = @(w) dvdw2*w + v0_2;

sides = [v1(wmin) u(wmin) wmin;...
    v1(wmax) u(wmax) wmax;...
    v12(wmin) u(wmin) wmin;...
    v12(wmax) u(wmax) wmax;...
    v2(wmin) u(wmin) wmin;...
    v2(wmax) u(wmax) wmax]';

sides_orthog = M32*sides;
s1x_vw_orthog = [sides_orthog(1,1) sides_orthog(1,2)];
s1y_vw_orthog = [sides_orthog(3,1) sides_orthog(3,2)];
s12x_vw_orthog = [sides_orthog(1,3) sides_orthog(1,4)];
s12y_vw_orthog = [sides_orthog(3,3) sides_orthog(3,4)];
s2x_vw_orthog = [sides_orthog(1,5) sides_orthog(1,6)];
s2y_vw_orthog = [sides_orthog(3,5) sides_orthog(3,6)];
s1x_uw_orthog = [sides_orthog(2,3) sides_orthog(2,4)];
s1y_uw_orthog = [sides_orthog(3,3) sides_orthog(3,4)];
%{
figure;
subplot(1,2,1);
hold on
plot(s12x_vw, s12y_vw,'-y','linewidth',2)
plot(s1x_vw, s1y_vw,'-b','linewidth',2)
plot(s2x_vw, s2y_vw,'-g','linewidth',2)
scatter(alpha_v,alpha_w,20,[.5 .5 .5],'filled');
xlabel('v'); ylabel('w'); daspect([1 1 1]);
subplot(1,2,2);
hold on
plot(s1x_uw, s1y_uw,'-y','linewidth',2)
scatter(alpha_u,alpha_w,20,[.5 .5 .5],'filled');
xlabel('u'); ylabel('w'); daspect([1 1 1]);
%}

if options_verbose && options_homography
    h = figure('position', [872   568   794   408]);
    subplot(1,2,1);
    hold on
    plot(s12x_vw_orthog, s12y_vw_orthog,'-y','linewidth',2)
    plot(s1x_vw_orthog, s1y_vw_orthog,'-b','linewidth',2)
    plot(s2x_vw_orthog, s2y_vw_orthog,'-g','linewidth',2)
    scatter(orthog_v(is_dense),orthog_w(is_dense),20,[.5 .5 .5],'filled');
    scatter(orthog_v(~is_dense),orthog_w(~is_dense),20,[.85 .85 .85],'filled');
    xlabel('v'); ylabel('w'); daspect([1 1 1]);
    subplot(1,2,2);
    hold on
    plot(s1x_uw_orthog, s1y_uw_orthog,'-y','linewidth',2)
    scatter(orthog_u(is_dense),orthog_w(is_dense),20,[.5 .5 .5],'filled');
    scatter(orthog_u(~is_dense),orthog_w(~is_dense),20,[.85 .85 .85],'filled');
    xlabel('u'); ylabel('w'); daspect([1 1 1]);
    suptitle('Transformed sides');
    pause
end

clear   s s_new s12x_vw_orthog s12y_vw_orthog_orthog
if options_verbose && options_homography; delete(h); end
% Outputs:  s1x_vw_orthog                  - Corrected side 1 x values
%           s2x_vw_orthog                  - Corrected side 2 x values
%           s1y_vw_orthog                  - Corrected side 1 y values
%           s2y_vw_orthog                  - Corrected side 2 y values
%% Compute radii of alpha shapes

% y values
y_step = [min([s1y_vw_orthog s2y_vw_orthog]) max([s1y_vw_orthog s2y_vw_orthog])];
% Slopes of left and right lines
m1 = (s1x_vw_orthog(2) - s1x_vw_orthog(1))/(s1y_vw_orthog(2) - s1y_vw_orthog(1));
m2 = (s2x_vw_orthog(2) - s2x_vw_orthog(1))/(s2y_vw_orthog(2) - s2y_vw_orthog(1));
% Intercepts of left and right lines
b1 = (-s1y_vw_orthog(1)*m1) + s1x_vw_orthog(1);
b2 = (-s2y_vw_orthog(1)*m2) + s2x_vw_orthog(1);
% x-values
s1_step = m1*y_step + b1;
s2_step = m2*y_step + b2;

% Compute radius and center
r_alpha = mean(abs([s1_step' s2_step']),2);
w_alpha = y_step';
n_step_alpha = numel(r_alpha);
u_alpha = s1x_uw_orthog(1) + r_alpha;
v_alpha = mean(s1x_vw_orthog)+r_alpha;

%% Compute radii using other methods
%{
% % Kasa method
% if options_kasa
%     [u_kasa,v_kasa,r_kasa,~,~] = ransac_circle([orthog_u orthog_v]',t_kasa);
% end
% % Pratt method
% if options_pratt
%     Par = CircleFitByPratt([orthog_u orthog_v]);
%     u_pratt = Par(1);
%     v_pratt = Par(2);
%     r_pratt = Par(3);
% end
% %{
% if options_hough_matlab;
%     radius_range = [40 60];
%     isize = 100;
%     image = zeros(isize,isize);
%     un_demean = round(isize*(un-min(un))/max((un-min(un))));
%     vn_demean = round(isize*(vn-min(vn))/max((vn-min(vn))));
%     for p = 1:numel(un);
%         if un_demean(p) ~=0 && vn_demean(p) ~=0;
%             image(un_demean(p),vn_demean(p)) = image(un_demean(p),vn_demean(p))+1;
%         end
%     end
%     [centers, radii] = imfindcircles(image,radius_range);
% end
% if options_hough;
%     [u_Hough, v_Hough, r_Hough] = Tree_Detection_fxn(un',vn');
% end
%
% %}
%
% if options_verbose && options_circzoom
%     h = figure('position',[1005 56 557 472]);
%     hold on
%     scatter(orthog_u,orthog_v,20,'b','filled');
%     for ws = 1:n_step_alpha;
%         rectangle('position',[u_alpha(ws)-r_alpha(ws),v_alpha(ws)-r_alpha(ws),r_alpha(ws)*2,r_alpha(ws)*2],...
%             'curvature',[1,1],'linestyle','-','edgecolor','g','linewidth',2);
%     end
%     axis(axis);
%     if options_kasa
%         rectangle('position',[u_kasa-r_kasa,v_kasa-r_kasa,r_kasa*2,r_kasa*2],...
%             'curvature',[1,1],'linestyle','-','edgecolor','black','linewidth',2);
%     end
%     if options_pratt
%         rectangle('position',[u_pratt-r_pratt,v_pratt-r_pratt,r_pratt*2,r_pratt*2],...
%             'curvature',[1,1],'linestyle','-','edgecolor','m','linewidth',2);
%     end
%     daspect([1 1 1]);
%     xlabel('u'); ylabel('v');
%     title('Detail circles');
%     pause
% end
%
% clear Par
% if options_verbose && options_circzoom; delete(h); end
% % Outputs:  r_kasa                  - Circle radii using kasa method
% %           u_kasa                  - Circle center (u) using kasa method
% %           v_kasa                  - Circle cetner (v) using kasa method
% %           r_pratt                 - Circle radii using pratt method
% %           u_pratt                 - Circle center (u) using pratt method
% %           v_pratt                 - Circle cetner (v) using pratt method
%}

%% Transform back to original coordinates and find circle inliers

M = M3*M2*M1;
wuv_alpha = [ w_alpha u_alpha v_alpha]';
wuv_left = [ w_alpha u_alpha v_alpha + r_alpha]';
wuv_right = [ w_alpha u_alpha v_alpha - r_alpha]';

xyz_alpha = M\wuv_alpha;
xyz_left  = M\wuv_left;
xyz_right  = M\wuv_right;
% Outputs:  xyz_alpha               - Circle coordinates from Alpha-shapes
%           M                       - Transformation back to xyz

%%
u_m = mean(u_alpha);
v_m = mean(v_alpha);
r_m = mean(r_alpha);

X2 = [orthog_u orthog_v]';
X2 = X2 - repmat([u_m; v_m],[1, n_all]);
E_alpha =  abs(sqrt(sum(X2.^2,1))-r_m)';
t_inlier = max(min(2*info_point_spacing,r_m),.2*r_m);
is_inlier_alpha = E_alpha<t_inlier;
is_front = (orthog_u<mean(u_alpha));
is_inlier_alpha = is_inlier_alpha & is_front;

z1 = xyz_alpha(:,end);
y1 = xyz_alpha(:,1);
z1_wuv = wuv_alpha(:,end);
y1_wuv = wuv_alpha(:,1);
z1_uvw = z1_wuv([2,3,1]);
y1_uvw = y1_wuv([2,3,1]);



if options_verbose && options_inliers
    figure;
    hold on
    scatter3(orthog_u(is_inlier_alpha), orthog_v(is_inlier_alpha), orthog_w(is_inlier_alpha),10,[.25 .25 .25],'filled');
    scatter3(orthog_u(~is_inlier_alpha), orthog_v(~is_inlier_alpha), orthog_w(~is_inlier_alpha),10,[.85 .85 .85],'filled');
    [~,~,~] = Cone(z1_uvw,y1_uvw,[r_alpha(2) r_alpha(1)],30,'r',0,0);
    [~,~,~] = Cone(z1_uvw,y1_uvw,[r_alpha(2)+t_inlier r_alpha(1)+t_inlier],30,'r',1,0);
    [~,~,~] = Cone(z1_uvw,y1_uvw,[r_alpha(2)-t_inlier r_alpha(1)-t_inlier],30,'g',1,0);
    alpha(.25);
    daspect([1 1 1]); xlabel('x'); ylabel('y'); zlabel('z');
    view(0,90);
end

clear u_m v_m r_m X2
if options_verbose && options_inliers; delete(h); end
% Outputs:  E_alpha                 - Error for all points from circle
%           is_inlier_alpha         - Inliers based on error threshold

%% Andrieu Image

if ~isempty(axis_a);
    axis_c = 1:numel(axis_a);
    axis_r = 1:numel(axis_e);
    [sub_a, sub_e,~] = cart2sph(all_x(is_inlier_alpha),all_y(is_inlier_alpha),all_z(is_inlier_alpha));
    sub_a = rad2deg(sub_a);
    sub_e = rad2deg(sub_e);
    is_a_neg = sub_a<0;
    sub_a(is_a_neg)=sub_a(is_a_neg)+360;
    sub_row= round(interp1(axis_e,axis_r,sub_e));
    sub_col = round(interp1(axis_a,axis_c,sub_a));
    I12ieq_c = I12ieq;
    is_remove = isnan(sub_row) | isnan(sub_col);
    sub_row = sub_row(~is_remove);
    sub_col = sub_col(~is_remove);
    for p = 1:numel(sub_row);
        val = I12ieq(sub_row(p),sub_col(p),3);
        I12ieq_c(sub_row(p),sub_col(p),1) = val/2;
        I12ieq_c(sub_row(p),sub_col(p),2) = val/2;
    end
    I12ieq_c(:,:,1) = fliplr(I12ieq_c(:,:,1));
    I12ieq_c(:,:,2) = fliplr(I12ieq_c(:,:,2));
    I12ieq_c(:,:,3) = fliplr(I12ieq_c(:,:,3));
    axis_c = 1:numel(axis_a);
    axis_r = 1:numel(axis_e);
    
    %points = [seg_z_center seg_y_center seg_z_left seg_y_left seg_z_right seg_y_right]
    xyz = [xyz_alpha xyz_left xyz_right];
    [seg_row, seg_col] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
        xyz(1,:), xyz(2,:),xyz(3,:));
end

%% Return if criterion are not satisfied

% Check taper
taper1 = atand((s1x_vw_orthog(1)-s1x_vw_orthog(2))/(s1y_vw_orthog(1)-s1y_vw_orthog(2)));
taper2 = -atand((s2x_vw_orthog(1)-s2x_vw_orthog(2))/(s2y_vw_orthog(1)-s2y_vw_orthog(2)));
taper = mean([taper1 taper2]);
if taper <= t_taper_min || taper> t_taper_max;
    return_str.taper = true;
end

% Check radius
if max(r_alpha) > t_r1_max || max(r_alpha) < t_r1_min;
    return_str.radius = true;
end

% Check coverage
bole_area = 2*(mean(r_alpha)+t_inlier)*t_r_adj; %block size
%bole_area = 2*mean(r_alpha)*sum(sqrt((z1-y1).^2));
bole_range = (sqrt(z1'*z1)+sqrt(y1'*y1))/2;
pt_spacing = t_sampling*bole_range;
n_pts_exp = bole_area/pt_spacing^2;
coverage = sum(is_inlier_alpha)/n_pts_exp;
if coverage<t_coverage;
    return_str.coverage = true;
end

% Check angle
a1 = xyz_alpha(:,end) - xyz_alpha(:,1);
a1norm = a1/norm(a1);
if a1norm(3)<0;
    a1norm = -a1norm;
end
theta = acosd(a1norm(3));
if theta > t_theta;
    return_str.angle = true;
end
%{
if (return_str.angle || return_str.coverage ||...
            return_str.radius || return_str.taper);
        options_verbose = true;
end
%}
if options_verbose && options_return
    if (return_str.angle || return_str.coverage ||...
            return_str.radius || return_str.taper);
        color = 'r';
        titlestr = 'Failure';
    else
        color = 'b';
        titlestr = 'Pass';
    end
    h1 = figure('position', [ 872   720   402   256]);
    subplot(1,3,1);
    hold on
    v_val = [-r_alpha r_alpha];
    plot(v_val',[w_alpha w_alpha]',color,'Linewidth',2);
    plot(v_val,[w_alpha w_alpha],color,'Linewidth',2);
    scatter(orthog_v(~is_dense),orthog_w(~is_dense),20,[.85 .85 .85],'filled');
    scatter(orthog_v(is_dense),orthog_w(is_dense),20,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    hold on
    scatter(orthog_u(~is_dense),orthog_w(~is_dense),20,[.85 .85 .95],'filled');
    scatter(orthog_u(is_dense),orthog_w(is_dense),20,[.5 .5 .5],'filled');
    u_val = [u_alpha u_alpha] + [-r_alpha r_alpha];
    plot(u_val',[w_alpha w_alpha]',color,'Linewidth',2);
    plot(u_val,[w_alpha w_alpha],color,'Linewidth',2);
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(orthog_u(~is_dense),orthog_v(~is_dense),20,[.85 .85 .85],'filled');
    scatter(orthog_u(is_dense),orthog_v(is_dense),20,[.5 .5 .5],'filled');
    for ws = 1:numel(r_alpha);
        rectangle('position',[u_alpha(ws)-r_alpha(ws),v_alpha(ws)-r_alpha(ws),r_alpha(ws)*2,r_alpha(ws)*2],...
            'curvature',[1,1],'linestyle','-','edgecolor',color,'linewidth',2);
    end
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle(titlestr);
    h2 = figure('position',[1284 720 384 242]);
    hold on;
    scatter3(all_x(1), all_y(1), all_z(1),10,[1 1 1],'filled');
    scatter3(all_x(1), all_y(1), all_z(1),10,[1 1 1],'filled');
    scatter3(all_x(1), all_y(1), all_z(1),10,[1 1 1],'filled');
    scatter3(all_x(1), all_y(1), all_z(1),10,[1 1 1],'filled');
    scatter3(all_x(is_inlier_alpha), all_y(is_inlier_alpha), all_z(is_inlier_alpha),10,[.5 .5 .5],'filled');
    scatter3(all_x(~is_inlier_alpha), all_y(~is_inlier_alpha), all_z(~is_inlier_alpha),10,[.85 .85 .85],'filled');
    %scatter3(all_x(is_dense), all_y(is_dense), all_z(is_dense),10,[.5 .5 .5],'filled');
    %scatter3(all_x(~is_dense), all_y(~is_dense), all_z(~is_dense),10,[.85 .85 .85],'filled');
    [~,~,~] = Cone(z1,y1,[r_alpha(2) r_alpha(1)],30,color,1,0);
    alpha(.5);
    daspect([1 1 1]); xlabel('x'); ylabel('y'); zlabel('z');
    campos([0 0 0])
    h_leg = legend(sprintf('Coverage = %3.1f%%',100*coverage),...
        sprintf('Angle= %3.1f deg',theta),...
        sprintf('Radius= %3.2f m', max(r_alpha)),...
        sprintf('Taper= %3.1f deg',taper));
    legtxt=findobj(h_leg,'type','text');
    if return_str.coverage; set(legtxt(4),'color',color);end
    if return_str.angle; set(legtxt(3),'color',color);end
    if return_str.radius; set(legtxt(2),'color',color);end
    if return_str.taper; set(legtxt(1),'color',color);end
    if ~isempty(axis_a);
    h3 = figure;
    imshow(I12ieq_c)
    [~, n_col,~] = size(I12ieq);
    hold on
    plot([n_col-seg_col(2),n_col-seg_col(1)],[seg_row(2) seg_row(1)],color,'linewidth', .75);
    plot([n_col-seg_col(4),n_col-seg_col(3)],[seg_row(4) seg_row(3)],color,'linewidth', .75);
    plot([n_col-seg_col(6),n_col-seg_col(5)],[seg_row(6) seg_row(5)],color,'linewidth', .75);
    colormap('gray');
    axis off; axis image
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(h3,'Units','pixels','Position',[900 375 n_col/2 n_row/2]);
    end
    %pause
end

if options_verbose && options_return
    delete(h1); delete(h2); 
    if ~isempty(axis_a); delete(h3); end;
end;
if return_str.angle || return_str.taper || ...
        return_str.radius || return_str.coverage;
    return
end
%%
clear y_step s1y_vw_orthog s2y_vw_orthog s1x_vw_orthog s2x_vw_orthog m1 m2 b1 b2 s1_step s2_step n_step_alpha
% Outputs:  r_alpha                 - Circle radii using alpha shapes
%           w_alpha                 - Cirle heights using alpha shapes
%           u_alpha                 - Circle center (u) using alpha shapes
%           v_alpha                 - Circle cetner (v) using alpha shapes
%% Find inliners to model_cyl function


%%
return_str.success = true;
end

