function [xyz_alpha, r_alpha, is_inlier_alpha,E_alpha,M, r_kasa, r_pratt,xyz_left, xyz_right ] = model_cyl( ...
    all_x,all_y,all_z,corners_xyz,t_r,...
    t_theta,t_sampling,t_density_cutoff,t_min_obj_size,t_min_pts,t_kasa )
%MODEL_CYL produces a cylinder model for a point subset 
%   

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

% Set options
options_verbose = true;

%
options_point_subset = false;
options_principal_direction = false;
options_orthogonalization_1 = false;
options_density = false;
options_remove_insufficient_pts = false;
options_best_alpha_shape =false;
options_alpha_shapes = true;
options_side_points = false;
options_RANSAC_axis_fit = false;
options_orthogonalization_2 = false;
options_homography = false;
options_circ = false;
options_circzoom = false;
options_inliers = false;
%}

%{
options_point_subset = true;
options_principal_direction = true;
options_orthogonalization_1 = true;
options_density = true;
options_best_alpha_shape =true;
options_alpha_shapes = true;
options_side_points = true;
options_RANSAC_axis_fit = true;
options_orthogonalization_2 = true;
options_homography = true;
options_circ = true;
options_circzoom = false;
options_inliers = true;
%}

options_kasa = false;
options_pratt = false;
options_hough = false;
options_hough_matlab = false;

% Initial computations
n_all = numel(all_z);
info_point_spacing = norm(mean([all_x,all_y,all_z]))*t_sampling;

% Plot point subset
if options_verbose && options_point_subset;
    h = figure;
    scatter3(all_x, all_y, all_z, 10, 'b', 'filled');
    hold on
    axis auto
    daspect([1 1 1 ]);
    campos([0 0 0]);
    title('Initial Point Subset');
     pause
end
%}

if options_verbose && options_point_subset; delete(h); end
% Outputs:      options             - user-defined run-type options 
%               n_all               - number of points (all)
%               info_point_spacing  - spacing between laser scans 
%% Determine principle direction 

% Fit a 3d line to data using RANSAC
if numel(all_x)<4;
    return
end
%{
% Just return minimal output if no plotting
if ~(options_verbose && options_principal_direction);
    [line_z1,line_y1,~] = ransac_3line([all_x, all_y, all_z]',info_point_spacing);
end



% Plot the results
if options_verbose && options_principal_direction;
    [line_z1,line_y1,all_is_line_inlier] = ransac_3line([all_x, all_y, all_z]',info_point_spacing);
    h = figure;
    hold on;
    plot3([line_y1(1) line_z1(1)],[line_y1(2) line_z1(2)],[line_y1(3) line_z1(3)],'-k','linewidth',4)
    scatter3(all_x(all_is_line_inlier),all_y(all_is_line_inlier),all_z(all_is_line_inlier),10,'b','filled')
    scatter3(all_x(~all_is_line_inlier),all_y(~all_is_line_inlier),all_z(~all_is_line_inlier),10,'r','filled')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis auto
    daspect([1 1 1]);
    campos([0 0 0])
    title('RANSAC Line fit');
    pause
end
%}
%}


line_z1 = [mean(all_x),mean(all_y),mean(all_z)]';
line_y1 = line_z1 + [0;0;1];
% Plot the results
if options_verbose && options_principal_direction;
    h = figure;
    hold on;
    plot3([line_y1(1) line_z1(1)],[line_y1(2) line_z1(2)],[line_y1(3) line_z1(3)],'-k','linewidth',4)
    scatter3(all_x,all_y,all_z,10,'r','filled')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis auto
    daspect([1 1 1]);
    campos([0 0 0])
    title('RANSAC Line fit');
    pause
end

% Check if stem is vertical enough
line_a1 = line_z1 - line_y1;
line_a1norm = line_a1/norm(line_a1);
if line_a1norm(3)<0;
    line_a1norm = -line_a1norm;
end
line_theta = acosd(line_a1norm(3));
if line_theta >t_theta
    return
end

clear all_is_line_inlier line_a1 line_theta
if options_verbose && options_principal_direction; delete(h); end
% Outputs:  line_z1                 - lower (in z) point on a 3D line 
%           line_y1                 - upper (in z) point on a 3D line
%           line_a1norm             - normal vector in growth direction
%% Gram Schmidt Orthogonalization 1

% Initial Gram Schmidt Orthogonalization
all_xyz = [all_x all_y all_z];
% First vector is principal growth direction
q1 = line_a1norm;
% Second vector is from sensor to object
p2 = (line_z1+line_y1)/2;
w2 = p2-dot(p2,q1)*q1;
q2 = w2/norm(w2);
% Third vector orthogonal 
p3 = [1 1 1]';
w3 = p3 - dot(p3,q1)*q1 - dot(p3,q2)*q2;
q3 = w3/norm(w3);
all_w = dot(all_xyz,repmat(q1',[n_all,1]),2); % Projection onto a1norm
all_u = dot(all_xyz,repmat(q2',[n_all,1]),2); % Projection onto second basis vector
all_v = dot(all_xyz,repmat(q3',[n_all,1]),2); % Projection onto third basis vector
M1 = [q1'; q2'; q3'];

% Show initial orthogonalization 
if options_verbose && options_orthogonalization_1
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1);
    scatter(all_v,all_w,20,'b','filled');
    hold on
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    scatter(all_u,all_w,20,'b','filled');
    hold on;
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(all_u,all_v,20,'b','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Orthogonalization based on RANSAC');
    pause
end

clear q1 p2 w2 q2 p3 w3 q3 all_xyz line_z1 line_y1 line_a1norm
if options_verbose && options_orthogonalization_1; delete(h); end
% Outputs:  all_w                   - Projection of all points onto a1norm
%           all_u                   - Projection of points onto sensor line
%           all_v                   - Projection of points onto 3rd basis 
%           M1                      - Transformation from inital G.S.O.
%% Density Estimation 
% (could just use alpha shapes)

corners_wuv = M1*corners_xyz;
info_density_step = t_min_obj_size;%6*info_point_spacing;
ctrs_u = min(corners_wuv(2,:)):info_density_step:max(corners_wuv(2,:));
ctrs_v = min(corners_wuv(3,:)):info_density_step:max(corners_wuv(3,:));
ctrs_w = min(corners_wuv(1,:)):info_density_step:max(corners_wuv(1,:));

%{
figure;
scatter3(all_u, all_v, all_w,10,'r','filled');
hold on
scatter3(corners_wuv(2,:),corners_wuv(3,:),corners_wuv(1,:),50,'b','filled')
%}

% Compute histogram and dilate 
[ X2d_exp] = local_point_density(M1, ctrs_u, ctrs_v, ctrs_w );
if isempty(X2d_exp);
    return
end
X2d_exp = rot90(X2d_exp,-1);

% Find maximum azimuth and elevation differences 

ctrs_uv{1} = ctrs_u;
ctrs_uv{2} = ctrs_v;
[X2d_count,~,bins] = hist3_kelbe([all_u all_v],ctrs_uv);
%X2d_count = rot90(X2d_count);

%{
% Compute histogram and dilate 
[ X2d_exp, spacing_u, spacing_v ] = local_point_density( all_u,all_v,all_w, M1, info_point_spacing );
if isempty(X2d_exp);
    return
end
X2d_exp = rot90(X2d_exp,-1);

% Find maximum azimuth and elevation differences 

ctrs{1} = spacing_u;
ctrs{2} = spacing_v;
[X2d_count,~,bins] = hist3_kelbe([all_u all_v],ctrs);
%X2d_count = rot90(X2d_count);
%}

%n_points_exp_max = 9*t_r/info_point_spacing;
X2d_ratio = X2d_count./X2d_exp;
%is_dense= X2d_ratio>t_density_cutoff;
X2d_iscount = X2d_ratio > .2;
X2d_dilate = imdilate(X2d_iscount,strel('square',3));
bin_ind = sub2ind([numel(ctrs_u), numel(ctrs_v)],  bins(:,1),bins(:,2));
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

%if sum(is_dense)<50
%    return
%end

if options_verbose && options_density
    h = figure('position', [872   568   794   408]);
    subplot(1,2,1);
    scatter3(all_u(is_dense),all_v(is_dense),all_w(is_dense),20,all_ratio(is_dense),'filled');
    hold on
    scatter3(all_u(~is_dense),all_v(~is_dense),all_w(~is_dense),20,[.5 .5 .5],'filled');
    %scatter3(all_u(is_dense),all_v(is_dense),all_w(is_dense),20,'b','filled');
    %hold on
    %scatter3(all_u(~is_dense),all_v(~is_dense),all_w(~is_dense),20,'r','filled');
    daspect([1 1 1]);
    grid on; xlabel('v'); ylabel('w');
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


clear ctrs X2d_count bins bin_ind X2d_dilate X2d_iscount
if options_verbose && options_density; delete(h); end
% Outputs:  dense_w                 - Points of sufficient density 
%           dense_u                 - Points of sufficient density 
%           dense_v                 - Points of sufficient density 
%% Remove points if insufficient number
if numel(dense_u)<.9*t_min_pts;
    if options_remove_insufficient_pts;
    h = figure('position', [872   568   794   408]);
    scatter3(all_u(is_dense),all_v(is_dense),all_w(is_dense),20,all_ratio(is_dense),'filled');
    hold on
    scatter3(all_u(~is_dense),all_v(~is_dense),all_w(~is_dense),20,[.5 .5 .5],'filled');
    %scatter3(all_u(is_dense),all_v(is_dense),all_w(is_dense),20,'b','filled');
    %hold on
    %scatter3(all_u(~is_dense),all_v(~is_dense),all_w(~is_dense),20,'r','filled');
    daspect([1 1 1]);
    grid on; xlabel('v'); ylabel('w');
    colorbar; caxis([0 1])
    campos([0 0 0]);
    title('Removed: Insufficient number of points');
    pause
    end
    if options_verbose && options_remove_insufficient_pts; delete(h); end
    return
end
clear is_dense
%% Best alpha shape for uv projection
% Alpha shape of dense points in uv projection ('nadir' view)
%if isempty(dense_u)
%    foo = 1;
%    return
%end
[~,Suv_dense] = alphavol([dense_u dense_v],info_point_spacing);

%{
if isempty(Suv_dense.bnd) || isempty(Suv_dense.bnd);
    % [~,Svw] = alphavol([va wa],point_spacing,1);
    % [~,Suw] = alphavol([ua wa],point_spacing,1);
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

% Alpha shapes with "best" points for vw and uw projections
[~,Svw_alpha] = alphavol([alpha_v alpha_w],info_point_spacing);
[~,Suw_alpha] = alphavol([alpha_u alpha_w],info_point_spacing);

if isempty(Svw_alpha.bnd) || isempty(Suw_alpha.bnd);
    return
end
% Find largest alpha shape for remaining projections 
[ bnd_vw_x, bnd_vw_y, ~ ] = alpha_maxcomp( Svw_alpha.bnd, alpha_v,alpha_w );
[ bnd_uw_x, bnd_uw_y, ~ ] = alpha_maxcomp( Suw_alpha.bnd, alpha_u,alpha_w );

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
    % plot(Svw.x(Svw.b)',Svw.y(Svw.b)','-r','LineWidth',2)
    plot(bnd_vw_x,bnd_vw_y,'-r','LineWidth',2)
    hold on
    scatter(dense_v(is_alpha),dense_w(is_alpha),10,'b','filled');
    scatter(dense_v(~is_alpha),dense_w(~is_alpha),10,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    plot(bnd_uw_x,bnd_uw_y,'-r','LineWidth',2)
    hold on
    scatter(dense_u(is_alpha),dense_w(is_alpha),10,'b','filled');
    scatter(dense_u(~is_alpha),dense_w(~is_alpha),10,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    plot(bnd_uv_x_alpha,bnd_uv_y_alpha,'-r','LineWidth',2)
    hold on
    scatter(dense_u(is_alpha),dense_v(is_alpha),10,'b','filled');
    scatter(dense_u(~is_alpha),dense_v(~is_alpha),10,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Alpha shapes');
    %pause
end

clear Suv_alpha Svw_alpha Suw_alpha bnd_uv_x_alpha bnd_uv_y_alpha
if options_verbose && options_alpha_shapes; delete(h); end
% Outputs:  bnd_vw_x                - x Boundary points for vw projection
%           bnd_vw_y                - y Boundary points for vw projection
%           bnd_uw_x                - x Boundary points for uw projection
%           bnd_uw_y                - y Boundary points for uw projection
%% Find side points 

f_base = 45; % Fixed value 360/8

% Array of adjacent edge points 
edge_vw_x = [bnd_vw_x(1:end-1) circshift(bnd_vw_x(1:end-1),[-1 0])];
edge_vw_y = [bnd_vw_y(1:end-1) circshift(bnd_vw_y(1:end-1),[-1 0])];
edge_uw_x = [bnd_uw_x(1:end-1) circshift(bnd_uw_x(1:end-1),[-1 0])];
edge_uw_y = [bnd_uw_y(1:end-1) circshift(bnd_uw_y(1:end-1),[-1 0])];
% Difference of adjacent edge points 
delta_vw_x = edge_vw_x(:,1) - edge_vw_x(:,2);
delta_vw_y = edge_vw_y(:,1) - edge_vw_y(:,2);
delta_vw = [delta_vw_x delta_vw_y];
delta_vw = delta_vw./repmat(sqrt(delta_vw(:,1).^2 + delta_vw(:,2).^2),[1,2]);
delta_uw_x = edge_uw_x(:,1) - edge_uw_x(:,2);
delta_uw_y = edge_uw_y(:,1) - edge_uw_y(:,2);
delta_uw = [delta_uw_x delta_uw_y];
delta_uw = delta_uw./repmat(sqrt(delta_uw(:,1).^2 + delta_uw(:,2).^2),[1,2]);
% Angle of edge segments  
theta_vw =  atan2d(delta_vw(:,1),delta_vw(:,2));
theta_uw = atan2d(delta_uw(:,1),delta_uw(:,2));
% Logical defining side for edge points 
is_sidevw1 = (theta_vw<f_base&theta_vw>-f_base);
is_sidevw2 = (theta_vw>180-f_base)|(theta_vw<-180+f_base);
is_sideuw1 = (theta_uw<f_base&theta_uw>-f_base);

% Array of side points 
side1_vw_x = edge_vw_x(is_sidevw1);
side1_vw_y = edge_vw_y(is_sidevw1);
side2_vw_x = edge_vw_x(is_sidevw2);
side2_vw_y = edge_vw_y(is_sidevw2);
side1_uw_x = edge_uw_x(is_sideuw1);
side1_uw_y = edge_uw_y(is_sideuw1);

if options_verbose && options_side_points
h = figure('position',[865 549 806 405]);
subplot(1,2,1);
hold on
scatter(alpha_v, alpha_w,50,[.5 .5 .5])
scatter(edge_vw_x(is_sidevw1,1), edge_vw_y(is_sidevw1,1),50,'r','filled')
is_noside = (~is_sidevw1 & ~is_sidevw2);
scatter(edge_vw_x(is_noside,1), edge_vw_y(is_noside,1),50,[.5 .5 .5],'filled')
scatter(edge_vw_x(is_sidevw2,1), edge_vw_y(is_sidevw2,1),50,'g','filled')
daspect([1 1 1]);
xlabel('v'); ylabel('w');
subplot(1,2,2);
hold on
scatter(alpha_u, alpha_w,50,[.5 .5 .5])
scatter(edge_uw_x(is_sideuw1,1), edge_uw_y(is_sideuw1,1),50,'r','filled')
scatter(edge_uw_x(~is_sideuw1,1), edge_uw_y(~is_sideuw1,1),50,[.5 .5 .5],'filled')
xlabel('u'); ylabel('w');
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
    %  fprintf('\n Returned because %3.0f Points for ransac \n', min(sum(is_sidevw1),sum(is_sidevw2)));
    return;
end
% RANSAC
if options_verbose && options_RANSAC_axis_fit
    [S1_a_vw,S1_b_vw, is_inlier_vw1,E_vw1] = ransac_2line([side1_vw_x side1_vw_y]',info_point_spacing/2);
    [S2_a_vw,S2_b_vw, is_inlier_vw2,E_vw2] = ransac_2line([side2_vw_x side2_vw_y]',info_point_spacing/2);
    [S1_a_uw,S1_b_uw, is_inlier_uw1,E_uw1] = ransac_2line([side1_uw_x side1_uw_y]',info_point_spacing/2);
else
    [S1_a_vw,S1_b_vw, ~,~] = ransac_2line([side1_vw_x side1_vw_y]',info_point_spacing/2);
    [S2_a_vw,S2_b_vw, ~,~] = ransac_2line([side2_vw_x side2_vw_y]',info_point_spacing/2);
    [S1_a_uw,S1_b_uw, ~,~] = ransac_2line([side1_uw_x side1_uw_y]',info_point_spacing/2);
end

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
    plot(s12x_vw, s12y_vw,'-m','linewidth',2)
    plot(s1x_vw, s1y_vw,'-m','linewidth',2)
    plot(s2x_vw, s2y_vw,'-m','linewidth',2)
    scatter(alpha_v,alpha_w,50,[.5 .5 .5]);
    %scatter(side1_vw_x(is_inlier_vw1),side1_vw_y(is_inlier_vw1),50,E_vw1(is_inlier_vw1),'filled');
    %scatter(side1_vw_x(~is_inlier_vw1),side1_vw_y(~is_inlier_vw1),50,E_vw1(~is_inlier_vw1));
    scatter(side1_vw_x(is_inlier_vw1),side1_vw_y(is_inlier_vw1),50,'b','filled');
    scatter(side1_vw_x(~is_inlier_vw1),side1_vw_y(~is_inlier_vw1),50,'r','filled');
    scatter(side2_vw_x(is_inlier_vw2),side2_vw_y(is_inlier_vw2),50,'b','filled');
    scatter(side2_vw_x(~is_inlier_vw2),side2_vw_y(~is_inlier_vw2),50,'r','filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,2,2)
    hold on
    plot(s1x_uw, s1y_uw,'-m','linewidth',2)
    scatter(alpha_u,alpha_w,20,[.5 .5 .5]);
    scatter(side1_uw_x(is_inlier_uw1),side1_uw_y(is_inlier_uw1),50,'b','filled');
    scatter(side1_uw_x(~is_inlier_uw1),side1_uw_y(~is_inlier_uw1),50,'r','filled');
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

%% Adjust side vectors of face-on (vw) view using Homography

% Find the homography H where b = H*x
x_vw = [all_v all_w ones(n_all,1)]'; %3xn
b_vw = [orthog_v orthog_w ones(numel(all_v),1)]'; %3xn
H_vw = b_vw/x_vw;
x_uw = [all_u all_w ones(n_all,1)]'; %3xn
b_uw = [orthog_u orthog_w ones(numel(all_u),1)]'; %3xn
H_uw = b_uw/x_uw;

% Concatenate side vector 
s_vw = [ s1x_vw s2x_vw s12x_vw ; s1y_vw s2y_vw s12y_vw; 1 1 1 1 1 1];
s_uw = [ s1x_uw ;s1y_uw; 1 1 ];
% New side vector
s_vw_new = H_vw*s_vw;
s_uw_new = H_uw*s_uw;
% Separate into x,y components 
s1x_vw = s_vw_new(1,1:2);
s2x_vw = s_vw_new(1,3:4);
s12x_vw = s_vw_new(1,5:6);
s1y_vw = s_vw_new(2,1:2);
s2y_vw = s_vw_new(2,3:4);
s12y_vw = s_vw_new(2,5:6);
s1x_uw = s_uw_new(1,1:2);
s1y_uw = s_uw_new(2,1:2);

if options_verbose && options_homography
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1);
    hold on
    plot(s12x_vw, s12y_vw,'-m','linewidth',2)
    scatter(orthog_v,orthog_w,20,'b','filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    hold on
    scatter(orthog_u,orthog_w,20,'b','filled');
    plot(s1x_uw, s1y_uw,'-m','linewidth',2)
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(orthog_u,orthog_v,20,'b','filled');
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Homography');
    pause
end

clear b H x s s_new s12x_vw s12y_vw
if options_verbose && options_homography; delete(h); end
% Outputs:  s1x_vw                  - Corrected side 1 x values 
%           s2x_vw                  - Corrected side 2 x values 
%           s1y_vw                  - Corrected side 1 y values 
%           s2y_vw                  - Corrected side 2 y values 
%% Compute radii of alpha shapes 

% y values
y_step = [min([s1y_vw s2y_vw]) max([s1y_vw s2y_vw])];
% Slopes of left and right lines 
m1 = (s1x_vw(2) - s1x_vw(1))/(s1y_vw(2) - s1y_vw(1));
m2 = (s2x_vw(2) - s2x_vw(1))/(s2y_vw(2) - s2y_vw(1));
% Intercepts of left and right lines 
b1 = (-s1y_vw(1)*m1) + s1x_vw(1);
b2 = (-s2y_vw(1)*m2) + s2x_vw(1);
% x-values 
s1_step = m1*y_step + b1;
s2_step = m2*y_step + b2;

% Compute radius and center 
r_alpha = mean(abs([s1_step' s2_step']),2);
w_alpha = y_step';
n_step_alpha = numel(r_alpha);
u_alpha = s1x_uw(1) + r_alpha;
v_alpha = mean(s1x_vw)+r_alpha;

if options_verbose && options_circ
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1);
    hold on
    v_val = [-r_alpha r_alpha];
    plot(v_val',[w_alpha w_alpha]','-b','Linewidth',2);
    plot(v_val,[w_alpha w_alpha],'-b','Linewidth',2);
    scatter(orthog_v,orthog_w,50,[.5 .5 .5],'filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    hold on
    scatter(orthog_u,orthog_w,50,[.5 .5 .5],'filled');    
    u_val = [u_alpha u_alpha] + [-r_alpha r_alpha];
    plot(u_val',[w_alpha w_alpha]','-b','Linewidth',2);
    plot(u_val,[w_alpha w_alpha],'-b','Linewidth',2);
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(orthog_u,orthog_v,50,[.5 .5 .5],'filled');
    for ws = 1:numel(r_alpha);
        rectangle('position',[u_alpha(ws)-r_alpha(ws),v_alpha(ws)-r_alpha(ws),r_alpha(ws)*2,r_alpha(ws)*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','b','linewidth',2);
    end
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Alpha shapes circle fitting');
    pause
end

clear y_step s1y_vw s2y_vw s1x_vw s2x_vw m1 m2 b1 b2 s1_step s2_step n_step_alpha
if options_verbose && options_circ; delete(h); end
% Outputs:  r_alpha                 - Circle radii using alpha shapes  
%           w_alpha                 - Cirle heights using alpha shapes  
%           u_alpha                 - Circle center (u) using alpha shapes 
%           v_alpha                 - Circle cetner (v) using alpha shapes 
%% Compute radii using other methods 

% Kasa method
if options_kasa
    [u_kasa,v_kasa,r_kasa,~,~] = ransac_circle([orthog_u orthog_v]',t_kasa);
end
% Pratt method 
if options_pratt
    Par = CircleFitByPratt([orthog_u orthog_v]);
    u_pratt = Par(1);
    v_pratt = Par(2);
    r_pratt = Par(3);
end 
%{
if options_hough_matlab;
    radius_range = [40 60];
    isize = 100;
    image = zeros(isize,isize);
    un_demean = round(isize*(un-min(un))/max((un-min(un))));
    vn_demean = round(isize*(vn-min(vn))/max((vn-min(vn))));
    for p = 1:numel(un);
        if un_demean(p) ~=0 && vn_demean(p) ~=0;
            image(un_demean(p),vn_demean(p)) = image(un_demean(p),vn_demean(p))+1;
        end
    end
    [centers, radii] = imfindcircles(image,radius_range);
end
if options_hough;
    [u_Hough, v_Hough, r_Hough] = Tree_Detection_fxn(un',vn');
end
    
%}

if options_verbose && options_circzoom
    h = figure('position',[1005 56 557 472]);
    hold on
    scatter(orthog_u,orthog_v,20,'b','filled');
    for ws = 1:n_step_alpha;
        rectangle('position',[u_alpha(ws)-r_alpha(ws),v_alpha(ws)-r_alpha(ws),r_alpha(ws)*2,r_alpha(ws)*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','g','linewidth',2);
    end
    axis(axis);
    if options_kasa
        rectangle('position',[u_kasa-r_kasa,v_kasa-r_kasa,r_kasa*2,r_kasa*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','black','linewidth',2);
    end
    if options_pratt
        rectangle('position',[u_pratt-r_pratt,v_pratt-r_pratt,r_pratt*2,r_pratt*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','m','linewidth',2);
    end
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    title('Detail circles');
    pause
end

clear Par 
if options_verbose && options_circzoom; delete(h); end
% Outputs:  r_kasa                  - Circle radii using kasa method  
%           u_kasa                  - Circle center (u) using kasa method
%           v_kasa                  - Circle cetner (v) using kasa method
%           r_pratt                 - Circle radii using pratt method  
%           u_pratt                 - Circle center (u) using pratt method
%           v_pratt                 - Circle cetner (v) using pratt method
%% Find inliners to model_cyl function  

u_m = mean(u_alpha);
v_m = mean(v_alpha);
r_m = mean(r_alpha);

X2 = [orthog_u orthog_v]';
X2 = X2 - repmat([u_m; v_m],[1, n_all]);
E_alpha =  abs(sqrt(sum(X2.^2,1))-r_m)';
is_inlier_alpha = E_alpha<2*info_point_spacing;

if options_verbose && options_inliers
    h = figure('position', [872   568   794   408]);
    subplot(1,3,1);
    hold on
    v_val = [-r_alpha r_alpha];
    plot(v_val',[w_alpha w_alpha]','-m','Linewidth',2);
    plot(v_val,[w_alpha w_alpha],'-m','Linewidth',2);
    scatter(orthog_v(is_inlier_alpha),orthog_w(is_inlier_alpha),50,'b','filled');
    scatter(orthog_v(~is_inlier_alpha),orthog_w(~is_inlier_alpha),50,'r','filled');
    daspect([1 1 1]);
    xlabel('v'); ylabel('w');
    subplot(1,3,2)
    hold on
    scatter(orthog_u(is_inlier_alpha),orthog_w(is_inlier_alpha),50,'b','filled');
    scatter(orthog_u(~is_inlier_alpha),orthog_w(~is_inlier_alpha),50,'r','filled');    
    u_val = [u_alpha u_alpha] + [-r_alpha r_alpha];
    plot(u_val',[w_alpha w_alpha]','-m','Linewidth',2);
    plot(u_val,[w_alpha w_alpha],'-m','Linewidth',2);
    daspect([1 1 1]);
    xlabel('u'); ylabel('w');
    subplot(1,3,3)
    hold on
    scatter(orthog_u(is_inlier_alpha),orthog_v(is_inlier_alpha),50,'b','filled');
    scatter(orthog_u(~is_inlier_alpha),orthog_v(~is_inlier_alpha),50,'r','filled');    for ws = 1:numel(r_alpha);
        rectangle('position',[u_alpha(ws)-r_alpha(ws),v_alpha(ws)-r_alpha(ws),r_alpha(ws)*2,r_alpha(ws)*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','m','linewidth',2);
    end
    daspect([1 1 1]);
    xlabel('u'); ylabel('v');
    suptitle('Alpha shapes circle fitting');
    pause
end

clear u_m v_m r_m X2 
if options_verbose && options_inliers; delete(h); end
% Outputs:  E_alpha                 - Error for all points from circle
%           is_inlier_alpha         - Inliers based on error threshold 
%% Transform back to original coordinates and find circle inliers 

M = M3*M2*M1;
xyz_alpha = M\[ w_alpha u_alpha v_alpha]';
xyz_left  = M\[ w_alpha u_alpha v_alpha + r_alpha]';
xyz_right  = M\[ w_alpha u_alpha v_alpha - r_alpha]';


% Outputs:  xyz_alpha               - Circle coordinates from Alpha-shapes
%           M                       - Transformation back to xyz
end

