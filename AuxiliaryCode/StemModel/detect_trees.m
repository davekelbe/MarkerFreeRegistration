function  [seg_xb, seg_yb,seg_zb, seg_z, seg_y, seg_r, seg_iter, seg_index,...
    seg_fill, seg_taper, seg_min_obj_size,...
    seg_z_left, seg_y_left, seg_z_right, seg_y_right, seg_status] = detect_trees(...
    data_x,data_y,data_z,data_z0,data_dem,...
    I12ieq,axis_a,axis_e, ...
    t_r, t_zmin0, t_zmax0, path_top,...
    t_theta,t_r1_max,t_taper_min, t_taper_max,t_coverage,...
    t_min_density)
%DETECT_TREES Wrapper function for model_cyl.m
%
% ************* Inputs **************
%
% data_x                - x values from lidar
% data_y                - y values from lidar
% data_z                - z values from lidar
% data_z0               - z0 values (height above ground)
% data_dem              - height of DEM at each xyz point
% I12ieq                - Intensity equalized Andrieu Projection Image
% axis_a                - Vector of azimuth values for Andrieu Projection
% axis_e                - Vector of elevation values for Andrieu Projection
% t_r                   - Voxel size [m]
% t_zmin0               - Minimum height above ground
% t_zmax0               - Maximum height above ground
% path_top              - Top-level directory for saving output
% t_theta               - maximum user-defined lean angle [degrees]
% t_r1_max              - maximum user_defined radius [m]
% t_r1_min              - minimum user-defined radius [m]
% t_taper_min           - minimum user-defined taper [degrees]
% t_taper_max           - maximum user-defined taper [degrees]
% t_coverage            - Minimum user-defined coverage [%]
% t_min_density         - minimum density is sub-voxel [%]
%
% ************* Outputs *************
% seg_xb                - x block number
% seg_yb                - y block number
% seg_zb                - z block number
% seg_z                 - Upper axis point of cylinder
% seg_y                 - Lower axis point of cylinder
% seg_r                 - Radius of tree cylinder
% seg_iter              - Iteration (1d block number) of fit
% seg_index             - Inlier points of cylinder fitting
% seg_fill              - Fill percentage/100
% seg_taper             - Taper of cylinder
% seg_min_obj_size      - Minimum object size at fitted range
% seg_z_left            - Upper left-axis point of cylinder
% seg_z_right           - Upper right-axis point of cylinder
% seg_y_left            - Lower left-axis point of cylinder
% seg_y_right           - Lower right-axis point of cylinder
%
% ************* Example ***************
% ************* Copyright *************
% Copyright (C) 2013, David Kelbe, Rochester Institute of Technology


%% Initialization
%
% Set options
options_print = false;
options_verbose = false;
options_save_eps = false;
options_save_fig = false;
options_save_png = false;
options_point_subset = false;
options_andrieu1 = false;
options_radius = false;
options_vertical = false;
options_extend = false;
options_outside = false;
options_bins = false;
options_vertical2 = false;
options_final = false;
options_andrieu = false;
%}

% Set options
%{
options_print = true;
options_verbose = true;
options_save_eps = false;
options_save_fig = false;
options_save_png = false;
options_point_subset = true;
options_andrieu1 = false;
options_radius = true;
options_vertical = true;
options_extend = true;
options_outside = true;
%options_bins = true;
options_vertical2 = true;
options_final = true;
options_andrieu = false;
%}

% Additional info
info_slash = '\';
info_tick = t_r/4; % Set three tick marks on plot
[info_sampling,~,~,~] = SICK_specs();

% Set paths
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
path_fig = sprintf('%s%s%s',path_top,'fig',info_slash);
path_eps = sprintf('%s%s%s',path_top,'eps',info_slash);
path_png = sprintf('%s%s%s',path_top,'png',info_slash);
path_eps_detect_trees = path_eps; % sprintf('%s%sdetect_trees%s',info_slash, path_eps,info_slash);
path_png_detect_trees = path_png; % sprintf('%s%sdetect_trees%s',info_slash, path_png,info_slash);

%Make non existant folders
if ~exist(path_mat,'dir');
    command = sprintf('mkdir %s',path_mat);
    [~,~] = system(command);
end
if ~exist(path_fig,'dir');
    command = sprintf('mkdir %s',path_fig);
    [~,~] = system(command);
end
if ~exist(path_eps,'dir');
    command = sprintf('mkdir %s',path_eps);
    [~,~] = system(command);
end
if ~exist(path_png,'dir');
    command = sprintf('mkdir %s',path_png);
    [~,~] = system(command);
end
if ~exist(path_eps_detect_trees,'dir');
    command = sprintf('mkdir %s',path_eps_detect_trees);
    [~,~] = system(command);
end
if ~exist(path_png_detect_trees,'dir');
    command = sprintf('mkdir %s',path_png_detect_trees);
    [~,~] = system(command);
end

clear command
% Outputs:      options             - user-defined run-type options
%               paths               - number of points (all)
%               info                - system specific info
%% Data preparation

% Initialize intensity image
I12ieq = double(I12ieq);
I12ieq = I12ieq./max(I12ieq(:));

% Subdivide data into blocks
block_x_lo = floor(min(data_x)):t_r:max(data_x);
block_x_hi = circshift(block_x_lo,[0,-1]);
block_x_hi(end) = block_x_lo(end) + t_r;
block_y_lo = floor(min(data_y)):t_r:max(data_y);
block_y_hi = circshift(block_y_lo,[0,-1]);
block_y_hi(end) = block_y_lo(end) + t_r;
block_z_lo = floor(min(data_z)):t_r:max(data_z);
block_z_hi = circshift(block_z_lo,[0,-1]);
block_z_hi(end) = block_z_lo(end) + t_r;
n_x_block = numel(block_x_lo);
n_y_block = numel(block_y_lo);
n_z_block = numel(block_z_lo);

% Initialize output variables
info_maxiter = 2*n_x_block*n_y_block*n_z_block; %Pad to 2x size to allow multiple stems
seg_xb = zeros(info_maxiter,1); % x block number
seg_yb = zeros(info_maxiter,1); % y block number
seg_zb = zeros(info_maxiter,1); % z block number
seg_z = zeros(info_maxiter,3); % 3-vector with lower center axis point
seg_y = zeros(info_maxiter,3); % 3-vector with upper center axis point
seg_r = zeros(info_maxiter,2); % radius
seg_iter = zeros(info_maxiter,1); % main iteration
seg_index = cell(info_maxiter,1); % index to inlier points
seg_fill = zeros(info_maxiter,1); % Fill percentage based on sampling
seg_taper = zeros(info_maxiter,1); % Taper angle
seg_min_obj_size = zeros(info_maxiter,1); % Minimum object size
seg_z_left = zeros(info_maxiter,3); % 3-vector with lower left axis point
seg_y_left = zeros(info_maxiter,3); % 3-vector with upper left axis point
seg_z_right = zeros(info_maxiter,3); % 3-vector with lower right axis point
seg_y_right = zeros(info_maxiter,3); % 3-vector with upper right axis point
seg_status = cell(info_maxiter,1); % index to inlier points

% Other preparation
n_pts = numel(data_x);
data_is_removed = false(n_pts,1);
data_index = (1:n_pts)'; %Index to data (within function)
s = 1; % Detected candidate tree number

clear t_maxiter
% Outputs:      block_x_lo          - Array of lower bounds for voxelization
%               block_x_hi          - Array of upper bounds for voxelization
%               block_y_lo          - Array of lower bounds for voxelization
%               block_y_hi          - Array of upper bounds for voxelization
%               block_z_lo          - Array of lower bounds for voxelization
%               block_z_hi          - Array of upper bounds for voxelization
%               n_x_block           - Number of blocks in x dimension
%               n_y_block           - Number of blocks in y dimension
%               n_z_block           - Number of blocks in z dimension
%               seg_xb              - Empty array to hold x block number
%               seg_yb              - Empty array to hold y block number
%               seg_zb              - Empty array to hold z block_number
%               seg_z               - Empty array to hold lower axis point
%               seg_y               - Empty array to hold upper axis point
%               seg_r               - Empty array to hold radius
%               seg_iter            - Empty array to hold iteration number
%               seg_index           - Empty array to hold index
%               seg_fill            - Empty array to hold fill percentage
%               n_pts               - Number of points
%               data_infxn_is_inlier- Logical array of inlier
%               data_infxn_index    - Index to points in original data
%               s                   - segment number (iterated)

%% Iterate through each block
%{
xyzb = [12 7 3;1 9 6; 4 7 4; 4 10 2; 5 11 2; 6 7 4; 7 4 2; 7 4 3; ...
    9 16 2; 10 3 5; 12 4 2; 14 15 2; 6 5 3; 9 4 4; 9 16 2; ...
    7 3 3; 7 4 3; 6 4 7; 8 13 5; 9 6 7; 10 8 4; 15 13 2; ...
    1 9 7; 2 11 3; 2 11 4; 2 11 6; 3 9 5; 3 9 6; 4 5 2; 4 5 7;...
    4 7 3; 6 11 2; 6 15 3; 7 3 6; 8 11 13;...
    9 12 13; 9 12 4; 10 7 7; 12 12 9; 14 2 1; 14 2 3; 14 3 6; 14 4 3; 15 9 7];
for i = 1:size(xyzb,1)
        xb = xyzb(i,1);
        yb = xyzb(i,2);
        zb = xyzb(i,3);
        %xb = 5; yb = 11; zb = 2;
%}
for xb = 1:n_x_block;
    for yb = 1:n_y_block; %9
        for zb =1:n_z_block; %6
            %xb = 4; yb =5; zb = 2;
            
            % Center boundaries
            x_loc = block_x_lo(xb);
            x_hic = block_x_hi(xb);
            y_loc = block_y_lo(yb);
            y_hic = block_y_hi(yb);
            z_loc = block_z_lo(zb);
            z_hic = block_z_hi(zb);
            % Outer (block) boundaries
            x_lob = max(block_x_lo(xb)-1,block_x_lo(1));
            x_hib = min(block_x_hi(xb)+1,block_x_hi(end));
            y_lob = max(block_y_lo(yb)-1,block_y_lo(1));
            y_hib = min(block_y_hi(yb)+1,block_y_hi(end));
            z_lob = max(block_z_lo(zb)-1,block_z_lo(1));
            z_hib = min(block_z_hi(zb)+1,block_z_hi(end));
            
            %% Continue in block
            
            % Continue until all potential trees are exhausted
            badfit = false;
            bn = 1; %Block number (iterated)
            block_iter = 1;
            while badfit==false;
                %% Reset general info
                % Reset EPS and PNG number for each block
                eps_n = 1; %EPS Image number
                png_n = 1; %PNG Image number
                cylcolor = 'g'; % Assume "fit"
                
                % badfit = true; %Uncomment to run each box only once
                %%  Subset data based on block
                data_is_sub = data_x>=x_lob & data_x<x_hib&...
                    data_y>=y_lob & data_y<y_hib &...
                    data_z>=z_lob & data_z<z_hib &...
                    data_z0>=t_zmin0&data_z0<t_zmax0&~data_is_removed;
                n_sub = sum(data_is_sub);
                
                % Adjust z_loc based on DEM
                if block_iter == 1;
                    t_r_adj = t_r;
                    mean_dem = mean(data_dem(data_is_sub));
                    if mean_dem+t_zmin0>z_loc;
                        z_loc = mean_dem+t_zmin0;
                        z_lob = mean_dem+t_zmin0;
                        t_r_adj = z_hic-z_loc;
                        t_taper_max = t_taper_max + 4;
                        if (t_r_adj <= 0);
                            badfit = true;
                            continue;
                        end
                    end
                end
                
                %mean_range = norm([(x_loc + x_hic)/2, (y_loc + y_hic)/2, (z_loc + z_hic)/2]);
                mean_range = mean(sqrt(sum([x_loc x_loc x_loc x_loc x_hic x_hic x_hic x_hic;...
                    y_loc y_loc y_hic y_hic y_loc y_loc y_hic y_hic;...
                    x_loc z_hic z_loc z_hic z_loc z_hic z_loc z_hic].^2,1)));
                [ ~, ~, t_min_obj_size, t_min_pts ] = SICK_specs( mean_range, t_r_adj);
                
                % Check if there are enough points
                if ~options_point_subset % If images need to be plotted, wait to continue until after plot
                    if n_sub < t_coverage*t_min_pts;
                        if options_print
                            fprintf('\nBlock (%2.0f,%2.0f,%2.0f) aborted due to insufficient # points\n',xb,yb,zb);
                        end
                        badfit = true;
                        continue
                    end
                end
                if options_print
                    fprintf('\nBlock (%2.0f,%2.0f,%2.0f) from [%2.0f,%2.0f);[%2.0f,%2.0f);[%2.0f,%2.0f)\n',xb,yb,zb,x_loc,x_hic,y_loc,y_hic,z_loc,z_hic);
                end
                % Subset data
                sub_x = data_x(data_is_sub);
                sub_y = data_y(data_is_sub);
                sub_z = data_z(data_is_sub);
                sub_index = data_index(data_is_sub); % Index to points in data (in function)
                
                % Find center of subset
                sub_iscenter = sub_x>=x_loc & sub_x<x_hic&...
                    sub_y>=y_loc & sub_y<y_hic&...
                    sub_z>=z_loc & sub_z<z_hic;
                n_center = sum(sub_iscenter);
                
                % Determine axes
                sub_axes = [x_loc x_hic y_loc y_hic z_loc z_hic];
                % Corners
                %{
                corners = [  x_loc y_loc z_loc;...
                        x_hic y_loc z_loc;...
                        x_hic y_hic z_loc;...
                        x_loc y_hic z_loc;...
                        x_loc y_loc z_hic;...
                        x_hic y_loc z_hic;...
                        x_hic y_hic z_hic;...
                        x_loc y_hic z_hic]';
                %}
                
                % Plot point subset
                if options_point_subset && n_center>10 ;%&& n_center<t_minpts;
                    if options_verbose
                        if n_center<t_coverage*t_min_pts;
                            h1 = figure('color','r','position',[910   650   682   284]);
                        else
                            h1 = figure('color','g','position',[910   650   682   284]);
                        end
                    elseif any([options_save_eps, options_save_png, options_save_fig]);
                        h1 = figure('visible','off','color','w');
                    end
                    scatter3(sub_x(sub_iscenter), sub_y(sub_iscenter), sub_z(sub_iscenter), 10, 'b', 'filled');
                    hold on
                    scatter3(sub_x(~sub_iscenter), sub_y(~sub_iscenter), sub_z(~sub_iscenter), 10, 'r', 'filled');
                    campos([0 0 0])
                    if n_center < t_coverage*t_min_pts;
                        pf = 'FAIL';
                    else
                        pf = 'PASS';
                    end
                    title(sprintf('(%g, %g, %g); Subset: %s: %g points',xb,yb,zb, pf, n_center));
                    xlabel('x'); ylabel('y'); zlabel('z');
                    legend('center','exterior','location','east')
                    daspect([1 1 1]);
                    
                    [ I12ieq_c ] = color_andrieu( I12ieq,axis_a, axis_e, sub_x, sub_y, sub_z, sub_iscenter );
                    h2 = figure('color','w');
                    imshow(I12ieq_c)
                    colormap('gray');
                    axis off; axis image
                    set(gca,'Units','normalized','Position',[0 0 1 1]);
                    [n_row, n_col,~] = size(I12ieq);
                    set(gcf,'Units','pixels','Position',[900 305 n_col/2 n_row/2]);
                    if options_save_eps;
                        if n_center < t_coverage*t_min_pts;
                            filename_eps_point_subset_reject = sprintf('%sx%02.0f-y%02.0f-z%02.0f_b%01.0f_i%01.0f_pointsubset-reject.eps',...
                                path_eps_detect_trees, xb, yb, zb, bn, eps_n);
                            print(h1,'-depsc','-painters',filename_eps_point_subset_reject)
                        else
                            filename_eps_point_subset_accept = sprintf('%sx%02.0f-y%02.0f-z%02.02f_b%01.0f_i%01.0f_pointsubset-accept.eps',...
                                path_eps_detect_trees, xb, yb, zb, bn, eps_n);
                            print(h1,'-depsc','-painters',filename_eps_point_subset_accept)
                        end
                        eps_n = eps_n + 1;
                    end
                    pause
                    if  any([options_save_eps, options_save_png, options_save_fig, options_verbose])
                        delete(h1); delete(h2);
                    end
                end
                
                clear mean_range
                %% Andrieu Image of subset
                
                % Color subset on Andrieu image
                %{
                if options_andrieu1 && n_center>3;
                    axis_c = 1:numel(axis_a);
                    axis_r = 1:numel(axis_e);
                    [sub_a, sub_e,~] = cart2sph(sub_x(sub_iscenter),sub_y(sub_iscenter),sub_z(sub_iscenter));
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
                    if options_verbose
                        h = figure('color','w');
                        imshow(I12ieq_c)
                        colormap('gray');
                        axis off; axis image
                        set(gca,'Units','normalized','Position',[0 0 1 1]);
                        [n_row, n_col,~] = size(I12ieq);
                        set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
                    end
                    if options_save_png;
                        filename_png_andrieu1 = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_andrieu1.png',...
                            path_png_detect_trees, xb, yb, bn, png_n);
                        imwrite(I12ieq_c,filename_png_andrieu1,'png');
                        png_n = png_n + 1;
                    end
                    if  any([options_save_fig, options_verbose])
                        delete(h);
                    end
                end
                %}
                
                if sum(sub_iscenter) < t_coverage*t_min_pts;
                    if options_print;
                        fprintf('\nBlock (%g,%g,%g) aborted due to insufficient # points\n',xb,yb,zb);
                    end
                    badfit = true;
                    continue
                end
                
                %% Initial RANSAC Cylinder fit to center points
                
                t_r1_min = .4*(t_min_obj_size/2);
                options_verbose_model_cyl = false;
                [cyl_ctr,rcirc,~,~,M,~,~,~,~,return_str,~,~,~] = model_cyl_allcomp(sub_x(sub_iscenter),sub_y(sub_iscenter),sub_z(sub_iscenter),...
                    t_r_adj,t_theta,t_r1_max,t_r1_min,t_taper_min,...
                    t_taper_max,info_sampling,t_min_density,t_min_obj_size,...
                    t_min_pts, t_coverage, sub_axes, options_verbose_model_cyl,[],...
                    [],[],[]);
                % Convert to coarse estimate for consistency with code
                if ~(return_str.empty_histogram || return_str.npts_dense ||...
                        return_str.empty_alpha_shape || return_str.side_points);
                    r1 = rcirc; % [bottom radius, top radius]
                    z1 = cyl_ctr(:,end); % top
                    y1 = cyl_ctr(:,1); % bottom
                    a1 = z1 - y1;
                    a1norm = a1/norm(a1);
                    t1 = atand((r1(1)-r1(2))/norm(z1-y1)); % positive values are normal
                else
                    badfit = true;
                    continue
                end
                
                %{
                n_seg = numel(rcirc)-1;
                figure;
                scatter3(sub_x(sub_iscenter),sub_y(sub_iscenter),sub_z(sub_iscenter),10,'r','filled');
                caxis manual
                hold on
                daspect([1 1 1]);
                for s = 1:n_seg
                    [Coneh,End1h,End2h] = Cone(cyl_ctr(:,s+1),cyl_ctr(:,s),[rcirc(s+1) rcirc(s)],30,cylcolor,1,0);
                    set(Coneh,'facealpha',.2);
                    set(End1h,'facealpha',.2);
                    set(End2h,'facealpha',.2);
                end
                %}
                
                %% Recompute error for all points
                X = [sub_x sub_y sub_z]';
                Xnew = M*X;
                znew = M*z1;
                X2 = [Xnew(2,:); Xnew(3,:)];
                X2 = X2 - repmat([znew(2); znew(3)],[1, n_sub]);
                E =  abs(sqrt(sum(X2.^2,1))-mean(r1))'; % Need to fix
                
                % Check if block needs to be extended
                %{
                % Check if block needs to be extended
                t_ext_max = 0.3; %Changed v4 from 0.3 to 0.2
                r_ext = max(4*r1,t_ext_max);
                x_lo_ext = min(z1(1),y1(1))-r_ext;
                x_hi_ext = max(z1(1),y1(1))+r_ext;
                y_lo_ext = min(z1(2),y1(2))-r_ext;
                y_hi_ext = max(z1(2),y1(2))+r_ext;
            
                % IF YES
                if x_lo_ext<x_lo || x_hi_ext>x_hi ||...
                    y_lo_ext<y_lo || y_hi_ext>y_hi;
                %}
                %% Extend data based on error
                sub_is_ext = (E<t_r1_max/2)& (sub_z>=z_loc) & (sub_z<z_hic);
                sub_x_ext = sub_x(sub_is_ext);
                sub_y_ext = sub_y(sub_is_ext);
                sub_z_ext = sub_z(sub_is_ext);
                %X = [sub_x_ext sub_y_ext sub_z_ext]';
                %sub_index_ext = data_index(is_sub_ext);
                sub_index_ext = sub_index(sub_is_ext);
                
                % Plot block extension
                if options_extend;
                    if options_verbose
                        h = figure('color','m');
                    elseif any([options_save_eps, options_save_png, options_save_fig]);
                        h = figure('visible','off','color','w');
                    end
                    hold on
                    scatter3(sub_x(sub_is_ext),sub_y(sub_is_ext),sub_z(sub_is_ext),...
                        10,'b','filled')
                    %  scatter3(sub_x(~sub_is_ext),sub_y(~sub_is_ext),sub_z(~sub_is_ext),...
                    %      10,'r','filled')
                    %{
                    [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,'g',1,0);
                    plot3([z1(1) y1(1)],[z1(2) y1(2)],[z1(3) y1(3)],'-k','linewidth',2)
                    set(Coneh,'facealpha',.2)
                    set(End1h,'facealpha',.2)
                    set(End2h,'facealpha',.2)
                    %}
                    xlabel('x'); ylabel('y'); zlabel('z');
                    campos([0 0 0]);
                    grid on; axis auto
                    daspect([1 1 1]);
                    title('Extend');
                    legend('Original','Extended');
                    if options_save_eps
                        filename_eps_extend = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_extend.eps',...
                            path_eps_detect_trees, xb, yb, bn, eps_n);
                        print(gcf,'-depsc','-opengl',filename_eps_extend)
                        eps_n = eps_n + 1;
                    end
                    if  any([options_save_eps, options_save_png, options_save_fig, options_verbose])
                        delete(h);
                    end
                end
                %}
                %% RANSAC Cylinder fit to extended data
                
                options_verbose_model_cyl = false;
                if ~options_verbose_model_cyl;
                    axis_a = [];
                    axis_e = [];
                    I12ieq = [];
                end
                [cyl_ctr,rcirc,is_inlier,E,~,~,~,cyl_left,cyl_right,return_str,coverage,t_inlier, status] = model_cyl_allcomp(sub_x_ext,sub_y_ext,sub_z_ext,...
                    t_r_adj,t_theta,t_r1_max,t_r1_min, t_taper_min,...
                    t_taper_max,info_sampling,t_min_density,t_min_obj_size,...
                    t_min_pts, t_coverage, sub_axes, options_verbose_model_cyl,a1norm,...
                    axis_a, axis_e, I12ieq); % was 0.5*t_min_density
                
                if return_str.success;
                    r1 = rcirc;
                    z1 = cyl_ctr(:,end);
                    y1 = cyl_ctr(:,1);
                    z_left = cyl_left(:,end);
                    y_left = cyl_left(:,1);
                    z_right = cyl_right(:,end);
                    y_right = cyl_right(:,1);
                    a1 = z1 - y1;
                    t1 = atand((r1(1)-r1(2))/norm(z1-y1)); % Expect small positive values
                else
                    badfit = true;
                    continue
                end
                
                %{
                h = figure;
                hold on
                scatter3(sub_x(sub_is_ext),sub_y(sub_is_ext),sub_z(sub_is_ext),...
                            10,'b','filled')
                scatter3(sub_x(~sub_is_ext),sub_y(~sub_is_ext),sub_z(~sub_is_ext),...
                            10,'r','filled');
                [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,'g',1,0);
                set(Coneh,'facealpha',.2)
                set(End1h,'facealpha',.2)
                set(End2h,'facealpha',.2)
                axis auto;
                daspect([1 1 1]);
                campos([0 0 0]);
                %}
                
                %% Check if fit was outside original block
                %{
                x_cent = (z1(1)+y1(1))/2;
                y_cent = (z1(2)+y1(2))/2;
                [theta,rho] = cart2pol(x_cent,y_cent);
                [x_cent,y_cent] = pol2cart(theta,rho);
                %}
                x_cent = mean(sub_x_ext(is_inlier));
                y_cent = mean(sub_y_ext(is_inlier));
                
                if (x_cent>x_loc-.05 && x_cent<x_hic+0.05 && ...
                        y_cent>y_loc-0.05 && y_cent<y_hic+0.05)
                else
                    if options_outside;
                        if options_verbose
                            h = figure('color','r');
                        elseif any([options_save_eps, options_save_png, options_save_fig]);
                            h = figure('visible','off','color','w');
                        end
                        hold on
                        scatter3(sub_x(sub_is_ext),sub_y(sub_is_ext),sub_z(sub_is_ext),...
                            10,'b','filled')
                        scatter3(sub_x(~sub_is_ext),sub_y(~sub_is_ext),sub_z(~sub_is_ext),...
                            10,'r','filled')
                        %{
                        [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,'g',1,0);
                        plot3([z1(1) y1(1)],[z1(2) y1(2)],[z1(3) y1(3)],'-k','linewidth',2)
                        set(Coneh,'facealpha',.2)
                        set(End1h,'facealpha',.2)
                        set(End2h,'facealpha',.2)
                            %}
                            xlabel('x'); ylabel('y'); zlabel('z');
                            view(0,90);
                            axis auto;
                            xlim([x_lob x_hib]);
                            ylim([y_lob y_hib]);
                            zlim([t_zmin t_zmax]);
                            xtick = round((x_lob:info_tick:x_hib)*4)/4;
                            ytick = round((y_lob:info_tick:y_hib)*4)/4;
                            mid = ceil(numel(xtick)/2);
                            n_xtick = numel(xtick);
                            xticklabel =  cell(n_xtick,1);
                            xticklabel{mid} = sprintf('%0.2f',xtick(mid));
                            n_ytick = numel(ytick);
                            yticklabel =  cell(n_ytick,1);
                            yticklabel{mid} = sprintf('%0.2f',ytick(mid));
                            ztick = round((t_zmin:info_tick:t_zmax)*4)/4;
                            set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
                            set(gca, 'YTick', ytick, 'YTickLabel', yticklabel);
                            set(gca, 'ZTick', ztick, 'ZTickLabel', sprintf('%0.2f|', ztick));
                            scatter3(x_cent,y_cent,t_zmax,50,'k','filled')
                            rectangle('position',[x_loc,y_loc,t_r,t_r]) % might get messed up from adjusting t_r for dem height
                            grid on
                            daspect([1 1 1]);
                            title('Failed: Outside');
                            if options_save_eps
                                filename_eps_outside = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_outside.eps',...
                                    path_eps_detect_trees, xb, yb, bn, eps_n);
                                print(gcf,'-depsc','-opengl',filename_eps_outside)
                                eps_n = eps_n + 1;
                            end
                            if  any([options_save_eps, options_save_png, options_save_fig, options_verbose])
                                delete(h);
                            end
                    end
                    badfit=true;
                    continue
                end
                %}
                
                %% Update data_is_inlier
                %is_remove = (E<t_remove);
                %info_point_spacing = norm(mean([sub_x,sub_y,sub_z]))*t_sampling;
                is_remove = E<3*t_inlier;%E<2*mean(r1);%5*info_point_spacing; %was is_inlier
                data_is_removed(sub_index_ext(is_remove)) = true; % Points excluded in future iterations
                index_inlier = data_index(sub_index_ext(is_inlier)); % Inlier points for coloring
                
                %{
                h1 = figure;
                is_only_remove = is_remove & ~is_inlier;
                is_only_inlier = is_inlier;
                is_only_other = ~(is_only_remove | is_only_inlier);
                hold on
                [~,~,~] = Cone(z1,y1,[r1(2) r1(1)],30,[.5 .5 .5],0,0);
                alpha(0.5);
                scatter3(sub_x_ext(is_only_remove), sub_y_ext(is_only_remove), sub_z_ext(is_only_remove),20,'r','filled');
                scatter3(sub_x_ext(is_only_inlier), sub_y_ext(is_only_inlier), sub_z_ext(is_only_inlier),20,'g','filled');
                scatter3(sub_x_ext(is_only_other), sub_y_ext(is_only_other), sub_z_ext(is_only_other),20,[.5 .5 .5],'filled');
                daspect([1 1 1]); xlabel('x'); ylabel('y'); zlabel('z'); campos([0 0 0]);
                legend('Cone','Removed', 'Inliers','Remaining');
                pause
                delete(h1);
                %}
                %% Plot the final cylinder results
                if options_final;
                    if options_verbose
                        if return_str.coverage;
                            h = figure('color','r');
                        else
                            h = figure('color','g');
                        end
                    elseif any([options_save_eps, options_save_png, options_save_fig]);
                        h = figure('visible','off','color','w');
                    end
                    sub_ixcenter = find(sub_is_ext);
                    
                    inlier_center = sub_ixcenter(is_inlier);
                    outlier_center = sub_ixcenter(~is_inlier);
                    is = false(numel(sub_x),1);
                    is_center_inlier = is;
                    is_center_inlier(inlier_center) = true;
                    is_center_outlier = is;
                    is_center_outlier(outlier_center) = true;
                    scatter3(sub_x(is_center_inlier), sub_y(is_center_inlier), sub_z(is_center_inlier), 20, 'b', 'filled');
                    hold on
                    scatter3(sub_x(is_center_outlier), sub_y(is_center_outlier), sub_z(is_center_outlier), 10, 'r', 'filled');
                    scatter3(sub_x(~sub_iscenter), sub_y(~sub_iscenter), sub_z(~sub_iscenter), 30, [.5 .5 .5]);
                    campos([0 0 0]);
                    daspect([1 1 1]);
                    xlim([x_lob x_hib]);
                    ylim([y_lob y_hib]);
                    zlim([z_lob z_hib]);
                    xtick = round((x_lob:info_tick:x_hib)*4)/4;
                    ytick = round((y_lob:info_tick:y_hib)*4)/4;
                    ztick = round((z_lob:info_tick:z_hib)*4)/4;
                    mid = ceil(numel(xtick)/2);
                    n_xtick = numel(xtick);
                    xticklabel =  cell(n_xtick,1);
                    xticklabel{mid} = sprintf('%0.2f',xtick(mid));
                    n_ytick = numel(ytick);
                    yticklabel =  cell(n_ytick,1);
                    yticklabel{mid} = sprintf('%0.2f',ytick(mid));
                    n_ztick = numel(ytick);
                    zticklabel =  cell(n_ztick,1);
                    zticklabel{mid} = sprintf('%0.2f',ztick(mid));
                    set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
                    set(gca, 'YTick', ytick, 'YTickLabel', yticklabel);
                    set(gca, 'ZTick', ztick, 'ZTickLabel', zticklabel);
                    set(gca,'zticklabelmode','auto')
                    set(gca,'yticklabelmode','auto')
                    set(gca,'xticklabelmode','auto')
                    [Coneh,End1h,End2h] = Cone(z1,y1,[r1(2) r1(1)],30,cylcolor,1,0);
                    set(Coneh,'facealpha',.4)
                    set(End1h,'facealpha',.4)
                    set(End2h,'facealpha',.4)
                    title(title_coverage);
                    xlabel('x'); ylabel('y'); zlabel('z');
                    set(gcf, 'PaperPositionMode', 'auto')
                    if options_save_eps;
                        if coverage < t_coverage;
                            filename_eps_final_reject = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_final_reject.eps',...
                                path_eps_detect_trees, xb, yb, bn, eps_n);
                            print(gcf,'-depsc','-opengl',filename_eps_final_reject)
                        else
                            filename_eps_final_accept = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_final_accept.eps',...
                                path_eps_detect_trees, xb, yb, bn, eps_n);
                            print(gcf,'-depsc','-opengl',filename_eps_final_accept)
                        end
                        eps_n = eps_n + 1;
                    end
                    if  any([options_save_eps, options_save_png, options_save_fig, options_verbose])
                        delete(h);
                    end
                end
                
                %% Plot tree axis on Andrieu image
                %
                if options_andrieu;
                    axis_c = 1:numel(axis_a);
                    axis_r = 1:numel(axis_e);
                    [sub_a, sub_e,~] = cart2sph(sub_x_ext,sub_y_ext,sub_z_ext);
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
                        if is_inlier(p);
                            I12ieq_c(sub_row(p),sub_col(p),1) = val/2;
                            I12ieq_c(sub_row(p),sub_col(p),2) = val/2;
                        else
                            I12ieq_c(sub_row(p),sub_col(p),2) = val/2;
                            I12ieq_c(sub_row(p),sub_col(p),3) = val/2;
                        end
                    end
                    
                    I12ieq_c(:,:,1) = fliplr(I12ieq_c(:,:,1));
                    I12ieq_c(:,:,2) = fliplr(I12ieq_c(:,:,2));
                    I12ieq_c(:,:,3) = fliplr(I12ieq_c(:,:,3));
                    
                    if options_verbose
                        h = figure('color','w');
                        imshow(I12ieq_c)
                        colormap('gray');
                        axis off; axis image
                        set(gca,'Units','normalized','Position',[0 0 1 1]);
                        [n_row, n_col,~] = size(I12ieq);
                        set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
                    end
                    if options_save_png;
                        if coverage < t_coverage;
                            filename_png_andrieu_reject = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_andrieu_reject.png',...
                                path_png_detect_trees, xb, yb, bn, png_n);
                            imwrite(I12ieq_c,filename_png_andrieu_reject,'png');
                        else
                            filename_png_andrieu_accept = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_andrieu_accept.png',...
                                path_png_detect_trees, xb, yb, bn, png_n);
                            imwrite(I12ieq_c,filename_png_andrieu_accept,'png');
                        end
                        png_n = png_n + 1;
                    end
                    if  options_verbose;
                        delete(h);
                    end
                end
                
                %% Save parameters
                if return_str.success
                    seg_xb(s) = xb;
                    seg_yb(s) = yb;
                    seg_zb(s) = zb;
                    seg_z(s,:) = z1;
                    seg_y(s,:) = y1;
                    seg_r(s,:) = r1;
                    seg_iter(s) = s;
                    seg_index{s} = index_inlier;
                    seg_fill(s) = coverage;
                    seg_taper(s) = t1;
                    seg_min_obj_size(s) = t_min_obj_size;
                    seg_z_left(s,:) = z_left;
                    seg_y_left(s,:) = y_left;
                    seg_z_right(s,:) = z_right;
                    seg_y_right(s,:) = y_right;
                    seg_status{s} = status;
                    
                    
                    % Update counter
                    s = s  + 1;
                    bn = bn + 1;
                    if options_print
                        fprintf('\nBlock (%g,%g,%g) continued to completion \n',xb,yb,zb);
                    end
                    %cylcolor = 'r';
                    %badfit = true;
                    %continue
                end
                block_iter = block_iter + 1;
            end %badfit
            %}
            
        end %zb
    end %yb
end %xb

%% Update segment arrays
is_seg = (seg_iter~=0);
seg_xb = seg_xb(is_seg);
seg_yb = seg_yb(is_seg);
%seg_zb = seg_zb(is_seg);
seg_z = seg_z(is_seg,:);
seg_y = seg_y(is_seg,:);
seg_r = seg_r(is_seg,:);
seg_iter = seg_iter(is_seg);
seg_index = seg_index(~cellfun('isempty',seg_index));
seg_fill = seg_fill(is_seg);
seg_taper = seg_taper(is_seg);
seg_min_obj_size = seg_min_obj_size(is_seg);
seg_z_left = seg_z_left(is_seg,:);
seg_y_left = seg_y_left(is_seg,:);
seg_z_right = seg_z_right(is_seg,:);
seg_y_right = seg_y_right(is_seg,:);
seg_status = seg_status(~cellfun('isempty',seg_status));
end % function

