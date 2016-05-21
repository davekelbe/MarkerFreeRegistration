function  [seg_xb, seg_yb,seg_zb, seg_z, seg_y, seg_r, seg_iter, seg_index, seg_fill] = detect_trees(data_x,...
    data_y,data_z,data_z0,...
    I12ieq,axis_a,axis_e, ...
    t_r, t_zmin0, t_zmax0, path_top,...
    t_minpts,t_error_line,t_error_circ,t_error_side,t_theta,t_r1_max,t_Eext,t_remove,t_coverage,...
    t_sampling,t_density_cutoff,info_z_step,t_kasa)
%DETECT_TREES Wrapper for model_cyl


%% Initialization
%{
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
options_print = true;
options_verbose = true;
options_save_eps = false;
options_save_fig = false;
options_save_png = false;
options_point_subset = false;
options_andrieu1 = false;
options_radius = true;
options_vertical = false;
options_extend = false;
options_outside = true;
%options_bins = true;
options_vertical2 = true;
options_final = true;
options_andrieu = false;
%}

% Additional info
info_slash = '\';
info_tick = t_r/4; % Set three tick marks on plot

% Set paths
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
path_fig = sprintf('%s%s%s',path_top,'fig',info_slash);
path_eps = sprintf('%s%s%s',path_top,'eps',info_slash);
path_png = sprintf('%s%s%s',path_top,'png',info_slash);
path_eps_detect_trees = sprintf('%sdetect_trees%s',path_eps,info_slash);
path_png_detect_trees = sprintf('%sdetect_trees%s',path_png,info_slash);

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
seg_r = zeros(info_maxiter,1); % radius
seg_iter = zeros(info_maxiter,1); % main iteration
seg_index = cell(info_maxiter,1); % index to inlier points
seg_fill = zeros(info_maxiter,1); % Fill percentage based on sampling

% Other preparation
n_pts = numel(data_x);
data_infxn_is_inlier = false(n_pts,1);
data_infxn_index = (1:n_pts)'; %Index to data (within function)
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
    for yb = 9:n_y_block;
        for zb = 6:n_z_block;
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
            
            fprintf('\nBlock (%g,%g,%g) from [%2.0f,%2.0f);[%2.0f,%2.0f);[%2.0f,%2.0f)\n',xb,yb,zb,x_loc,x_hic,y_loc,y_hic,z_loc,z_hic);
            
            %% Continue in block
            
            % Continue until all potential trees are exhausted'
            badfit = false;
            bn = 1; %Block number (iterated)
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
                    data_z0>=t_zmin0&data_z0<t_zmax0&~data_infxn_is_inlier;
                n_sub = sum(data_is_sub);
                
                % Check if there are enough points
                if ~options_point_subset % If images need to be plotted, wait to continue until after plot
                    if n_sub < t_minpts;
                        if options_print
                            fprintf('\nBlock (%g,%g,%g) aborted due to insufficient # points\n',xb,yb,zb);
                        end
                        badfit = true;
                        continue
                    end
                end
                
                % Subset data
                sub_x = data_x(data_is_sub);
                sub_y = data_y(data_is_sub);
                sub_z = data_z(data_is_sub);
                sub_infxn_index = data_infxn_index(data_is_sub); % Index to points in data (in function)
                
                % Find center of subset
                sub_iscenter = sub_x>=x_loc & sub_x<x_hic&...
                    sub_y>=y_loc & sub_y<y_hic&...
                    sub_z>=z_loc & sub_z<z_hic;
                n_center = sum(sub_iscenter);
                
                % Determine axes
                dem_val = mean(sub_z-data_z0(data_is_sub));
                t_zmin = t_zmin0 + dem_val;
                t_zmax = t_zmax0 + dem_val;
                block_axes = [x_lob x_hib y_lob y_hib z_lob z_hib];
                
                % Plot point subset
                if options_point_subset && n_center>10 ;%&& n_center<t_minpts;
                    if options_verbose
                        if n_center<t_minpts;
                            h = figure('color','r');
                        else
                            h = figure('color','g');
                        end
                    elseif any([options_save_eps, options_save_png, options_save_fig]);
                        h = figure('visible','off','color','w');
                    end
                    scatter3(sub_x(sub_iscenter), sub_y(sub_iscenter), sub_z(sub_iscenter), 10, 'b', 'filled');
                    hold on
                    scatter3(sub_x(~sub_iscenter), sub_y(~sub_iscenter), sub_z(~sub_iscenter), 10, 'r', 'filled');
                    campos([0 0 0])
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
                    if n_center < t_minpts;
                        pf = 'FAIL';
                    else
                        pf = 'PASS';
                    end
                    title(sprintf('(%g, %g, %g); Subset: %s: %g points',xb,yb,zb, pf, n_center));
                    xlabel('x'); ylabel('y'); zlabel('z');
                    legend('center','exterior','location','east')
                    daspect([1 1 1]);
                    if options_save_eps;
                        if n_center < t_minpts;
                            filename_eps_point_subset_reject = sprintf('%sx%02.0f-y%02.0f-z%02.0f_b%01.0f_i%01.0f_pointsubset-reject.eps',...
                                path_eps_detect_trees, xb, yb, zb, bn, eps_n);
                            print(gcf,'-depsc','-painters',filename_eps_point_subset_reject)
                        else
                            filename_eps_point_subset_accept = sprintf('%sx%02.0f-y%02.0f-z%02.02f_b%01.0f_i%01.0f_pointsubset-accept.eps',...
                                path_eps_detect_trees, xb, yb, zb, bn, eps_n);
                            print(gcf,'-depsc','-painters',filename_eps_point_subset_accept)
                        end
                        eps_n = eps_n + 1;
                    end
                    if  any([options_save_eps, options_save_png, options_save_fig, options_verbose])
                        delete(h);
                    end
                end
                %% Andrieu Image of subset
                
                % Color subset on Andrieu image
                %
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
                
                if sum(sub_iscenter) < t_minpts;
                    if options_print;
                        fprintf('\nBlock (%g,%g,%g) aborted due to insufficient # points\n',xb,yb,zb);
                    end
                    badfit = true;
                    continue
                end
                
                %% Initial RANSAC Cylinder fit to center points
                
                [cyl_ctr,rcirc,~,~,M,~,~] = model_cyl(sub_x(sub_iscenter),sub_y(sub_iscenter),sub_z(sub_iscenter),...
                    t_error_line,t_error_circ,t_error_side,t_theta,t_sampling,t_density_cutoff,info_z_step,t_kasa);
                % Convert to coarse estimate for consistency with code
                if ~isempty(rcirc);
                    r1 = mean(rcirc);
                    z1 = cyl_ctr(:,end);
                    y1 = cyl_ctr(:,1);
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
                
                %% Check if radius is too large
                if r1>t_r1_max
                    if options_print
                        fprintf('\nBlock (%g,%g,%g) aborted due to large radius \n',xb,yb,zb);
                    end
                    badfit = true;
                    if options_radius;
                        if options_verbose
                            h = figure('color','r');
                        elseif any([options_save_eps, options_save_png, options_save_fig]);
                            h = figure('visible','off','color','w');
                        end
                        scatter3(sub_x(sub_iscenter), sub_y(sub_iscenter), sub_z(sub_iscenter), 10, 'b', 'filled');
                        hold on
                        scatter3(sub_x(~sub_iscenter), sub_y(~sub_iscenter), sub_z(~sub_iscenter), 10, 'r', 'filled');
                        [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,cylcolor,1,0);
                        set(Coneh,'facealpha',.2)
                        set(End1h,'facealpha',.2)
                        set(End2h,'facealpha',.2)
                        campos([0 0 0])
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
                        n_ztick = numel(ztick);
                        zticklabel =  cell(n_ztick,1);
                        zticklabel{mid} = sprintf('%0.2f',ztick(mid));
                        set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
                        set(gca, 'YTick', ytick, 'YTickLabel', yticklabel);
                        set(gca, 'ZTick', ytick, 'ZTickLabel', zticklabel);
                        if r1>t_r1_max
                            pf = 'FAIL';
                        else
                            pf = 'PASS';
                        end
                        title(sprintf('(%g, %g, %g); Radius: %s: %2.1f m',xb,yb,zb,pf,r1));
                        xlabel('x'); ylabel('y'); zlabel('z');
                        daspect([1 1 1]);
                        % saveas(gca,filename_eps_cylfit,'epsc');
                        if options_save_eps
                            filename_eps_vertical_reject = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_radius_reject.eps',...
                                path_eps_detect_trees, xb, yb, bn, eps_n);
                            print(gcf,'-depsc','-painters',filename_eps_vertical_reject)
                            eps_n = eps_n + 1;
                        end
                        if  any([options_save_eps, options_verbose])
                            delete(h);
                        end
                    end
                    continue
                end
                
                %% Check if the fit was vertical
                if isnan(r1);
                    if options_print
                        fprintf('\nBlock (%g,%g,%g) aborted due to non-vertical normal \n',xb,yb,zb);
                    end
                    badfit = true;
                    if options_vertical; % Plot the results
                        theta = acosd(a1norm(3));
                        if options_verbose
                            h = figure('color','r');
                        elseif any([options_save_eps, options_save_png, options_save_fig]);
                            h = figure('visible','off','color','w');
                        end
                        scatter3(sub_x(sub_iscenter), sub_y(sub_iscenter), sub_z(sub_iscenter), 10, 'b', 'filled');
                        hold on
                        scatter3(sub_x(~sub_iscenter), sub_y(~sub_iscenter), sub_z(~sub_iscenter), 10, 'r', 'filled');
                        plot3([z1(1) y1(1)],[z1(2) y1(2)],[z1(3) y1(3)],'-k','LineWidth',2)
                        campos([0 0 0])
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
                        %ztick = round((t_zmin:tick:t_zmax)*4)/4;
                        n_ztick = numel(ztick);
                        zticklabel =  cell(n_ztick,1);
                        zticklabel{mid} = sprintf('%0.2f',ztick(mid));
                        %ztick = round((t_zmin:tick:t_zmax)*4)/4;
                        set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
                        set(gca, 'YTick', ytick, 'YTickLabel', yticklabel);
                        set(gca, 'ZTick', ytick, 'ZTickLabel', zticklabel);
                        %set(gca, 'ZTick', ztick, 'ZTickLabel', sprintf('%0.2f|', ztick));
                        if isempty(r1)
                            pf = 'FAIL';
                        else
                            pf = 'PASS';
                        end
                        title(sprintf('(%g, %g); Angle: %s: %2.1f deg',xb,yb,pf,theta));
                        xlabel('x'); ylabel('y'); zlabel('z');
                        daspect([1 1 1]);
                        % saveas(gca,filename_eps_cylfit,'epsc');
                        if options_save_eps
                            filename_eps_vertical_reject = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_vertical_reject.eps',...
                                path_eps_detect_trees, xb, yb, bn, eps_n);
                            print(gcf,'-depsc','-painters',filename_eps_vertical_reject)
                            eps_n = eps_n + 1;
                        end
                        if  any([options_save_eps, options_verbose])
                            delete(h);
                        end
                    end
                    continue
                end
                
                %% Recompute error for all points
                X = [sub_x sub_y sub_z]';
                Xnew = M*X;
                znew = M*z1;
                X2 = [Xnew(2,:); Xnew(3,:)];
                X2 = X2 - repmat([znew(2); znew(3)],[1, n_sub]);
                E =  abs(sqrt(sum(X2.^2,1))-r1)';
                
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
                sub_is_ext = (E<t_Eext)& (sub_z>=z_loc) & (sub_z<z_hic);
                sub_x_ext = sub_x(sub_is_ext);
                sub_y_ext = sub_y(sub_is_ext);
                sub_z_ext = sub_z(sub_is_ext);
                %X = [sub_x_ext sub_y_ext sub_z_ext]';
                %sub_index_ext = data_index(is_sub_ext);
                sub_infxn_index_ext = sub_infxn_index(sub_is_ext);
                
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
                        campos([0 0 0]);
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
                        grid on
                        daspect([1 1 1]);
                        title('Extend');
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
                
                [cyl_ctr,rcirc,is_inlier,E,~,~,~] = model_cyl(sub_x(sub_iscenter),sub_y(sub_iscenter),sub_z(sub_iscenter),...
                    t_error_line,t_error_circ,t_error_side,t_theta,t_sampling,t_density_cutoff,info_z_step,t_kasa);
                if ~isempty(rcirc);
                    r1 = mean(rcirc);
                    z1 = cyl_ctr(:,end);
                    y1 = cyl_ctr(:,1);
                    a1 = z1 - y1;
                end
                
                % Check if the fit was vertical
                if isempty(rcirc);
                    if options_print
                        fprintf('\nBlock (%g,%g,%g) aborted due to non-vertical normal \n',xb,yb,zb);
                    end
                    badfit = true;
                    continue
                end
                
                %{
                h = figure;
                hold on
                scatter3(sub_x(sub_is_ext),sub_y(sub_is_ext),sub_z(sub_is_ext),...
                            10,'b','filled')
                %scatter3(sub_x(~sub_is_ext),sub_y(~sub_is_ext),sub_z(~sub_is_ext),...
                %            10,'r','filled')
                %
                %{
                [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,'r',1,0);
                set(Coneh,'facealpha',.2)
                set(End1h,'facealpha',.2)
                set(End2h,'facealpha',.2)
                %}
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
                            rectangle('position',[x_loc,y_loc,t_r,t_r])
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
                
                %% Check if the fit was vertical
                if isempty(r1);
                    if options_print
                        fprintf('\nBlock (%g,%g) aborted due to non-vertical normal \n',xb,yb);
                    end
                    badfit = true;
                    if options_vertical2; % Plot the results
                        a1 = z1 - y1;
                        a1norm = a1/norm(a1);
                        if a1norm(3)<0;
                            a1norm = -a1norm;
                        end
                        theta = acosd(a1norm(3));
                        if options_verbose
                            h = figure('color','r');
                        elseif any([options_save_eps, options_save_png, options_save_fig]);
                            h = figure('visible','off','color','w');
                        end
                        scatter3(sub_x_ext, sub_y_ext, sub_z_ext, 10, E, 'filled');
                        hold on
                        colorbar
                        plot3([z1(1) y1(1)],[z1(2) y1(2)],[z1(3) y1(3)],'-k','LineWidth',2)
                        campos([0 0 0])
                        xlim_new = get(gca, 'xlim');
                        ylim_new = get(gca, 'ylim');
                        zlim_new = get(gca, 'zlim');
                        x_lo_new = xlim_new(1);
                        x_hi_new = xlim_new(2);
                        y_lo_new = ylim_new(1);
                        y_hi_new = ylim_new(2);
                        t_zmin_new = zlim_new(1);
                        t_zmax_new = zlim_new(2);
                        xtick_new = round((x_lo_new:info_tick:x_hi_new)*4)/4;
                        ytick_new = round((y_lo_new:info_tick:y_hi_new)*4)/4;
                        n_xtick_new = numel(xtick_new);
                        xticklabel =  cell(n_xtick_new,1);
                        midx = ceil(numel(xtick_new)/2);
                        xticklabel{midx} = sprintf('%0.2f',xtick(midx));
                        n_ytick_new = numel(ytick);
                        yticklabel =  cell(n_ytick_new,1);
                        midy = ceil(numel(ytick_new)/2);
                        yticklabel{midy} = sprintf('%0.2f',ytick(midy));
                        ztick_new = round((t_zmin_new:info_tick:t_zmax_new)*4)/4;
                        set(gca, 'XTick', xtick_new, 'XTickLabel', xticklabel);
                        set(gca, 'YTick', ytick_new, 'YTickLabel', yticklabel);
                        set(gca, 'ZTick', ztick_new, 'ZTickLabel', sprintf('%0.2f|', ztick_new));
                        title(sprintf('(%g, %g); Angle = %2.1f deg',xb,yb,theta));
                        xlabel('x'); ylabel('y'); zlabel('z');
                        daspect([1 1 1]);
                        if options_save_eps
                            filename_eps_vertical2_reject = sprintf('%sx%02.0f-y%02.0f_b%01.0f_i%01.0f_vertical2.eps',...
                                path_eps_detect_trees, xb, yb, bn, eps_n);
                            print(gcf,'-depsc','-painters',filename_eps_vertical2_reject)
                            eps_n = eps_n + 1;
                        end
                        if  any([options_save_eps, options_save_png, options_save_fig, options_verbose])
                            delete(h);
                        end
                    end
                    continue
                end
                
                % end % End extension
                
                %% Update data_is_inlier
                is_remove = (E<t_remove);
                data_infxn_is_inlier(sub_infxn_index_ext(is_remove)) = true;
                index_inlier = data_infxn_index(sub_infxn_index_ext(is_remove));
                
                %% Check confidence metric
                bole_area = 2*r1*sum(sqrt((z1-y1).^2));
                bole_range = (sqrt(z1'*z1)+sqrt(y1'*y1))/2;
                pt_spacing = t_sampling*bole_range;
                n_pts_exp = bole_area/pt_spacing^2;
                coverage = sum(is_inlier)/n_pts_exp;
                
                if coverage < t_coverage;
                    pf = 'FAIL';
                else
                    pf = 'PASS';
                end
                
                title_coverage = sprintf('Coverage: (%g,%g,%g) %s %2.0f%%',xb,yb,zb,pf, 100*coverage);
                
                %% Plot the final cylinder results
                if options_final;
                    if options_verbose
                        if coverage < t_coverage;
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
                    [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,cylcolor,1,0);
                    set(Coneh,'facealpha',.4)
                    set(End1h,'facealpha',.4)
                    set(End2h,'facealpha',.4)
                    title(title_coverage);
                    xlabel('x'); ylabel('y'); zlabel('z');
                    %{
                    hold on
                    scatter3(sub_x_ext(is_remove),sub_y_ext(is_remove),sub_z_ext(is_remove),...
                        10,'b','filled')
                    scatter3(sub_x_ext(~is_remove),sub_y_ext(~is_remove),sub_z_ext(~is_remove),...
                        10,'r','filled')
                    %scatter3(sub_x(is_inlier),sub_y(is_inlier),sub_z(is_inlier),...
                    %    10,'b','filled')
                    %scatter3(sub_x(~is_inlier),sub_y(~is_inlier),sub_z(~is_inlier),...
                    %    10,'r','filled')
                    [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,cylcolor,1,0);
                    set(Coneh,'facealpha',.2)
                    set(End1h,'facealpha',.2)
                    set(End2h,'facealpha',.2)
                    xlim_new = get(gca, 'xlim');
                    ylim_new = get(gca, 'ylim');
                    zlim_new = get(gca, 'zlim');
                    x_lo_new = xlim_new(1);
                    x_hi_new = xlim_new(2);
                    y_lo_new = ylim_new(1);
                    y_hi_new = ylim_new(2);
                    z_lo_new = zlim_new(1);
                    z_hi_new = zlim_new(2);
                    xtick = round((x_lo_new:tick:x_hi_new)*4)/4;
                    ytick = round((y_lo_new:tick:y_hi_new)*4)/4;
                    ztick = round((z_lo_new:tick:z_hi_new)*4)/4;
                    %ztick = round((t_zmin:tick:t_zmax)*4)/4;
                    n_xtick = numel(xtick);
                    xticklabel =  cell(n_xtick,1);
                    if numel(xtick)>1
                        xticklabel{2} = sprintf('%0.2f',xtick(2));
                    else
                        xticklabel{1} = sprintf('%0.2f',xtick(1));
                    end
                    n_ytick = numel(ytick);
                    yticklabel =  cell(n_ytick,1);
                    if numel(ytick)>1
                        yticklabel{2} = sprintf('%0.2f',ytick(2));
                    else
                        yticklabel{1} = sprintf('%0.2f',ytick(1));
                    end
                    n_ztick = numel(ztick);
                    zticklabel =  cell(n_ztick,1);
                    if numel(ztick)>1
                        zticklabel{2} = sprintf('%0.2f',ztick(2));
                    else
                        zticklabel{1} = sprintf('%0.2f',ztick(1));
                    end
                    set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
                    set(gca, 'YTick', ytick, 'YTickLabel', yticklabel);
                    set(gca, 'ZTick', ytick, 'ZTickLabel', zticklabel);
                    %set(gca, 'ZTick', ztick, 'ZTickLabel', sprintf('%0.2f|', ztick));
                    xlabel('x'); ylabel('y'); zlabel('z');
                    campos([0 0 0]);
                    daspect([1 1 1]); grid on
                    %}
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
                seg_xb(s) = xb;
                seg_yb(s) = yb;
                seg_zb(s) = zb;
                seg_z(s,:) = z1;
                seg_y(s,:) = y1;
                seg_r(s) = r1;
                seg_iter(s) = s;
                seg_index{s} = index_inlier;
                seg_fill(s) = coverage;
                % Update counter
                s = s  + 1;
                bn = bn + 1;
                
                if coverage < t_coverage;
                    if options_print
                        fprintf('\nBlock (%g,%g,%g) continued to completion \n',xb,yb,zb);
                    end
                    cylcolor = 'r';
                    badfit = true;
                    continue
                end
            end %badfit
            %}
            
        end %zb 
    end %yb
end %xb

%% Update segment arrays 
is_seg = (seg_r~=0);
seg_xb = seg_xb(is_seg);
seg_yb = seg_yb(is_seg);
seg_z = seg_z(is_seg,:);
seg_y = seg_y(is_seg,:);
seg_r = seg_r(is_seg);
seg_iter = seg_iter(is_seg);
seg_index = seg_index(~cellfun('isempty',seg_index));
seg_fill = seg_fill(is_seg);

end % function

