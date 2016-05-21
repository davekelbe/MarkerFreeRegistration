function [  ] = kelbe_registration_combine_dijkstraposewcs_param_analyses( aux )
% Only load valid plots
%info_valid_plot = {'13','15'};
%      info_valid_plot = cell(1,25);
%      for s = 1:5;
%          info_valid_plot{s} = sprintf('%02.0f', s);
%      end
%      info_site = 31;
options_showfig = false;
info_valid_plot = aux.info_valid_plot;
info_site = aux.info_site;
path_mat = aux.path_mat;
    aux.P_LCS = aux.P_LCSn;
    aux.P_rad = aux.P_radn;
tree = 1;
%% KELBE REGISTRATION
% Clear without breakpoints
% Default settings
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = false;
path_save = 'Z:\Desktop\RegistrationOutput\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_matvar = 'D:\Users\djk2312\Documents\matvar\';
options_loadvar =false;
%if ~options_loadvar;
% Initialize program-level options and paths
% Outputs
% options               options
% paths                 paths
clear D ctr d filepath_tree plot t tree
%% Parameters
t_rad = 0.2; %0.06^2; % Trees can be matched if their radius is within threshold
t_coll = 0.05; % Collinnearity threshold: don't consider triangles that are collinear
t_eig_error = 1e1;
t_RANSAC_xyz = 0.4^2;
t_RANSAC_nsearch = 100;
t_flagr = 20; % deg
t_flagt = 1; % m


% Initialize program-level options and paths
% Outputs
% t_*               thresholds
%% Load points
% Input to function are stem maps derived from lidar.m
fprintf('\nLoad points\n');

% Filename Lookup for Harvard Data
info_exp = aux.info_exp;
info_suffix = aux.info_suffix;
info_slash = '\';
path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);

n_S = numel(info_valid_plot);
P_LCS = aux.P_LCS;
[~,n_tree] = cellfun(@size, P_LCS);
P_rad = aux.P_rad;
P_plot = 1:n_S;
% path_ply = cell(n_S,1);
% filepath_ply = cell(n_S,1);
% for d = 1:numel(D);
%     path_ply{ctr} = sprintf('%s%s%sply%s',path_site,info_plot,info_slash,info_slash);
%     filepath_ply{ctr} = sprintf('%spoints_full_%03.0f-%02.0f.ply', path_ply{ctr}, info_site, plot);
% end

isvalid = (P_plot~=0);
P_LCS = P_LCS(isvalid);
P_rad = P_rad(isvalid);
P_plot = P_plot(isvalid);
n_S = sum(isvalid);

%[P_LCS, P_rad, P_plot, n_tree, n_S ] = generate_test_data;
%{
i_xmin = min(cellfun(@(x) min(x(1,:)),P_LCS));
i_xmax = max(cellfun(@(x) max(x(1,:)),P_LCS));
i_ymin = min(cellfun(@(x) min(x(2,:)),P_LCS));
i_ymax = max(cellfun(@(x) max(x(2,:)),P_LCS));
i_zmin = min(cellfun(@(x) min(x(3,:)),P_LCS));
i_zmax = max(cellfun(@(x) max(x(3,:)),P_LCS));
options_axesval = [i_xmin i_xmax i_ymin i_ymax i_zmin i_zmax];
%}
% Colormap for sensors
P_color = jet(n_S);

% Individual camera views
%{
if false;%options_verbose && options_imagepoints;
    for s = 1:n_S;
        clear legend_str
        figure
        hold on
        plot3(0,0,0,'^k','markersize',10,...
            'markerfacecolor',P_color(s,:));
        hdummy = plot3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),'ok','markersize',5,...
            'markerfacecolor',P_color(s,:));
        set(hdummy, 'visible', 'off');
        for t = 1:numel(P_rad{s});
            h = filledCircle([P_LCS{s}(1,t); P_LCS{s}(2,t)]',P_rad{s}(t),1000,P_color(s,:));
        end
        %scatter3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),30,...
        %    color_P_index(truth_P_index{s},:),'filled');
        %axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
        axis auto
        axisval = axis;
        %set(gca, 'xtick',
        xlabel('x Position relative to plot center [m]');
        ylabel('y Position relative to plot center [m]');
        zlabel('z Position relative to plot center [m]');
        view(0,90);
        grid on
        %titlestr = sprintf('Scan %g',P_plot(s));
        %title(titlestr);
        legend_str{1} = sprintf('Scanner %g',P_plot(s));
        legend_str{2} = sprintf('Stem map in LCS %g', P_plot(s));
        legend(legend_str,'location','northeast');
        legend boxoff       % Hides the legend's axes (legend border and background)
        set(gca, 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        % filepath_save = sprintf('%sLCS_%02.0f.eps',path_save, P_plot(s));
        % saveas(gcf,filepath_save,'psc2')
    end
end
%}
%}
% Outputs
% P_LCS                 points in local coordinate system (cell)
% P_rad                 radius of points
% P_plot                plot number
% P_n                   number of points
% i_xmin                minimum x value of data
% i_zmax                maximum z value of data
% color                 colormap for scans
% info                  program level information
% options               program level options
% path                  paths
% n                     number of given variable
clear D ctr d filepath_tree plot t tree i_xmax i_xmin i_ymin i_ymax i_zmin i_zmax
%if false;
%% Determine Pose Consistency in WCS
% Load variables
filepath_match12_RMSE = sprintf('%s%s',path_mat, 'match12_RMSE.mat');
load(filepath_match12_RMSE);
filepath_match12_R = sprintf('%s%s',path_mat, 'match12_R.mat');
load(filepath_match12_R);
filepath_match12_t = sprintf('%s%s',path_mat, 'match12_t.mat');
load(filepath_match12_t);
filepath_match1_R = sprintf('%s%s',path_mat, 'match1_R.mat');
load(filepath_match1_R);
filepath_match_R = sprintf('%s%s',path_mat, 'match_R.mat');
if exist(filepath_match_R,'file')
load(filepath_match_R);
else
    match_R = match1_R;
end
filepath_match1_t = sprintf('%s%s',path_mat, 'match1_t.mat');
load(filepath_match1_t);
filepath_match_t = sprintf('%s%s',path_mat, 'match_t.mat');
if exist(filepath_match_t,'file')
load(filepath_match_t);
else
    match_t = match1_t;
end
filepath_match2_R = sprintf('%s%s',path_mat, 'match2_R.mat');
load(filepath_match2_R);
filepath_match2_t = sprintf('%s%s',path_mat, 'match2_t.mat');
load(filepath_match2_t);
%match_R = match1_R; % Fixed - now load real reestimate
%match_t = match1_t; %
foo = 1;
%% Check if abort 
if sum(isnan(match12_RMSE(:))) == sum(sum(~eye(size(match12_RMSE))));
    return
end
if size(match12_RMSE,1)~=n_S;
    return
end
%% Remove very bad matches

match12_RMSE(match12_RMSE>1) = nan;
%% Effective (WCS) R and t using Dijkstra's
%
% Shortest path to WCS
% Weighted undirected adjacency matrix
G_path = cell(n_S,n_S);
G_dist = cell(n_S,n_S);
%G_pred = cell(n_S,n_S);
isbad = isnan(match12_RMSE);
match12_RMSE(isbad) = inf;
match12h_RMSE = match12_RMSE./2;
for i = 1:n_S
    for j = 1:n_S;
      %  fprintf('\n%d,%d\n',i,j);
      %  if i==1&&j==22&&info_site==6;
      %      foo = 1;
      %  end
        [G_dist{i,j},G_path{i,j},~] = graphshortestpath(sparse(match12_RMSE'),j,i);
    end
end

% 4x4 matrics 
match1_R4 = cellfun(@(x,y) cat(2,x,y),match_R,match_t, 'uniformoutput',false);
match1_R4 = cellfun(@(x) cat(1,x,[0 0 0 1]),match1_R4, 'uniformoutput',false);
isvalid1 = ~cellfun(@isempty, match_R);
match1_R4(~isvalid1) = {[]};
match2_R4 = cellfun(@(x,y) cat(2,x,y),match2_R,match2_t, 'uniformoutput',false);
match2_R4 = cellfun(@(x) cat(1,x,[0 0 0 1]),match2_R4, 'uniformoutput',false);
isvalid2 = ~cellfun(@isempty, match2_R);
match2_R4(~isvalid2) = {[]};

% Effective R and t back to WCS
G_R = cell(n_S,n_S);
G_R4_12 = cell(n_S,n_S);
G_t = cell(n_S,n_S);
G_Rrev = cell(1,n_S);
G_trev = cell(1,n_S);
G_RMSE = inf(n_S,n_S);
G_RMSE_sum = inf(n_S,n_S);
G_RMSE_geo = inf(n_S,n_S);
G_RMSE_mult = inf(n_S,n_S);
G_RMSE_max = inf(n_S,n_S);
G_RMSEh_sum = inf(n_S,n_S);
G_RMSEh_geo = inf(n_S,n_S);
G_RMSEh_mult = inf(n_S,n_S);
G_RMSEh_max = inf(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        path = G_path{i,j};
        path2 = fliplr(G_path{i,j});
        Rtemp = eye(3);
        ttemp = zeros(3,1);
        Rtemprev = eye(3);
        ttemprev = zeros(3,1);
        if numel(path) >0
            R4temp1 = eye(4);
            R4temp2 = eye(4);
            RMSEtemp = 0;
            RMSEtemp_sum = 0;
            RMSEtemp_geo = 0;
            RMSEtemp_mult = 1;
            RMSEtemp_max = 0;
            RMSEtemph_sum = 0;
            RMSEtemph_geo = 0;
            RMSEtemph_mult = 1;
            RMSEtemph_max = 0;
            invalid_R = false;
            invalid_R2 = false;
            for k = 1:numel(path)-1;
                if any(size(match_R{path(k+1),path(k)})~=size(Rtemp));
                    invalid_R=true;
                    break
                end
                if any(size(match2_R{path2(k+1),path2(k)})~=size(Rtemp));
                    invalid_R2=true;
                    foo = 1;
                   % break
                end
                R4temp1 = match1_R4{path(k+1),path(k)}*R4temp1;
                R4temp2 = match2_R4{path2(k),path2(k+1)}*R4temp2;
                Rtemp = match_R{path(k+1),path(k)}*Rtemp;
                ttemp = match_R{path(k+1),path(k)}*(ttemp)+ match_t{path(k+1),path(k)};
                Rtemprev = match2_R{path2(k),path2(k+1)}*Rtemprev;
                ttemprev = match2_R{path2(k),path2(k+1)}*(ttemprev)+ match2_t{path(k+1),path(k)};
                RMSEtemp = RMSEtemp + match12_RMSE(path(k+1),path(k));
                RMSEtemp_sum = RMSEtemp_sum + match12_RMSE(path(k+1),path(k));
                RMSEtemp_geo = RMSEtemp_geo + (match12_RMSE(path(k+1),path(k)))^2;
                RMSEtemp_mult = RMSEtemp_mult * (match12_RMSE(path(k+1),path(k)))*100;
                RMSEtemp_max = max(RMSEtemp_max, match12_RMSE(path(k+1),path(k)));
                RMSEtemph_sum = RMSEtemph_sum + match12h_RMSE(path(k+1),path(k));
                RMSEtemph_geo = RMSEtemph_geo + (match12h_RMSE(path(k+1),path(k)))^2;
                RMSEtemph_mult = RMSEtemph_mult * (match12h_RMSE(path(k+1),path(k)))*100;
                RMSEtemph_max = max(RMSEtemph_max, match12h_RMSE(path(k+1),path(k)));
            end
            if ~invalid_R;
                G_R4_12{i,j} = R4temp2*R4temp1;
                G_R{i,j} = Rtemp;
                G_t{i,j} = ttemp;
                %G_Rrev{i,j} = Rtemprev;
                %G_trev{i,j} = ttemprev;
                G_RMSE(i,j) = RMSEtemp;
                G_RMSE_sum(i,j) = RMSEtemp_sum;
                G_RMSE_geo(i,j) = sqrt(RMSEtemp_geo);
                G_RMSE_mult(i,j) = RMSEtemp_mult./100;
                G_RMSE_max(i,j) = RMSEtemp_max;
                G_RMSEh_sum(i,j) = RMSEtemph_sum;
                G_RMSEh_geo(i,j) = sqrt(RMSEtemph_geo);
                G_RMSEh_mult(i,j) = RMSEtemph_mult./100;
                G_RMSEh_max(i,j) = RMSEtemph_max;
            end
        end
    end
end

% Path RMSE 
isvalid = ~cellfun(@isempty,G_R4_12);
G_RMSE_path = zeros(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
       if isvalid(i,j);
     P_LCSt = G_R4_12{i,j}*[P_LCS{i}; ones(1, size(P_LCS{i},2))];
     P_LCSt = P_LCSt(1:3,:);
     G_RMSE_path(i,j) = sqrt(nanmean(sum((P_LCSt - P_LCS{i}).^2)));
       end
    end
end

%{
    % Declare match from i-i identity
    for i = 1:n_S;
        match_R{i,i} = eye(3);
        match_t{i,i} = zeros(3,1);
        match_i{i,i} = 1:numel(P_LCS{i});
        match_j{i,i} = 1:numel(P_LCS{i});
    end
%}

% Fig: Dijkstra for each reference node
%{
 for i = 1:n_S;
        filepath_fig = sprintf('%sDijstra_%03.0f-%02.0f.png',...
            path_save, info_site, i);
    if exist(filepath_fig, 'file');
        continue
    end
    RMSE_color = vec2cmap(G_RMSE(i,:),'jet',0,1);

        x_axis = [0 0 0; 1 0 0]';
        y_axis = [0 0 0; 0 1 0]';
        z_axis = [0 0 0; 0 0 1]';
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
       % err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        figure;
        scatter3([0 0],[0 0],[0,0],1, [0,1]);%[min(G_RMSE(i,:)) max(G_RMSE(i,:))])
        colormap('jet');
        c = colorbar;
        c.Label.String = 'RMSE Error [m]';
        title(sprintf('Dijstra Spanning Tree for Site %3.0f-%g', info_site,i));
        hold on;
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        axis equal
        for j= 1:n_S;
            hold on
            if isempty(G_R{i,j});
                continue
            end
            x_axist = G_R{i,j}*x_axis + repmat(G_t{i,j},[1,2]);
            y_axist = G_R{i,j}*y_axis + repmat(G_t{i,j},[1,2]);
            z_axist = G_R{i,j}*z_axis + repmat(G_t{i,j},[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            
            % Plot connections with RMSE
            plot3([0 G_t{i,j}(1)],[0 G_t{i,j}(2)],[0 G_t{i,j}(3)],...
                'color', RMSE_color(j,:),'linewidth',2);
            %h = surf(x_sph+G_t{i,j}(1), y_sph+G_t{i,j}(2), ...
            %    z_sph+G_t{i,j}(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
        view(0,90);
        set(gcf, 'color', 'white')
        export_fig(filepath_fig, '-m1');
        close(gcf);
 end
%}

G_rx =  180*cellfun(@decompose_rotation_rx,G_R)/pi;
G_ry =  180*cellfun(@decompose_rotation_ry,G_R)/pi;
G_rz =  180*cellfun(@decompose_rotation_rz,G_R)/pi;
G_tx =  cellfun(@decompose_translation_tx,G_t);
G_ty =  cellfun(@decompose_translation_ty,G_t);
G_tz =  cellfun(@decompose_translation_tz,G_t);

filepath_G_R = sprintf('%s%s',path_mat, 'G_R.mat');
if ~exist(filepath_G_R, 'file');
    save(filepath_G_R, 'G_R');
    filepath_G_t = sprintf('%s%s',path_mat, 'G_t.mat');
    save(filepath_G_t, 'G_t');
    filepath_G_tx = sprintf('%s%s',path_mat, 'G_tx.mat');
    save(filepath_G_tx, 'G_tx');
    filepath_G_ty = sprintf('%s%s',path_mat, 'G_ty.mat');
    save(filepath_G_ty, 'G_ty');
    filepath_G_tz = sprintf('%s%s',path_mat, 'G_tz.mat');
    save(filepath_G_tz, 'G_tz');
    filepath_G_rx = sprintf('%s%s',path_mat, 'G_rx.mat');
    save(filepath_G_rx, 'G_rx');
    filepath_G_ry = sprintf('%s%s',path_mat, 'G_ry.mat');
    save(filepath_G_ry, 'G_ry');
    filepath_G_rz = sprintf('%s%s',path_mat, 'G_rz.mat');
    save(filepath_G_rz, 'G_rz');
end
clear i j path Rtemp ttemp k Gsparse
% Outputs
% G_path               (j) best path to i
% G_R                  (j) effective Rotation to i
% G_t                  (j) effective translation to i
% G_rx
% G_tx
%}
%%
% Fig: Dijkstra for each reference node with piecewise edge
%{
 for i = 1:1;%n_S;
        filepath_fig = sprintf('%sDijstra_piecewise_%03.0f-%02.0f.png',...
            path_save, info_site, i);
    if exist(filepath_fig, 'file');
        continue
    end
        x_axis = [0 0 0; 1 0 0]';
        y_axis = [0 0 0; 0 1 0]';
        z_axis = [0 0 0; 0 0 1]';
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
       % err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        figure;
        scatter3([0 0],[0 0],[0,0],1, [min(G_RMSE(i,:)) max(G_RMSE(i,:))])
        colormap('jet');
        c = colorbar;
        c.Label.String = 'RMSE Error [m]';
        title(sprintf('Dijstra Spanning Tree for Site %3.0f-%g', info_site,i));
        hold on;
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        axis equal
        for j= 1:n_S;
            hold on
            if isempty(G_R{i,j});
                continue
            end
          
            x_axist = G_R{i,j}*x_axis + repmat(G_t{i,j},[1,2]);
            y_axist = G_R{i,j}*y_axis + repmat(G_t{i,j},[1,2]);
            z_axist = G_R{i,j}*z_axis + repmat(G_t{i,j},[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            
            % Plot connections with RMSE
         %   plot3([0 G_t{i,j}(1)],[0 G_t{i,j}(2)],[0 G_t{i,j}(3)],...
         %       'color', RMSE_color(j,:),'linewidth',2);
            %h = surf(x_sph+G_t{i,j}(1), y_sph+G_t{i,j}(2), ...
            %    z_sph+G_t{i,j}(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
        % Collate RMSE values for colormap
        RMSE_all = [];
        RMSE_temp = cell(n_S,1);
        for j = 1:n_S;
            RMSE_temp{j} = cell(numel(G_path{i,j})-1,1);
            for p = 1:numel(G_path{i,j})-1;
                RMSE_all = [RMSE_all, match12_RMSE(G_path{i,j}(p+1), G_path{i,j}(p))];
                RMSE_temp{j}{p} = match12_RMSE(G_path{i,j}(p+1), G_path{i,j}(p));
            end
        end
        if isempty(RMSE_all);
            continue
        end
        color_all = vec2cmap(RMSE_all,'jet');
        color_temp= RMSE_temp;
        for j = 1:n_S;
            for p = 1:numel(G_path{i,j})-1;
                color_temp{j}{p} = vec2cmap(RMSE_temp{j}{p},'jet', 0, 1);
            end
        end
        
        scatter3([0 0],[0 0],[0,0],1, [0,1]);%[min(RMSE_all) max(RMSE_all)])
        c = colorbar;
        c.Label.String = 'RMSE Error [m]';
        for j = 1:n_S;
            for p = 1:numel(G_path{i,j})-1;
                node1 = G_t{i,G_path{i,j}(p)};
                node2 = G_t{i,G_path{i,j}(p+1)};
                if isempty(node1) || isempty(node2);
                    continue
                end
                plot3([node1(1) node2(1)], [node1(2) node2(2)],...
                    [node1(3) node2(3)],'color', color_temp{j}{p}, 'linewidth',2)
            end
        end
        view(0,90);
        set(gcf, 'color', 'white')
       % export_fig(filepath_fig, '-m1');
        close(gcf);
 end
 
%}
%% Determine root node

% based on minimizing path distances
%node_length = sqrt(G_tx.^2 + G_ty.^2 + G_tz);
%node_meanlength = mean(node_length, 2);
%[~,ix_wcs] = min(node_meanlength);

% based on minimizing RMSE
isvalid = ~isinf(G_RMSE);
n_scans_connected = sum(isvalid,2);
G_RMSE(~isvalid)= nan;
isvalid_rows = (n_scans_connected==max(n_scans_connected));
G_RMSE(~isvalid_rows,:) = nan;
node_length = nansum(G_RMSE,2);
node_length(~isvalid_rows) = inf;
%node_length(sum(isvalid,2)==1) = inf;

[~,ix_wcs] = min(node_length);
%% Best Dijkstra from minimum total RMSE

% Minimal Spanning Tree Consensus Location
G_color = vec2cmap(G_RMSE(ix_wcs,:), 'jet');
filepath_fig = sprintf('%sMST_%03.0f.png',...
    path_mat, info_site);
if ~exist(filepath_fig,'file')&&options_showfig;
    figure;
    title(sprintf('Minimal Spanning Tree for Site %3.0f', info_site));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    if ~isempty(G_RMSE(ix_wcs,:));
        temp = G_RMSE(ix_wcs,:);
        temp = temp(~isinf(temp));
        h1 = scatter3([0 0],[0 0], [0 0], 10, [min(temp) max(temp)], 'filled');
    end
    colormap('jet')
    c = colorbar;
    c.Label.String = 'RMSE Error [m]';
    %set(h1, 'visible', 'off');
    for j = 1:n_S
        if isempty(G_R{ix_wcs,j});
            continue
        end
        x_axist = G_R{ix_wcs,j}*x_axis + repmat(G_t{ix_wcs,j},[1,2]);
        y_axist = G_R{ix_wcs,j}*y_axis + repmat(G_t{ix_wcs,j},[1,2]);
        z_axist = G_R{ix_wcs,j}*z_axis + repmat(G_t{ix_wcs,j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            G_color(j,:), 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            G_color(j,:), 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            G_color(j,:), 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
    end
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    close(gcf);
end
%% Save estimates for MST
G_R_MST = G_R(ix_wcs,:)';
G_t_MST = G_t(ix_wcs,:)';
G_RMSE_MST = G_RMSE(ix_wcs,:)';
G_RMSE_sum_MST = G_RMSE_sum(ix_wcs,:)';
G_RMSE_geo_MST = G_RMSE_geo(ix_wcs,:)';
G_RMSE_mult_MST = G_RMSE_mult(ix_wcs,:)';
G_RMSE_max_MST = G_RMSE_max(ix_wcs,:)';
G_RMSE_path_MST = G_RMSE_path(ix_wcs,:)';
G_RMSEh_sum_MST = G_RMSEh_sum(ix_wcs,:)';
G_RMSEh_geo_MST = G_RMSEh_geo(ix_wcs,:)';
G_RMSEh_mult_MST = G_RMSEh_mult(ix_wcs,:)';
G_RMSEh_max_MST = G_RMSEh_max(ix_wcs,:)';
filepath_G_R_MST = sprintf('%s%s',path_mat, 'G_R_MST.mat');
%if ~exist(filepath_G_R_MST, 'file');
save(filepath_G_R_MST, 'G_R_MST');
filepath_G_t_MST = sprintf('%s%s',path_mat, 'G_t_MST.mat');
save(filepath_G_t_MST, 'G_t_MST');
filepath_G_RMSE_MST = sprintf('%s%s',path_mat, 'G_RMSE_MST.mat');
save(filepath_G_RMSE_MST, 'G_RMSE_MST');
filepath_G_RMSE_sum_MST = sprintf('%s%s',path_mat, 'G_RMSE_sum_MST.mat');
save(filepath_G_RMSE_sum_MST, 'G_RMSE_sum_MST');
filepath_G_RMSE_geo_MST = sprintf('%s%s',path_mat, 'G_RMSE_geo_MST.mat');
save(filepath_G_RMSE_geo_MST, 'G_RMSE_geo_MST');
filepath_G_RMSE_mult_MST = sprintf('%s%s',path_mat, 'G_RMSE_mult_MST.mat');
save(filepath_G_RMSE_mult_MST, 'G_RMSE_mult_MST');
filepath_G_RMSE_max_MST = sprintf('%s%s',path_mat, 'G_RMSE_max_MST.mat');
save(filepath_G_RMSE_max_MST, 'G_RMSE_max_MST');
filepath_G_RMSE_path_MST = sprintf('%s%s',path_mat, 'G_RMSE_path_MST.mat');
save(filepath_G_RMSE_path_MST, 'G_RMSE_path_MST');
% half 
filepath_G_RMSEh_sum_MST = sprintf('%s%s',path_mat, 'G_RMSEh_sum_MST.mat');
save(filepath_G_RMSEh_sum_MST, 'G_RMSEh_sum_MST');
filepath_G_RMSEh_geo_MST = sprintf('%s%s',path_mat, 'G_RMSEh_geo_MST.mat');
save(filepath_G_RMSEh_geo_MST, 'G_RMSEh_geo_MST');
filepath_G_RMSEh_mult_MST = sprintf('%s%s',path_mat, 'G_RMSEh_mult_MST.mat');
save(filepath_G_RMSEh_mult_MST, 'G_RMSEh_mult_MST');
filepath_G_RMSEh_max_MST = sprintf('%s%s',path_mat, 'G_RMSEh_max_MST.mat');
save(filepath_G_RMSEh_max_MST, 'G_RMSEh_max_MST');
%end
%% Declare G_RMSE as G_RMSE_path for WMF and SVD 
G_RMSE = G_RMSE_sum;
%% Convert all Dijkstra's into WCS
G_R_WCS = cell(n_S,n_S);
G_t_WCS = cell(n_S,n_S);
w = ix_wcs;
% for i = 1:n_S;
%     for j = 1:n_S;
%         G_R_WCS{i,j} = G_R{w,i}*G_R{i,j};
%         G_t_WCS{i,j} = G_R{w,i}*G_t{i,j} + G_t{w,i};
%     end
% end


for i = 1:n_S;
    for j = 1:n_S;
        if isempty(G_R{i,w}) || isempty(G_R{i,j});
            continue
        end
        G_R_WCS{i,j} = G_R{i,w}'*G_R{i,j};
        G_t_WCS{i,j} = G_R{i,w}'*G_t{i,j} - G_R{i,w}'*G_t{i,w};
    end
end

G_RMSE_WCS = G_RMSE;
G_RMSE_WCS(:,w)= 0;
for i = 1:n_S;
    for j = 1:n_S;
        if j ~=w
            G_RMSE_WCS(i,j) = G_RMSE_WCS(i,j) + G_RMSE(i,w);
        end
    end
end

%All Dijkstra in WCS red green blue axes
%{
filepath_fig = sprintf('%sDijkstra_WCS_%03.0f.png',...
    path_mat, info_site);
if ~exist(filepath_fig,'file');
    figure;
    title(sprintf('Graph Locations for Site %3.0f in WCS', info_site));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    for i = 1:n_S
        hold on
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
        % err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        for j= 1:n_S;
            if isempty(G_R_WCS{i,j});
                continue
            end
            x_axist = G_R_WCS{i,j}*x_axis + repmat(G_t_WCS{i,j},[1,2]);
            y_axist = G_R_WCS{i,j}*y_axis + repmat(G_t_WCS{i,j},[1,2]);
            z_axist = G_R_WCS{i,j}*z_axis + repmat(G_t_WCS{i,j},[1,2]);
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            %h = surf(x_sph+G_t{i,j}(1), y_sph+G_t{i,j}(2), ...
            %    z_sph+G_t{i,j}(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
    end
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    close(gcf);
end
%}

%% Apply errors in sequence to estimate path-error

G_RMSE_WCS_vec = reshape(G_RMSE_WCS, 1,numel(G_RMSE_WCS));
color = vec2cmap(G_RMSE_WCS_vec,'jet');
G_color = reshape(color, n_S,n_S,3);
%
filepath_fig = sprintf('%sDijkstra_WCS_%03.0f.png',...
    path_mat, info_site);
if ~exist(filepath_fig, 'file')&&options_showfig;
    figure;
    title(sprintf('Graph Locations for Site %3.0f-%g', info_site,i));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    if ~isempty(G_RMSE_WCS_vec);
        temp = G_RMSE_WCS_vec;
        temp = temp(~isinf(temp));
        h1 = scatter3([0 0],[0 0], [0 0], 10, [min(temp) max(temp)], 'filled');
    end
    colormap('jet')
    c = colorbar;
    c.Label.String = 'RMSE Error [m]';
    for i = 1:n_S
        hold on
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
        % err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        for j= 1:n_S;
            if isempty(G_R_WCS{i,j});
                continue
            end
            x_axist = G_R_WCS{i,j}*x_axis + repmat(G_t_WCS{i,j},[1,2]);
            y_axist = G_R_WCS{i,j}*y_axis + repmat(G_t_WCS{i,j},[1,2]);
            z_axist = G_R_WCS{i,j}*z_axis + repmat(G_t_WCS{i,j},[1,2]);
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
                squeeze(G_color(i,j,:)), 'linewidth',2);
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
                squeeze(G_color(i,j,:)), 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
                squeeze(G_color(i,j,:)), 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            foo = 1;
            %h = surf(x_sph+G_t{i,j}(1), y_sph+G_t{i,j}(2), ...
            %    z_sph+G_t{i,j}(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
    end
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    close(gcf);
end
%}
%% Weighted consensus
G_rx_WCS =  180*cellfun(@decompose_rotation_rx,G_R_WCS)/pi;
G_ry_WCS =  180*cellfun(@decompose_rotation_ry,G_R_WCS)/pi;
G_rz_WCS =  180*cellfun(@decompose_rotation_rz,G_R_WCS)/pi;
G_tx_WCS =  cellfun(@decompose_translation_tx,G_t_WCS);
G_ty_WCS =  cellfun(@decompose_translation_ty,G_t_WCS);
G_tz_WCS =  cellfun(@decompose_translation_tz,G_t_WCS);

%{
G_rx_M = median(G_rx_WCS,1);
G_ry_M = median(G_ry_WCS,1);
G_rz_M = median(G_rz_WCS,1);
G_tx_M = median(G_tx_WCS,1);
G_ty_M = median(G_ty_WCS,1);
G_tz_M = median(G_tz_WCS,1);

G_rx_Mdiff = abs(G_rx_WCS - repmat(G_rx_M,[n_S,1]));
G_ry_Mdiff = abs(G_ry_WCS - repmat(G_ry_M,[n_S,1]));
G_rz_Mdiff = abs(G_rz_WCS - repmat(G_rz_M,[n_S,1]));
G_tx_Mdiff = abs(G_tx_WCS - repmat(G_tx_M,[n_S,1]));
G_ty_Mdiff = abs(G_ty_WCS - repmat(G_ty_M,[n_S,1]));
G_tz_Mdiff = abs(G_tz_WCS - repmat(G_tz_M,[n_S,1]));
%}

a = 1;
r = 5; 
%t = linspace(0,max(G_RMSE_WCS(~isinf(G_RMSE_WCS))),100);
wmf_fun = @(a,r,t) a*(exp(-r.*t));
%{
t = linspace(0,2,100);
%wmf_fun = @(t,xi)
figure; 
plot(t,wmf_fun(a,r,t)); 
hold on
plot(t,wmf_fun(a,0,t)); 
plot(t, 1./t);
ylim([0 1])
legend_str = {'exp -5x', 'exp{x}', '1/x'};
legend(legend_str);
%}

if aux.sensitivity_r == true;
    rarr = 0:10;
    SENS_ANAL_G_rx_WMF = nan(n_S,numel(rarr));
    SENS_ANAL_G_ry_WMF = nan(n_S,numel(rarr));
    SENS_ANAL_G_rz_WMF = nan(n_S,numel(rarr));
    SENS_ANAL_G_tx_WMF = nan(n_S,numel(rarr));
    SENS_ANAL_G_ty_WMF = nan(n_S,numel(rarr));
    SENS_ANAL_G_tz_WMF = nan(n_S,numel(rarr));
    SENS_ANAL_G_RMSE_WMF = nan(n_S,numel(rarr));
    for it = 1:numel(rarr)
        r = rarr(it);
            G_weight_WCS = wmf_fun(a,r,G_RMSE_WCS);
            G_rx_WMF = (nansum(G_rx_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_rx_WMF(:,it) = G_rx_WMF;
            G_ry_WMF = (nansum(G_ry_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_ry_WMF(:,it) = G_ry_WMF;
            G_rz_WMF = (nansum(G_rz_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_rz_WMF(:,it) = G_rz_WMF;
            G_tx_WMF = (nansum(G_tx_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_tx_WMF(:,it) = G_tx_WMF;
            G_ty_WMF = (nansum(G_ty_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_ty_WMF(:,it) = G_ty_WMF;
            G_tz_WMF = (nansum(G_tz_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_tz_WMF(:,it) = G_tz_WMF;
            G_RMSE_WMF = (nansum(G_RMSE_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
            SENS_ANAL_G_RMSE_WMF(:,it) = G_RMSE_WMF;
    end
    filepath_SENS_ANAL_G_rx_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_rx_WMF.mat');
    save(filepath_SENS_ANAL_G_rx_WMF, 'SENS_ANAL_G_rx_WMF');
    filepath_SENS_ANAL_G_ry_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_ry_WMF.mat');
    save(filepath_SENS_ANAL_G_ry_WMF, 'SENS_ANAL_G_ry_WMF');
    filepath_SENS_ANAL_G_rz_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_rz_WMF.mat');
    save(filepath_SENS_ANAL_G_rz_WMF, 'SENS_ANAL_G_rz_WMF');
    filepath_SENS_ANAL_G_tx_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_tx_WMF.mat');
    save(filepath_SENS_ANAL_G_tx_WMF, 'SENS_ANAL_G_tx_WMF');
    filepath_SENS_ANAL_G_ty_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_ty_WMF.mat');
    save(filepath_SENS_ANAL_G_ty_WMF, 'SENS_ANAL_G_ty_WMF');
    filepath_SENS_ANAL_G_tz_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_tz_WMF.mat');
    save(filepath_SENS_ANAL_G_tz_WMF, 'SENS_ANAL_G_tz_WMF');
    filepath_SENS_ANAL_G_RMSE_WMF = sprintf('%s%s',path_mat, 'SENS_ANAL_G_RMSE_WMF.mat');
    save(filepath_SENS_ANAL_G_RMSE_WMF, 'SENS_ANAL_G_RMSE_WMF');    
    return
end

if aux.sensitivity_r ==false;
G_weight_WCS = wmf_fun(a,r,G_RMSE_WCS);
G_rx_WMF = (nansum(G_rx_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
G_ry_WMF = (nansum(G_ry_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
G_rz_WMF = (nansum(G_rz_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
G_tx_WMF = (nansum(G_tx_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
G_ty_WMF = (nansum(G_ty_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
G_tz_WMF = (nansum(G_tz_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
G_RMSE_WMF = (nansum(G_RMSE_WCS.*G_weight_WCS,1)./nansum(G_weight_WCS,1))';
%else 
end

G_R_WMF = cell(n_S,1);
G_t_WMF = cell(n_S,1);
for j = 1:n_S;
    G_R_WMF{j} = compose_rotation(pi*G_rx_WMF(j)./180,...
        pi*G_ry_WMF(j)./180,pi*G_rz_WMF(j)./180);
    G_t_WMF{j} = [G_tx_WMF(j) G_ty_WMF(j) G_tz_WMF(j)]';
end

%{
filepath_fig = sprintf('%sDijkstra_WMF_%03.0f.png',...
    path_mat, info_site);
if ~exist(filepath_fig, 'file');
    figure;
    title(sprintf('WCS Dijkstra with WMF for %3.0f', info_site));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    for j = 1:n_S
        x_axist = G_R_WMF{j}*x_axis + repmat(G_t_WMF{j},[1,2]);
        y_axist = G_R_WMF{j}*y_axis + repmat(G_t_WMF{j},[1,2]);
        z_axist = G_R_WMF{j}*z_axis + repmat(G_t_WMF{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'k', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'k', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'k', 'linewidth',2)
    end
    for i = 1:n_S
        hold on
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
        % err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        for j= 1:n_S;
            if isempty(G_R_WCS{i,j});
                continue
            end
            x_axist = G_R_WCS{i,j}*x_axis + repmat(G_t_WCS{i,j},[1,2]);
            y_axist = G_R_WCS{i,j}*y_axis + repmat(G_t_WCS{i,j},[1,2]);
            z_axist = G_R_WCS{i,j}*z_axis + repmat(G_t_WCS{i,j},[1,2]);
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
                squeeze(G_color(i,j,:)), 'linewidth',2);
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
                squeeze(G_color(i,j,:)), 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
                squeeze(G_color(i,j,:)), 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            foo = 1;
            %h = surf(x_sph+G_t{i,j}(1), y_sph+G_t{i,j}(2), ...
            %    z_sph+G_t{i,j}(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
    end
    view(0,90);
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    close(gcf);
end
%}

% Weighted Majority Consensus Location
G_color_WMF = vec2cmap(G_RMSE_WMF, 'jet');
filepath_fig = sprintf('%sWMF_%03.0f.png',...
    path_mat, info_site);
if ~exist(filepath_fig,'file')&&options_showfig;
    figure;
    title(sprintf('Weighted Consensus Locations for Site %3.0f', info_site));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    if ~isempty(G_RMSE_WMF);
        temp = G_RMSE_WMF;
        temp = temp(~isinf(temp));
        h1 = scatter3([0 0],[0 0], [0 0], 10, [min(temp) max(temp)], 'filled');
    end
    colormap('jet')
    c = colorbar;
    c.Label.String = 'RMSE Error [m]';
    %set(h1, 'visible', 'off');
    for j = 1:n_S
        if isempty(G_R_WMF{j});
            continue
        end
        x_axist = G_R_WMF{j}*x_axis + repmat(G_t_WMF{j},[1,2]);
        y_axist = G_R_WMF{j}*y_axis + repmat(G_t_WMF{j},[1,2]);
        z_axist = G_R_WMF{j}*z_axis + repmat(G_t_WMF{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            G_color_WMF(j,:), 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            G_color_WMF(j,:), 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            G_color_WMF(j,:), 'linewidth',2)
    end
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    close(gcf);
end

%% Save variables
filepath_G_R_WMF = sprintf('%s%s',path_mat, 'G_R_WMF.mat');
%if ~exist(filepath_G_R_WMF, 'file');
save(filepath_G_R_WMF, 'G_R_WMF');
filepath_G_t_WMF = sprintf('%s%s',path_mat, 'G_t_WMF.mat');
save(filepath_G_t_WMF, 'G_t_WMF');
filepath_G_RMSE_WMF = sprintf('%s%s',path_mat, 'G_RMSE_WMF.mat');
save(filepath_G_RMSE_WMF, 'G_RMSE_WMF');
%{
        filepath_G_tx = sprintf('%s%s',path_mat, 'G_tx.mat');
        save(filepath_G_tx, 'G_tx');
        filepath_G_ty = sprintf('%s%s',path_mat, 'G_ty.mat');
        save(filepath_G_ty, 'G_ty');
        filepath_G_tz = sprintf('%s%s',path_mat, 'G_tz.mat');
        save(filepath_G_tz, 'G_tz');
        filepath_G_rx = sprintf('%s%s',path_mat, 'G_rx.mat');
        save(filepath_G_rx, 'G_rx');
        filepath_G_ry = sprintf('%s%s',path_mat, 'G_ry.mat');
        save(filepath_G_ry, 'G_ry');
        filepath_G_rz = sprintf('%s%s',path_mat, 'G_rz.mat');
        save(filepath_G_rz, 'G_rz');
%}
%end
%% Global optimization
% Have correspondences, just need reestimation
G_RMSE_SVD = zeros(n_S,1);
G_R_SVD = cell(n_S,1);
G_t_SVD = cell(n_S,1);

triW = nan(n_S,3);
w = ix_wcs;
for j = 1:n_S;
    if isempty(G_t{w,j});
        continue
    end
    triW(j,:) = G_t{w,j};
end

triI = cell(n_S,1);
triIt= cell(n_S,1);
isvalidw = ~isnan(triW(:,1));
for i = 1:n_S;
    triI{i} = nan(n_S,3);
    for j = 1:n_S;
        if isempty(G_t{i,j});
            continue
        end
        triI{i}(j,:) = G_t{i,j};
    end
    is_validi = ~isnan(triI{i}(:,1));
    is_valid = isvalidw & is_validi;
    if sum(is_valid)==1;
        continue
    end
    [G_R_SVD{i},G_t_SVD{i}] = rigid_transform_3D(triI{i}(is_valid,:),...
        triW(is_valid,:));
    triIt{i} = ((G_R_SVD{i}*triI{i}')+ repmat(G_t_SVD{i},1,n_S))'; %Note change
    temp = G_RMSE(:);
    temp(temp==inf)=nan;
    G_RMSE_SVD(i) = nanmean(temp);
    % Sanity check
    %{
    temp1 = triI{i}(is_valid,:);
    temp2 = triW(is_valid,:);
    temp3 = triIt{i};
    figure;
    scatter3(temp1(:,1), temp1(:,2),temp1(:,3),10,'r','filled');
    hold on
    scatter3(temp2(:,1), temp2(:,2),temp2(:,3),30,'b','filled');
    scatter3(temp3(:,1), temp3(:,2),temp3(:,3),10,'g','filled');
    legend_str = {'S_i','S_w','S_w>i',};
    legend(legend_str);
    for t = 1:size(temp1,1)
    plot3([temp1(t,1) temp2(t,1)], [temp1(t,2) temp2(t,2)],...
        [temp1(t,3) temp2(t,3)], '-k');
    end
    axis equal;
    view(0,90);
    foo = 1;
    %}
end

%

G_color_SVD = vec2cmap(G_RMSE_SVD, 'jet');
filepath_fig = sprintf('%sSVD_all_%03.0f.png',...
    path_mat, info_site);
if ~exist(filepath_fig,'file')&&options_showfig;
    figure;
    title(sprintf('SVD Locations for Site %3.0f', info_site));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    if ~isempty(G_RMSE_SVD);
        temp = G_RMSE_SVD;
        temp = temp(~isinf(temp));
        h1 = scatter3([0 0],[0 0], [0 0], 10, [min(temp) max(temp)], 'filled');
    end
    colormap('jet')
    c = colorbar;
    c.Label.String = 'RMSE Error [m]';
    %set(h1, 'visible', 'off');
    for j = 1:n_S
        if isempty(G_R_SVD{j});
            continue
        end
        x_axist = G_R_SVD{j}*x_axis + repmat(G_t_SVD{j},[1,2]);
        y_axist = G_R_SVD{j}*y_axis + repmat(G_t_SVD{j},[1,2]);
        z_axist = G_R_SVD{j}*z_axis + repmat(G_t_SVD{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            G_color_SVD(j,:), 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            G_color_SVD(j,:), 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            G_color_SVD(j,:), 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
    end
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    close(gcf);
end
%}

%% Save variables
filepath_G_R_SVD= sprintf('%s%s',path_mat, 'G_R_SVD.mat');
%if ~exist(filepath_G_R_SVD, 'file');
save(filepath_G_R_SVD, 'G_R_SVD');
filepath_G_t_SVD = sprintf('%s%s',path_mat, 'G_t_SVD.mat');
save(filepath_G_t_SVD, 'G_t_SVD');
filepath_G_RMSE_SVD = sprintf('%s%s',path_mat, 'G_RMSE_SVD.mat');
save(filepath_G_RMSE_SVD, 'G_RMSE_SVD');
%{
        filepath_G_tx = sprintf('%s%s',path_mat, 'G_tx.mat');
        save(filepath_G_tx, 'G_tx');
        filepath_G_ty = sprintf('%s%s',path_mat, 'G_ty.mat');
        save(filepath_G_ty, 'G_ty');
        filepath_G_tz = sprintf('%s%s',path_mat, 'G_tz.mat');
        save(filepath_G_tz, 'G_tz');
        filepath_G_rx = sprintf('%s%s',path_mat, 'G_rx.mat');
        save(filepath_G_rx, 'G_rx');
        filepath_G_ry = sprintf('%s%s',path_mat, 'G_ry.mat');
        save(filepath_G_ry, 'G_ry');
        filepath_G_rz = sprintf('%s%s',path_mat, 'G_rz.mat');
        save(filepath_G_rz, 'G_rz');
%}
%end
% Save w
filepath_w = sprintf('%s%s',path_mat, 'w.mat');
if ~exist(filepath_w, 'file');
    save(filepath_w, 'w');
end
%% Point Cloud Vector Field
%{
% Initialize volume
interval = 1;
vec_x = transpose(floor(options_axesval(1)):interval:ceil(options_axesval(2)));
vec_y = transpose(floor(options_axesval(3)):interval:ceil(options_axesval(4)));
vec_z = transpose(-5:interval:options_axesval(2));
[vec_X, vec_Y, vec_Z] = meshgrid(vec_x, vec_y, vec_z);
vec_xx = vec_X(:);
vec_yy = vec_Y(:);
vec_zz = vec_Z(:);
vec_xxyyzz = [vec_xx vec_yy vec_zz]';
n_vec = numel(vec_xx);

% Compute error vector
vec_uuvvww = cell(n_S,1);
for j = 2:n_S;
    vec_t =  G12_R{j}*vec_xxyyzz+ repmat(G12_T{j},1,n_vec);
    vec_uuvvww{j} = vec_t - vec_xxyyzz;
end
vec_mag = cellfun(@(x) sqrt(sum(x.^2,1)), vec_uuvvww, 'uniformOutput', false);

% Plot
%{
j = 2;
max_mag = max(vec_mag{j});
n_quivercolor = 25;
quiverstep = linspace(0,max_mag,n_quivercolor);
quivercolor = jet(n_quivercolor-1);
figure;
hold on
legend_str = cell(n_quivercolor-1,1);
for c = 1:n_quivercolor-1;
    is = vec_mag{j}>=quiverstep(c) & vec_mag{j}<=quiverstep(c+1);
    quiver3(vec_xx(is), vec_yy(is), vec_zz(is),...
        vec_uuvvww{j}(1,is)', vec_uuvvww{j}(2,is)', vec_uuvvww{j}(3,is)',...
        'color', quivercolor(c,:))
    legend_str{c} = sprintf('%2.2f < e < %2.2f', quiverstep(c), quiverstep(c+1));
end
xlabel('x relative [m]');
ylabel('y relative [m]');
zlabel('z relative [m]');
legend(legend_str,'location','bestoutside');
%}
clear interval vec_x vec_y vec_z vec_X vec_Y vec_Z vec_xx vec_yy vec_zz
clear vec_xxyyzz n_vec vec_uuvvww j vec_t
clear max_mag n_quivercolor quiverstep quivercolor legend_str is
% Outputs
% vector field figure
%}
%% Per-point Error
%{
filepath_plyte = cell(n_S,1);
for j = 1:n_S;
    fprintf('\nWriting ply %g of %g\n',j,n_S);
    filepath_plyte{j} = sprintf('%ste.ply',filepath_ply{j}(1:end-4));
    [vertex, ~] = read_ply(filepath_ply{j});
    data_x = vertex(:,1);
    data_y = vertex(:,2);
    data_z = vertex(:,3);
    n_data = numel(data_x);
    data2_xyz = (G_R{j}*[data_x data_y data_z]') + ...
        repmat(G_t{j},1,n_data);
    data12_xyz = (G12_R{j}*[data_x data_y data_z]') + ...
        repmat(G12_T{j},1,n_data);
    data_e = sqrt(sum((data12_xyz - [data_x data_y data_z]').^2,1));
    color = plyintensity2color(data_e, 'jet');
    if max(color(:))<.001;
        color = repmat([ 0 0 127.5], n_data,1);
    end
    write2ply(filepath_plyte{j},data2_xyz',color);
end
clear j vertex dataJ_x dataJ_y dataJ_z n_data
clear data2_xyz data12_xyz data_e color
% Outputs
% filepath_plyte        filepath for ply transformed with color = error
% PLY with color according to error written to disk
%}
%% Per-point Vector Field
%{
j = 2;
[vertex, ~] = read_ply(filepath_ply{j});
data_x = vertex(:,1);
data_y = vertex(:,2);
data_z = vertex(:,3);
n_data = numel(data_x);
data12_xyz = (G12_R{j}*[data_x data_y data_z]') + ...
    repmat(G12_T{j},1,n_data);
data_uvw = data12_xyz - [data_x data_y data_z]';
data_e = sqrt(sum((data_uvw).^2,1));

% Plot
figure;
hold on
max_mag = max(data_e);
n_quivercolor = 25;
quiverstep = linspace(0,max_mag,n_quivercolor);
quivercolor = jet(n_quivercolor-1);
legend_str = cell(n_quivercolor-1,1);
for c = 1:n_quivercolor-1;
    is = data_e>=quiverstep(c) & data_e<=quiverstep(c+1);
    quiver3(data_x(is), data_y(is), data_z(is),...
        data_uvw(1,is)', data_uvw(2,is)', data_uvw(3,is)',...
        'color', quivercolor(c,:))
    legend_str{c} = sprintf('%2.2f < e < %2.2f', quiverstep(c), quiverstep(c+1));
end
xlabel('x relative [m]');
ylabel('y relative [m]');
zlabel('z relative [m]');
legend(legend_str,'location','bestoutside');
%
clear j vertex data_x data_y data_z n_data data12_xyz
clear data_uv data_e n_quivercolor quiverstep quivercolor legend_str
clear c is
% Outputs
% Outputs
% per point vector field figure
%}

end



