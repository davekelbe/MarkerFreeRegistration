%% KELBE REGISTRATION WIGGLE
% Adds noise from normal distribution to stem locations
% Clear without breakpoints
tmp = dbstatus;
save('tmp.mat','tmp');
clear all; close all; clc;
load('tmp.mat');
dbstop(tmp);
clear tmp;
delete('tmp.mat');
% Default settings
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = false;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
options_loadvar = false;
%if ~options_loadvar;
% Initialize program-level options and paths
% Outputs
% options               options
% paths                 paths
clear D ctr d filepath_tree plot t tree
if ~options_loadvar;
    %% Parameters
    t_rad = 0.2; %0.06^2; % Trees can be matched if their radius is within threshold
    t_coll = 0.1; % Collinnearity threshold: don't consider triangles that are collinear
    t_eig_error = 1e1;
    t_RANSAC_rad = 0.2;
    t_RANSAC_xyz = 0.4^2;
    t_RANSAC_nsearch = 16;
    
    % Initialize program-level options and paths
    % Outputs
    % t_*               thresholds
    %% Load points
    % Input to function are stem maps derived from lidar.m
    % fprintf('\nLoad points\n');
    
    % Filename Lookup for Harvard Data
    info_exp = 'Harvard';
    info_suffix = '03-01';
    info_slash = '\';
    info_site = 31;
    path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
    D = dir(path_site);
    ctr = 1;
    % Only load valid plots
    info_valid_plot = {'08', '08n'};
    n_S = numel(info_valid_plot);
    n_tree = zeros(n_S,1);
    P_LCS = cell(n_S,1);
    P_rad = cell(n_S,1);
    P_plot = zeros(n_S,1);
    path_ply = cell(n_S,1);
    filepath_ply = cell(n_S,1);
    
    for d = 1:numel(D);
        if strcmp(D(d).name(1), '.')
            continue
        end
        info_plot = D(d).name;
        if ~any(strcmp(info_plot, info_valid_plot))
            continue
        end
        path_mat = sprintf('%s%s%smat%s',path_site,info_plot,info_slash,info_slash);
        filepath_tree = sprintf('%stree.mat',path_mat);
        if ~exist(filepath_tree,'file')
            continue
        end
        load(filepath_tree);
        n_tree(ctr) = numel(tree);
        plot = str2double(info_plot);
        P_LCS{ctr} = nan(3,n_tree(ctr));
        P_rad{ctr} = nan(n_tree(ctr),1);
        P_plot(ctr) = plot;
        for t = 1:n_tree(ctr);
            P_LCS{ctr}(:,t) = tree(t).loc(:,1);
            P_rad{ctr}(t) = tree(t).r(1);
        end
        path_ply{ctr} = sprintf('%s%s%sply%s',path_site,info_plot,info_slash,info_slash);
        filepath_ply{ctr} = sprintf('%spoints_full_%03.0f-%02.0f.ply', path_ply{ctr}, info_site, plot);
        ctr = ctr + 1;
    end
    
    %[P_LCS, P_rad, P_plot, n_tree, n_S ] = generate_test_data;
    
    i_xmin = min(P_LCS{1}(1,:));
    i_xmax = max(P_LCS{1}(1,:));
    i_ymin = min(P_LCS{1}(2,:));
    i_ymax = max(P_LCS{1}(2,:));
    i_zmin = min(P_LCS{1}(3,:));
    i_zmax = max(P_LCS{1}(3,:));
    options_axesval = [i_xmin i_xmax i_ymin i_ymax i_zmin i_zmax];
    
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
    %%
    P_rad{2} = P_rad{1};
    n_tree(2) = n_tree(1);
    
    %aux.P_color = P_color;
    %aux.P_plot = P_plot;
    aux.P_rad = P_rad;
    %aux.filepath_ply = filepath_ply;
    %aux.info_exp = info_exp;
    %aux.info_plot = info_plot;
    %aux.info_site = info_site;
    %aux.info_slash = info_slash;
    %aux.info_suffix = info_suffix;
    %aux.info_valid_plot = info_valid_plot;
    aux.n_S = n_S;
    aux.n_tree = n_tree;
    %aux.options_axesval = options_axesval;
    %aux.options_imagepoints = options_imagepoints;
    %aux.options_initialmatch = options_initialmatch;
    %aux.options_loadmatch = options_loadmatch;
    %aux.options_loadvar = options_loadvar;
    %aux.options_unique = options_unique;
    %aux.options_verbose = options_verbose;
    %aux.p_std = p_std;
    %aux.path_mat = path_mat;
    %aux.path_ply = path_ply;
    %aux.path_save = path_save;
    %aux.path_site = path_site;
    %aux.path_tikz = path_tikz;
    aux.t_RANSAC_nsearch = t_RANSAC_nsearch;
    %aux.t_RANSAC_rad = t_RANSAC_rad;
    aux.t_RANSAC_xyz = t_RANSAC_xyz;
    aux.t_coll = t_coll;
    aux.t_eig_error = t_eig_error;
    aux.t_rad = t_rad;
    
    % Set up trial
    s_int = 5;
    s_min = 3;
    s_max = n_tree(1);
    step_sub = s_min:s_int:s_max;
    n_std = numel(step_sub);
    n_trial = 5;
    match_RMSE = zeros(n_std, n_trial);
    % Run each trial
    for s = 1:n_std;
        fprintf('\n Working on subset %g of %g\n',s,n_std);
        for t = 1:n_trial;
            % Add noise
            n_tree(2) = step_sub(s);
            ix = randperm(n_tree(1),step_sub(s));
            P_LCS{2} = P_LCS{1}(:,ix);
            P_rad{2} = P_rad{1}(ix);
            aux.P_LCS = P_LCS;
            aux.P_rad = P_rad;
            aux.n_tree = n_tree;
            % Perform registration
            [match_R, match_t] = kelbe_registration_function( aux );
            if ~isempty(match_R);
            P_LCSt =  (match_R*P_LCS{2})+ repmat(match_t,1,n_tree(2));
            P_dist = sqrt(sum((P_LCSt - P_LCS{1}(:,ix)).^2));
            match_RMSE(s,t) = mean(P_dist);
            end
        end
    end
    
    %save('kelbe_registration_wiggle_214.mat');
else
    load('kelbe_registration_subset_214.mat');
end

mean_RMSE = mean(match_RMSE,2);

figure;
hold on;
scatter(step_sub, mean_RMSE,'filled');
xlabel_str = 'Subset of stem locations';
ylabel_str = 'RMSE of registration';
xlabel(xlabel_str);
ylabel(ylabel_str);
axis equal;
axis([0 1.2*s_max 0 max(mean_RMSE)+1]);
clear s_std s_max n_std n_trial match_RMSE P_noise match_R match_t P_LCSt
clear P_dist


