function [ match_R, match_t ] = kelbe_analysis_rmse_fun( aux )
%UNTITLED Summary of this function goes here
%% KELBE REGISTRATION TRANSFORM
% Adds noise from normal distribution to stem locations
% Clear without breakpoints
tree = 1;
n_trial = aux.n_trial;
info_valid_plot = aux.info_valid_plot;
info_site = aux.info_site;
info_analysis = 'pose';

%% Load points
% Input to function are stem maps derived from lidar.m
% fprintf('\nLoad points\n');

% Filename Lookup for Harvard Data
info_exp = 'Harvard';
info_suffix = 'reg';
info_slash = '\';
path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
D = dir(path_site);
ctr = 1;
% Only load valid plots
n_S = numel(info_valid_plot);
n_tree = zeros(n_S,1);
P_LCS = cell(n_S,1);
P_rad = cell(n_S,1);
P_plot = zeros(n_S,1);
path_ply = cell(n_S,1);
filepath_ply = cell(n_S,1);

% First set
info_plot = info_valid_plot{1};
path_mat = sprintf('%s%s%smat%s',path_site,info_plot,info_slash,info_slash);
filepath_anal_rmse_step = sprintf('%sanal_rmse_step.mat', path_mat);
if exist(filepath_anal_rmse_step, 'file');
    return
end
filepath_tree = sprintf('%stree.mat',path_mat);
if ~exist(filepath_tree,'file')
    return
end
load(filepath_tree);
 if isempty(tree(1).loc);
     return
 end
n_tree(1) = numel(tree);
plot_dbl = str2double(info_plot);
P_LCS{1} = nan(3,n_tree(1));
P_rad{1} = nan(n_tree(1),1);
P_plot(1) = plot_dbl;
for t = 1:n_tree(1);
    P_LCS{1}(:,t) = tree(t).loc(:,1);
    P_rad{1}(t) = tree(t).r(1);
end
path_ply{1} = sprintf('%s%s%sply%s',path_site,info_plot,info_slash,info_slash);
filepath_ply{1} = sprintf('%spoints_full_%03.0f-%02.0f.ply', path_ply{ctr}, info_site, plot_dbl);

% Second set

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
P_plot(2) = P_plot(1);

%aux.P_color = P_color;
aux.P_plot = P_plot;
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
aux.info_analysis = info_analysis;

%% Add noise
% Set up trial
s_int = 0.025;
s_max = .5;
step_std = 0:s_int:s_max;
n_std = numel(step_std);
match_RMSE = nan(n_std, n_trial);
match_RMSEin = nan(n_std, n_trial);
match_noise = nan(n_std, n_trial);

% Run each trial
for s = 1:n_std;
 %   fprintf('\n Working on sigma %g of %g\n',s,n_step);
    for t = 1:n_trial;
        % Add noise
        P_noise = step_std(s).*randn(3,n_tree(1));
        P_LCS{2} = P_LCS{1} + P_noise;
        match_noise(s,t) = step_std(s);
        match_RMSEin(s,t) = sqrt(nanmean(sum((P_LCS{2} - P_LCS{1}).^2)));
        aux.P_LCS = P_LCS;
        % Perform registration
        [match_R, match_t] = kelbe_registration_combine_dis3fun_analyses( aux );
        match_R = match_R{1,2};
        match_t = match_t{1,2};
        if ~isempty(match_R);
            P_LCSt =  (match_R*P_LCS{2})+ repmat(match_t,1,n_tree(2));
            P_dist = sqrt(sum((P_LCSt - P_LCS{1}).^2));            
            match_RMSE(s,t) = mean(P_dist);
        end
        %}
    end
end

% Find relationship between sigma and rmse 
%{
figure;
plot(mean(match_noise,2), mean(match_RMSEin,2))
p =  polyfit(mean(match_noise,2), mean(match_RMSEin,2),1);
yfit = polyval(p,0.2);
%}

match_RMSEout = match_RMSE;
filepath_anal_RMSEin_RMSEin = sprintf('%sanal_noise_RMSEin.mat', path_mat);
filepath_anal_RMSEin_RMSEout = sprintf('%sanal_noise_RMSEout.mat', path_mat);
save(filepath_anal_RMSEin_RMSEin, 'match_RMSEin');
save(filepath_anal_RMSEin_RMSEout, 'match_RMSEout');

% Show results
%{
mean_RMSE = mean(match_RMSE,2);
figure;
hold on;
scatter(rz_step, mean_RMSE,'filled');
xlabel_str = 'Roll, Pitch, Yaw mean, Sigma = interval/4';
ylabel_str = 'RMSE of registration';
xlabel(xlabel_str);
ylabel(ylabel_str);
axis equal;
%}

% Show coordinate transformations used

%% Plot pose positions
%{
x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

figure;
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);

for s = 1:n_step;
    for t = 1:n_trial;
        R = compose_rotation(rx(s,t),ry(s,t),rz(s,t));
    T = [tx(s,t) ty(s,t) tz(s,t)]';
    x_axis_R = (R*x_axis')' + repmat(T', [2,1]);
    y_axis_R = (R*y_axis')' + repmat(T', [2,1]);
    z_axis_R = (R*z_axis')' + repmat(T', [2,1]);
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
    end
end
title('TRANSLATION,ROLL/PITCH/YAW');

clear s_std s_max n_std n_trial match_RMSE P_noise match_R match_t P_LCSt
clear P_dist
%}
end



