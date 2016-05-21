function [  ] = kelbe_registration_notx( aux )
    % Only load valid plots
    %info_valid_plot = {'13','15'};
%      info_valid_plot = cell(1,25);
%      for s = 1:5;
%          info_valid_plot{s} = sprintf('%02.0f', s);
%      end
%      info_site = 31;
    info_valid_plot = aux.info_valid_plot;
    info_site = aux.info_site;
    path_mat = aux.path_mat;
    aux.P_LCS = aux.P_LCSn;
    aux.P_rad = aux.P_radn;
    
%% KELBE REGISTRATION
% Clear without breakpoints
%tmp = dbstatus;
%save('tmp.mat','tmp');
%clear all; close all; clc;
%load('tmp.mat');
%dbstop(tmp);
%clear tmp;
%delete('tmp.mat');
% Default settings
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = false;
path_save = 'D:\Users\djk2312\Documents\Figures\graphanalysis\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_matvar = 'D:\Users\djk2312\Documents\matvar\';
options_loadvar =false;
tree = 0;
if ~options_loadvar;
    % Initialize program-level options and paths
    % Outputs
    % options               options
    % paths                 paths
    clear D ctr d filepath_tree plot t tree
    %% Parameters
    t_rad = 0.2; %0.06^2; % Trees can be matched if their radius is within threshold
    t_coll = 0.1; % Collinnearity threshold: don't consider triangles that are collinear
    t_eig_error = 1e1;
    t_RANSAC_xyz = 0.4^2;
    t_RANSAC_nsearch = 1000;
    %t_flagr = 20; % deg
    %t_flagt = 1; % m
    
    
    % Initialize program-level options and paths
    % Outputs
    % t_*               thresholds
    %% Load points

    % Input to function are stem maps derived from lidar.m
    %fprintf('\nLoad points\n');
    
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
    
    i_xmin = min(cellfun(@(x) min(x(1,:)),P_LCS));
    i_xmax = max(cellfun(@(x) max(x(1,:)),P_LCS));
    i_ymin = min(cellfun(@(x) min(x(2,:)),P_LCS));
    i_ymax = max(cellfun(@(x) max(x(2,:)),P_LCS));
    i_zmin = min(cellfun(@(x) min(x(3,:)),P_LCS));
    i_zmax = max(cellfun(@(x) max(x(3,:)),P_LCS));
    options_axesval = [i_xmin i_xmax i_ymin i_ymax i_zmin i_zmax];
    
    % Colormap for sensors
    P_color = jet(n_S);
    
    % Individual camera views
    if false;%options_verbose && options_imagepoints;
        for s = 1:5:n_S;
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
            axis equal
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
n_T = aux.n_T;
P_ix = aux.P_ix;
PP_LCS = cell(n_S,1);
for i = 1:n_S;
    PP_LCS{i} = nan(n_T,3);
    PP_LCS{i}(P_ix{i},:) = P_LCS{i}';
end

match12_RMSE = nan(n_S,n_S);
match_ij = false(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
      %  if ~isnan(match12_er(i,j));
            match12_RMSE(i,j) = sqrt(nanmean(sum((PP_LCS{j} - PP_LCS{i}).^2,2)));
       % end
       if i==j;
           match_ij(i,j) = true;
       end 
    end
end
match_R  = cell(n_S,n_S);
match_t = cell(n_S,n_S);
match12_R  = cell(n_S,n_S);
match12_t = cell(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        match_R{i,j} = eye(3);
        match_t{i,j} = zeros(3,1);
        match12_R{i,j} = eye(3);
        match12_t{i,j} = zeros(3,1);
    end
end
        
%% Save variables 

    %filepath_nit = sprintf('%s%s',path_mat, 'match_nit.mat');
    %save(filepath_nit, 'match_nit');
    filepath_match12_RMSE = sprintf('%s%s',path_mat, 'match12_RMSE.mat');
    save(filepath_match12_RMSE, 'match12_RMSE');
    filepath_match_ij = sprintf('%s%s',path_mat, 'match_ij.mat');
    save(filepath_match_ij, 'match_ij');
    filepath_match12_R = sprintf('%s%s',path_mat, 'match12_R.mat');
    save(filepath_match12_R, 'match12_R');    
    filepath_match12_t = sprintf('%s%s',path_mat, 'match12_t.mat');
    save(filepath_match12_t, 'match12_t');
    filepath_match_R = sprintf('%s%s',path_mat, 'match_R.mat');
    save(filepath_match_R, 'match_R');    
    filepath_match_t = sprintf('%s%s',path_mat, 'match_t.mat');
    save(filepath_match_t, 'match_t');
    %filepath_match1_R = sprintf('%s%s',path_mat, 'match1_R.mat');
    %save(filepath_match1_R, 'match1_R');    
    %filepath_match1_t = sprintf('%s%s',path_mat, 'match1_t.mat');
    %save(filepath_match1_t, 'match1_t');
    %filepath_match2_R = sprintf('%s%s',path_mat, 'match2_R.mat');
    %save(filepath_match2_R, 'match2_R');    
    %filepath_match2_t = sprintf('%s%s',path_mat, 'match2_t.mat');
    %save(filepath_match2_t, 'match2_t'); 


end



