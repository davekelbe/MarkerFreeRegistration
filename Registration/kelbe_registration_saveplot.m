function [  ] = kelbe_registration_saveplot( aux )
    % Only load valid plots
    %info_valid_plot = {'13','15'};
%      info_valid_plot = cell(1,25);
%      for s = 1:5;
%          info_valid_plot{s} = sprintf('%02.0f', s);
%      end
%      info_site = 31;
    info_valid_plot = aux.info_valid_plot;
    info_site = aux.info_site;
    
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
path_save = 'D:\Users\djk2312\Documents\Figures\registration\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_matvar = 'D:\Users\djk2312\Documents\matvar\';
options_loadvar =false;
tree = 0;
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
    t_flagr = 20; % deg
    t_flagt = 1; % m
    
    
    % Initialize program-level options and paths
    % Outputs
    % t_*               thresholds
    %% Load points
    
    % Input to function are stem maps derived from lidar.m
    %fprintf('\nLoad points\n');
    
    % Filename Lookup for Harvard Data
    info_exp = 'Harvard';
    info_suffix = 'reg';
    info_slash = '\';
    path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
    D = dir(path_site);
    ctr = 1;

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
        if numel(tree)==1;
            if isempty(tree.r)
                continue
            end
        end
        if numel(tree)<=6;
                continue
        end
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
    isvalid = (P_plot~=0);
    P_LCS = P_LCS(isvalid);
    P_rad = P_rad(isvalid);
    P_plot = P_plot(isvalid);
    n_plot = numel(P_plot);
    match_plotI = repmat(P_plot, [1, n_plot]);
    match_plotJ = repmat(P_plot', [n_plot, 1]);
    filepath_match_plotI = sprintf('%s%s',path_mat, 'match_plotI.mat');
    save(filepath_match_plotI, 'match_plotI');
    filepath_match_plotJ = sprintf('%s%s',path_mat, 'match_plotJ.mat');
    save(filepath_match_plotJ, 'match_plotJ');   
end
    