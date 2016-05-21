%% Validation: transformation parameters 
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
options_loadmatch = true;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_temp = 'D:\Users\djk2312\Documents\MATLAB\temp\';
load('tempvalidation.mat');
% Initialize program-level options and paths
% Outputs
% G_npath               number of nodes in each path
% G_rx                  array of x angles [degrees]
% G_ry                  array of y angles [degrees]
% G_rz                  array of z angles [degrees]
% G_tx                  array of x translations [m]
% G_ty                  array of y translations [m]
% G_tz                  array of z translations [m]
% info                  program level information
% ************************
% options               options
% paths                 paths
%% Rxyz vs # nodes 
options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_rxyz_vs_node(G_npath,G_rx, G_ry, G_rz,info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );


%% Txyz vs. # nodes 
options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_txyz_vs_node(G_npath,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% Rxyz vs path distance 

options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_rxyz_vs_node(G_npath,G_rx, G_ry, G_rz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% Txyz vs path distance 

options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_txyz_vs_node(G_npath,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% Rxyz vs path distance 

options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_rxyz_vs_node(G_npath,G_rx, G_ry, G_rz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% tloc vs paths
options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_tloc_vs_node(all_npath,all_tree_exy, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% tloc vs dist
options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_tloc_vs_dist(all_dist,all_tree_exy, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

