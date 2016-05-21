%% Validation: tree locations 
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
%% Compare tree locations 

G_tree_hat = cell(n_g,1);
G_tree_true = cell(n_g,1);
GR4_eye = cellfun(@double,repmat({[eye(3) zeros(3,1); zeros(1,3) 1]},n_g,1),'Un',0);
for s = 1:n_S;
    is_valid = (G_start==s);
    P4_temp = [P_LCS{s}; ones(1,size(P_LCS{s},2))];
    GR4_temp = G_R4(is_valid);
    Geye_temp = GR4_eye(is_valid);
    G_tree_hat(is_valid) = cellfun(@(x) x*P4_temp, GR4_temp, 'uniformoutput',false);
    G_tree_true(is_valid) =  cellfun(@(x) x*P4_temp, Geye_temp, 'uniformoutput',false);
end

n_inG = cellfun(@(x) size(x,2), G_tree_hat);
n_all = sum(n_inG);

all_tree_hat = zeros(n_all,3);
all_tree_true = zeros(n_all,3); 
all_npath = zeros(n_all,1);
all_dist = zeros(n_all,1);
all_minmatch = zeros(n_all,1);
start_ix = [0; cumsum(n_inG)]+1;
end_ix = circshift(start_ix, [-1,0])-1;
start_ix = start_ix(1:end-1);
end_ix = end_ix(1:end-1);

for g = 1:n_g;
    all_tree_hat(start_ix(g):end_ix(g),:) = G_tree_hat{g}(1:3,:)';
    all_tree_true(start_ix(g):end_ix(g),:) = G_tree_true{g}(1:3,:)';
    all_npath(start_ix(g):end_ix(g)) = repmat(G_npath(g),[n_inG(g),1]);
    all_dist(start_ix(g):end_ix(g)) = repmat(G_dist(g),[n_inG(g),1]);
    all_minmatch(start_ix(g):end_ix(g)) = repmat(G_min_nmatch(g),[n_inG(g),1]);
end

all_tree_exy = sqrt(sum((all_tree_hat(:,1:2) - all_tree_true(:,1:2)).^2,2)); 
all_tree_ex = all_tree_hat(:,1) - all_tree_true(:,1); 
all_tree_ey = all_tree_hat(:,2) - all_tree_true(:,2); 
all_tree_ez = all_tree_hat(:,3) - all_tree_true(:,3); 

 