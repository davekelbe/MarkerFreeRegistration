function [  ] = register_wrapper( info_site, info_plot1,info_plot2 )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%% Initialization
close all
%profile on

% Print current status
fprintf('\n**************\n');
fprintf('info_site  = %3.0f \n', info_site);
fprintf('info_plot1 = %3.0f \n', info_plot1);
fprintf('info_plot2 = %3.0f \n', info_plot2);


% Set Internal Matlab Parameters
set(0,'DefaultFigureRenderer','OpenGL')
set(0,'defaultfigureposition',[966  322  557 472]);
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:legend:IgnoringExtraEntries');
set(0,'DefaultFigureColor', 'White');
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontName', 'Calibri')
set(0,'DefaultAxesFontSize', 14)

% Set Run-type options
%options_cyl_alpha = 0.5;
%options_ptsize = 5;
options_verbose = true;
options_verbose_progress = true;
options_verbose_output = true;
options_show_fig = false;
%options_save_eps = false;
%options_save_fig = false;
%options_save_png = false;
options_write_ply = true;
options_rotation = false;

% Info
info_exp = 'Harvard';
info_suffix = '03-01';
info_slash = '\';

% Make directories
path_common = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, 'Common',info_slash);
path_top1 = sprintf('%s%s%s%s%s%03.0f%s%02.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash,info_plot1,info_slash);
path_top2 = sprintf('%s%s%s%s%s%03.0f%s%02.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash,info_plot2,info_slash);

path_mat = sprintf('%s%s%s',path_top2,'mat',info_slash);
path_fig = sprintf('%s%s%s',path_top2,'fig',info_slash);
path_eps = sprintf('%s%s%s',path_top2,'eps',info_slash);
path_ply1 = sprintf('%s%s%s',path_top1,'ply',info_slash);
path_ply2 = sprintf('%s%s%s',path_top2,'ply',info_slash);
path_local = 'Z:\Desktop\Local\';
path_png = sprintf('%s%s%s',path_top2,'png',info_slash);
path_latex = 'Z:\Desktop\Paper\results\';

% Radius difference
t_rad = 0.1;
% Eigenvalue difference
%t_eig = inf - just sorting now
% Distance between tree centers for match in RANSAC
t_xyz = 0.4.^2;
% Number of RANSAC iterations
n_search = 8;


% Outputs:
%       options               - user-defined run-type options
%       info                  - experiment-specific info

%% Load data

path_mat1 = sprintf('%s%s%s%s',path_top1,'mat',info_slash);
path_mat2 = sprintf('%s%s%s%s',path_top2,'mat',info_slash);

% Filepaths for saved matlab variables
filepath_data_ieq1  = sprintf('%s%s',path_mat1,'data_ieq.mat');
filepath_data_x1  = sprintf('%s%s',path_mat1,'data_x.mat');
filepath_data_y1  = sprintf('%s%s',path_mat1,'data_y.mat');
filepath_data_z1  = sprintf('%s%s',path_mat1,'data_z.mat');
% If matlab variables exist, just load them
if options_verbose_progress;
    fprintf('\n %% Loading data \n');
end
load(filepath_data_ieq1);
load(filepath_data_x1);
load(filepath_data_y1);
load(filepath_data_z1);
data_ieq1 = data_ieq;
data_x1 = data_x;
data_y1 = data_y;
data_z1 = data_z;



% Filepaths for saved matlab variables
filepath_data_ieq2  = sprintf('%s%s',path_mat2,'data_ieq.mat');
filepath_data_x2  = sprintf('%s%s',path_mat2,'data_x.mat');
filepath_data_y2  = sprintf('%s%s',path_mat2,'data_y.mat');
filepath_data_z2  = sprintf('%s%s',path_mat2,'data_z.mat');
% If matlab variables exist, just load them
if options_verbose_progress;
    fprintf('\n %% Loading data \n');
end
load(filepath_data_ieq2);
load(filepath_data_x2);
load(filepath_data_y2);
load(filepath_data_z2);
data_ieq2 = data_ieq;
data_x2 = data_x;
data_y2 = data_y;
data_z2 = data_z;

%% Run registration

filepath_R = sprintf('%sR_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);
filepath_t = sprintf('%st_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);
filepath_m1_best = sprintf('%sm1_best_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);
filepath_m2_best = sprintf('%sm2_best_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);

if exist(filepath_R,'file');
    load(filepath_R);
    load(filepath_t);
else
    [R, t,m1_best, m2_best] = registration1(info_site, info_plot1, info_plot2);
    save(filepath_R,'R');
    save(filepath_t,'t');
    save(filepath_m1_best,'m1_best');
    save(filepath_m2_best,'m2_best');   
end

%% Write ply file

xyz2t = (R*[data_x2 data_y2 data_z2]') + repmat(t,1,numel(data_x2));

if options_write_ply;
    % Just 1st returns
    filepath_ply2t  = sprintf('%s%s_%03.0f-%02.0f-%02.0f.ply',path_ply2,'points',info_site,info_plot2, info_plot1);
    if ~exist(filepath_ply2t,'file');
        color = plyintensity2color(data_ieq2,'autumn');
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply2t);
        end
        write2ply(filepath_ply2t,xyz2t', color);
    end
end


if options_write_ply;
    % Just 1st returns
    filepath_ply1  = sprintf('%s%s_%03.0f-%02.0f.ply',path_ply1,'points',info_site,info_plot2);
    if ~exist(filepath_ply1,'file');
        color = plyintensity2color(data_ieq1,'winter');
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply1);
        end
        write2ply(filepath_ply1,[data_x1 data_y1 data_z1], color);
    end
end

%% Write tree file
% Load tree
tree = 1; % Initalize as variable
filepath_tree1 = sprintf('%s%s',path_mat1,'tree.mat');
load(filepath_tree1);
tree1 = tree;
filepath_tree2 = sprintf('%s%s',path_mat2,'tree.mat');
load(filepath_tree2);
tree2 = tree;
tree2t = tree2;

n_tree1 = numel(tree1);
n_tree2 = numel(tree2);

filepath_tree2t_mat = sprintf('%stree_%03.0f-%02.0f-%02.0f.mat',path_mat2,info_site,info_plot2, info_plot1);
if exist(filepath_tree2t_mat,'file');
    load(filepath_tree2t_mat);
else
    for t2 = 1:n_tree2;
        n_seg = numel(tree2(t2).r);
        tree2t(t2).loc = (R*tree2(t2).loc) + repmat(t,1,n_seg);
    end
    save(filepath_tree2t_mat,'tree2t')
end

filepath_tree2t_ply = sprintf('%stree_%03.0f-%02.0f-%02.0f.ply',path_ply2,info_site,info_plot2, info_plot1);
if ~exist(filepath_tree2t_ply,'file');
    tree2ply(filepath_tree2t_ply,tree2t,12);
end

filepath_tree1_ply = sprintf('%stree_%03.0f-%02.0f.ply',path_ply1,info_site,info_plot1);
if ~exist(filepath_tree1_ply,'file');
    tree2ply(filepath_tree1_ply,tree1,12);
end


%% ICP Optimization

filepath_R_icp = sprintf('%sR_icp_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);
filepath_t_icp = sprintf('%st_icp_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);

%if exist(filepath_R_icp,'file');
%    load(filepath_R_icp);
%    load(filepath_t_icp);
%else
    [ R_icp,t_icp ] = icp_refine( info_site, info_plot1, info_plot2 ); 
%    save(filepath_R_icp,'R_icp');
%    save(filepath_t_icp,'t_icp');  
%end

%% Write ply file

R_eff = R_icp*R;
t_eff = R_icp*t + t_icp;

xyz2i = (R_eff*[data_x2 data_y2 data_z2]') + repmat(t_eff,1,numel(data_x2));
%data_x2t = xyz2t(1,
if options_write_ply;
    % Just 1st returns
    filepath_ply2i  = sprintf('%s%s_%03.0f-%02.0f-%02.0f.ply',path_ply2,'points-icp',info_site,info_plot2, info_plot1);
    if ~exist(filepath_ply2i,'file');
        color = plyintensity2color(data_ieq2,'spring');
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply2i);
        end
        write2ply(filepath_ply2i,xyz2i', color);
    end
end


%% Write tree file
% Load tree
tree = 1; % Initalize as variable
filepath_tree1 = sprintf('%s%s',path_mat1,'tree.mat');
load(filepath_tree1);
tree1 = tree;
filepath_tree2 = sprintf('%s%s',path_mat2,'tree.mat');
load(filepath_tree2);
tree2 = tree;
tree2t = tree2;

n_tree1 = numel(tree1);
n_tree2 = numel(tree2);

filepath_tree2t_mat = sprintf('%stree_%03.0f-%02.0f-%02.0f.mat',path_mat2,info_site,info_plot2, info_plot1);
if exist(filepath_tree2t_mat,'file');
    load(filepath_tree2t_mat);
else
    for t2 = 1:n_tree2;
        n_seg = numel(tree2(t2).r);
        tree2t(t2).loc = (R*tree2(t2).loc) + repmat(t,1,n_seg);
    end
    save(filepath_tree2t_mat,'tree2t')
end

filepath_tree2t_ply = sprintf('%stree_%03.0f-%02.0f-%02.0f.ply',path_ply2,info_site,info_plot2, info_plot1);
if ~exist(filepath_tree2t_ply,'file');
    tree2ply(filepath_tree2t_ply,tree2t,12);
end

filepath_tree1_ply = sprintf('%stree_%03.0f-%02.0f.ply',path_ply1,info_site,info_plot1);
if ~exist(filepath_tree1_ply,'file');
    tree2ply(filepath_tree1_ply,tree1,12);
end
end

