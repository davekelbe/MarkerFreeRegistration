function [ R,t ] = icp_refine( info_site, info_plot1, info_plot2 )
%UNTITLED11 Summary of this function goes here
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

% Outputs:
%       options               - user-defined run-type options
%       info                  - experiment-specific info

%% Load data

path_mat1 = sprintf('%s%s%s%s',path_top1,'mat',info_slash);
path_mat2 = sprintf('%s%s%s%s',path_top2,'mat',info_slash);

% Filepaths for saved matlab variables
filepath_tree2t_mat = sprintf('%stree_%03.0f-%02.0f-%02.0f.mat',path_mat2,info_site,info_plot2, info_plot1);
load(filepath_tree2t_mat);
tree2 = tree2t;

filepath_tree1_mat = sprintf('%stree.mat',path_mat1);
tree= 0;
load(filepath_tree1_mat);
tree1 = tree;

filepath_m1_best = sprintf('%sm1_best_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);
filepath_m2_best = sprintf('%sm2_best_%03.0f-%02.0f-%02.0f.mat',path_mat2, info_site, info_plot2,info_plot1);
load(filepath_m1_best);
load(filepath_m2_best);

%% Resample to finer resolution

n_tree1 = numel(tree1);
n_tree2 = numel(tree2);

points1x = [];
points1y = [];
points1z = [];
is_bound = [];
t_resamp = 50;
for t1 = 1:n_tree1;
    if ~any(ismember(m1_best,t1))
        continue
    end
    points1x_tree = [];
    points1y_tree = [];
    points1z_tree = [];
    n_seg = numel(tree1(t1).r);
    for s = 1:n_seg-1;
    newx = linspace(tree1(t1).loc(1,s),tree1(t1).loc(1,s+1),t_resamp);
    newy = linspace(tree1(t1).loc(2,s),tree1(t1).loc(2,s+1),t_resamp);
    newz = linspace(tree1(t1).loc(3,s),tree1(t1).loc(3,s+1),t_resamp);
    points1x_tree = [points1x_tree newx(1:end-1)];
    points1y_tree = [points1y_tree newy(1:end-1)];
    points1z_tree = [points1z_tree newz(1:end-1)];
    end   
    new_is_bound = false([1,numel(points1z_tree)]);
    new_is_bound(1) = true;
    new_is_bound(end) = true;
    points1x = [points1x points1x_tree];
    points1y = [points1y points1y_tree];
    points1z = [points1z points1z_tree];
    is_bound = [is_bound new_is_bound];
    foo = 1; 
    %points1x = [points1x tree1(t1).loc
end 
is_bound = logical(is_bound);

points2x = [];
points2y = [];
points2z = [];
for t2 = 1:n_tree2;
    if ~any(ismember(m2_best,t2))
        continue
    end
    n_seg = numel(tree2(t2).r);
    for s = 1:n_seg-1;
    newx = linspace(tree2(t2).loc(1,s),tree2(t2).loc(1,s+1),t_resamp);
    newy = linspace(tree2(t2).loc(2,s),tree2(t2).loc(2,s+1),t_resamp);
    newz = linspace(tree2(t2).loc(3,s),tree2(t2).loc(3,s+1),t_resamp);
    points2x = [points2x newx(1:end-1)];
    points2y = [points2y newy(1:end-1)];
    points2z = [points2z newz(1:end-1)];
    end   
    %points1x = [points1x tree1(t1).loc
end
%{
figure; 
hold on
scatter3(points1x(is_bound), points1y(is_bound),points1z(is_bound),50,'r','filled');
scatter3(points1x(~is_bound), points1y(~is_bound),points1z(~is_bound),10,'b','filled');
scatter3(points2x, points2y,points2z,10,'g','filled');
%}
%% ICP 
q = [points1x; points1y; points1z];
p = [points2x; points2y; points2z];
bound = find(is_bound);
[R,t] = icp(q,p,50,'EdgeRejection', true, 'Boundary', bound, 'Matching', 'kDtree');
pt = R*p + repmat(t,[1,size(p,2)]);
points2tx = pt(1,:);
points2ty = pt(2,:);
points2tz = pt(3,:);

%{
figure;
scatter3(points1x, points1y, points1z,10,'r','filled');
hold on
scatter3(points2x, points2y, points2z,10,'b','filled');
scatter3(points2tx, points2ty, points2tz,10,'g','filled');
axis equal 
legend('Plot 1','Plot 2','Plot 2 transformed');
    %}
end

