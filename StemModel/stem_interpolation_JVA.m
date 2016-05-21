function [  ] = stem_interpolation_JVA(path_up,info_site, info_plot, info_exp, info_suffix )

% This version just computes intersection with terrain for generating tie
% poitns. No interpolation of full stem model is performed. 

info_slash = '/';
path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
    info_exp, info_slash, info_suffix,info_slash,info_site, info_slash,info_plot,info_slash);
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
path_ply = sprintf('%s%s%s',path_top,'ply',info_slash);
path_png = sprintf('%s%s%s',path_top,'png',info_slash);

filepath_seg_row = sprintf('%s%s%s',path_mat,'seg_row','.mat');
filepath_seg_col = sprintf('%s%s%s',path_mat,'seg_col','.mat');
filepath_seg_iter = sprintf('%s%s%s',path_mat,'seg_iter','.mat');
filepath_seg_z = sprintf('%s%s%s',path_mat,'seg_z','.mat');
filepath_seg_y = sprintf('%s%s%s',path_mat,'seg_y','.mat');
filepath_seg_r = sprintf('%s%s%s',path_mat,'seg_r','.mat');
filepath_seg_index = sprintf('%s%s%s',path_mat,'seg_index','.mat');
filepath_seg_fill = sprintf('%s%s%s',path_mat,'seg_fill','.mat');
filepath_seg_taper = sprintf('%s%s%s',path_mat,'seg_taper','.mat');
filepath_seg_lean = sprintf('%s%s%s',path_mat,'seg_lean','.mat');
filepath_seg_anorm = sprintf('%s%s%s',path_mat,'seg_anorm','.mat');


filepath_I12ieq = sprintf('%s%s_%03.0f-%03.0f%s',path_png,'I12ieq',info_site,info_plot,'.png');
filepath_I12r = sprintf('%s%s%s',path_mat,'I12r','.mat');
filepath_axis_a = sprintf('%s%s%s',path_mat,'axis_a','.mat');
filepath_axis_e = sprintf('%s%s%s',path_mat,'axis_e','.mat');

filepath_dem_qx = sprintf('%s%s%s',path_mat,'dem_qx','.mat');
filepath_dem_qy = sprintf('%s%s%s',path_mat,'dem_qy','.mat');
filepath_dem_qz = sprintf('%s%s%s',path_mat,'dem_qz','.mat');

load(filepath_dem_qx);
load(filepath_dem_qy);
load(filepath_dem_qz);


load(filepath_seg_row);
load(filepath_seg_col);
load(filepath_seg_iter);
load(filepath_seg_z);
load(filepath_seg_y);
load(filepath_seg_r);
load(filepath_seg_index);
load(filepath_seg_fill);
load(filepath_seg_taper);
load(filepath_seg_lean);
load(filepath_seg_anorm);

load(filepath_I12r);
I12ieq = imread(filepath_I12ieq);
load(filepath_axis_a);
load(filepath_axis_e);












%% old code 
is_valid = true(size(seg_anorm));

filepath_tree = sprintf('%s%s%s',path_mat,'tree','.mat');

t_rsearch = 0.15;
tree  = seg2tree_Jan( seg_z, seg_y, seg_r, seg_anorm, seg_iter,...
    is_valid,t_rsearch, path_mat );
save(filepath_tree,'tree');

%
%path_local = 'Z:\Desktop\Local\';
filepath_tree_ply = sprintf('%stree_%03.0f-%02.0f.ply',path_ply,info_site,info_plot);
if ~exist(filepath_tree_ply,'file');
    tree2ply(filepath_tree_ply,tree,20);
end
%}

filepath_tree_row = sprintf('%s%s%s',path_mat,'tree_row','.mat');
filepath_tree_col = sprintf('%s%s%s',path_mat,'tree_col','.mat');
filepath_tree_iter = sprintf('%s%s%s',path_mat,'tree_iter','.mat');
if exist(filepath_tree_row,'file');
    load(filepath_tree_row);
    load(filepath_tree_col);
    load(filepath_tree_iter);
else
    [ tree_row, tree_col,tree_iter ] = andrieu_tree_edges( I12r, axis_a, axis_e, tree );
    save(filepath_tree_row, 'tree_row');
    save(filepath_tree_col, 'tree_col');
    save(filepath_tree_iter, 'tree_iter');
end

%filepath_tree_final = sprintf('%stree_final_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
filepath_tree_final = sprintf('%stree_final_%03.0f-%03.0f.png',path_png,info_site,info_plot);
if ~exist(filepath_tree_final,'file');
    figure;
    imagesc(I12ieq);
    hold on
    for s = 1:numel(tree_col);
        plot([tree_col{s} tree_col{s}(1)],[tree_row{s} tree_row{s}(1)],'-r')
    end
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_tree_final,'png');
end
close all
end
