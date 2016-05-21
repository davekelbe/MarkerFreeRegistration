function [  ] = transform_ply( path_up, info_experiment, info_suffix, ...
    info_plot, info_site, p  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

info_slash = '/';


path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
    info_experiment, info_slash, info_suffix,info_slash,info_site, info_slash,info_plot(p),info_slash);
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
path_ply = sprintf('%s%s%s',path_top,'ply',info_slash);

path_top_WCS = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
    info_experiment, info_slash, info_suffix,info_slash,info_site, info_slash,info_plot(1),info_slash);
path_mat_WCS = sprintf('%s%s%s',path_top_WCS,'mat',info_slash);


filepath_G_R_SVD_WMF = sprintf('%s%s',path_mat_WCS, 'G_R_SVD_WMF.mat');
load(filepath_G_R_SVD_WMF);
filepath_G_t_SVD_WMF = sprintf('%s%s',path_mat_WCS, 'G_t_SVD_WMF.mat');
load(filepath_G_t_SVD_WMF);
filepath_G_RMSE_SVD = sprintf('%s%s',path_mat_WCS, 'G_RMSE_SVD.mat');
load(filepath_G_RMSE_SVD);

% Load ply
filepath_ply_in = sprintf('%spoints_full_%03.0f-%03.0f.ply',path_ply, info_site,info_plot(p));
filepath_ply_out = sprintf('%spoints_full_%03.0f-%03.0fWCS.ply',path_ply, info_site,info_plot(p));
[vertex, ~] = read_ply(filepath_ply_in);
data_x = vertex(:,1);
data_y = vertex(:,2);
data_z = vertex(:,3);
data_i = vertex(:,4:6);
n_data = numel(data_x);
data2_xyz = (G_R_SVD_WMF{p}*[data_x data_y data_z]') + ...
    repmat(G_t_SVD_WMF{p},1,n_data);
write2ply(filepath_ply_out,data2_xyz',data_i);

