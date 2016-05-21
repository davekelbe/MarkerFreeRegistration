function [ aux ] = register_prep( path_up, info_experiment, info_suffix, info_site_curr, info_plot_register )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

info_slash = '/';
%info_site = info_sites(s);

n_plot = numel(info_plot_register);
%info_site_curr = info_sites(s);
info_valid_plot = cell(n_plot,1);
for p = 1:n_plot;
    info_valid_plot{p} = sprintf('%02.0f', info_plot_register(p)); % changed from just p ?
end
aux.info_valid_plot = info_valid_plot;
aux.info_site = info_site_curr;


P_LCS = cell(n_plot,1);
P_rad = cell(n_plot,1);
for p = 1:n_plot;
    path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
        info_experiment, info_slash, info_suffix,info_slash,info_site_curr, info_slash,info_plot_register(p),info_slash);
    path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
    path_exp = sprintf('%s%s%s',path_up,info_experiment, info_slash);
    
    
    filepath_seg_PLCS = sprintf('%s%s%s',path_mat,'seg_PLCS','.mat');
    filepath_seg_Prad = sprintf('%s%s%s',path_mat,'seg_Prad','.mat');
    load(filepath_seg_PLCS);
    load(filepath_seg_Prad);
    
    P_LCS{p} = seg_PLCS';
    P_rad{p} = seg_Prad;
    
end

%P_LCS = P_LCS;
%P_rad = P_rad; 

% Load output from stem detection
filepath_P_LCS = sprintf('%s%s%s%03.0f%s%03.0f%s%s', path_exp, info_suffix, info_slash, info_site_curr, ...
    '/', info_plot_register(1),'/mat/','P_LCS.mat');
filepath_P_rad = sprintf('%s%s%s%03.0f%s%03.0f%s%s', path_exp, info_suffix, info_slash, info_site_curr, ...
    '/', info_plot_register(1),'/mat/','P_rad.mat');
save(filepath_P_LCS,'P_LCS');
save(filepath_P_rad,'P_rad');
aux.P_LCSn = P_LCS;
aux.P_radn = P_rad;

path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
    info_experiment, info_slash, info_suffix,info_slash,info_site_curr, info_slash,info_plot_register(1),info_slash);
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);

aux.path_mat = path_mat;

aux.info_plot = info_plot_register; 
end

