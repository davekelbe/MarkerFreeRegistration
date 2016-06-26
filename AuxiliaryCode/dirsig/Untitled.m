% Unify stem models for Martin 

options_showfig = true;
aux.info_exp = 'reg';
path_exp = 'D:\Users\djk2312\Documents\Harvard\reg\';
filename_csv = 'D:\Users\djk2312\Documents\thisshouldbecyclone\2012-08-Harvard.csv';
fid = fopen(filename_csv);
%C = textscan(fid, '%u8%u8%s%s%s%u16%u16%u16%u16%u16%u16%u16%u16%u16',...
%    'delimiter', ',',...
%    'headerlines',8);
C = textscan(fid, '%u8%u8%s%s%s%s%s%s%s%s%s%s%s%s',...
    'delimiter', ',',...
    'headerlines',8);
% 1 Site
% 2 Plot
% 3 Lidar Filename
% 4 Date
% 5 Lidar Orientation

n_scans = numel(C{1});
warning('off', 'arguments:exteriordata');

site_unique = unique(C{1});
n_site = numel(site_unique);
info_site = 31;
site = site_unique==info_site;
    info_site = site_unique(site);
    plot_t = C{2}(is_site);
    n_plot = numel(plot_t);
    info_valid_plot = cell(n_plot,1);
    for p = 1:n_plot;
        info_valid_plot{p} = sprintf('%02.0f', p);
    end
    path_reg = sprintf('%s%03.0f%s', path_exp, info_site, '\25\mat\');
    filepath_R_MST = sprintf('%sG_R_MST.mat',path_reg);
    filepath_t_MST = sprintf('%sG_t_MST.mat',path_reg);
    for p = 1:n_plot;
        % Load DEM 
        % Load Stem Model 
        path_ply = sprintf('%s%03.0f%s%%02.0f%s', path_exp, info_site, '\', info_plot,'\ply\');
        filepath_tree = sprintf('%stree_%03.0f-%02.0f.ply',path_ply, info_site, info_plot);
    foo = 1;
    end
    
        