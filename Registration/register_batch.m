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

aux.info_exp = 'Harvard';
aux.info_suffix = 'reg';

site_unique = unique(C{1});
n_site = numel(site_unique);

path_exp = sprintf('%s%s%s%s%s', 'D:\Users\djk2312\Documents\', aux.info_exp,'\', aux.info_suffix, '\');
%% Create P_LCS 
%{
for site = 1:n_site;
    info_site = site_unique(site);
    P_LCS = cell(n_plot,1);
    P_rad = cell(n_plot,1);
    for p = 1:n_plot;
        info_plot = plot_t(p);
        tree = [];
        path_plot = sprintf('%s%03.0f%s%02.0f%s', path_exp, info_site, ...
            '\', plot_t(p),'\mat\');
        filepath_tree = sprintf('%stree.mat',path_plot);
        if ~exist(filepath_tree,'file')
            continue
        end
        load(filepath_tree);
        n_tree = numel(tree);
        P_LCS{info_plot} = nan(3,n_tree);
        P_rad{info_plot} = nan(n_tree,1);
        for t = 1:n_tree;
            if isempty(tree(t).loc);
                continue
            end
            P_LCS{info_plot}(:,t) = tree(t).loc(:,1);
            P_rad{info_plot}(t) = tree(t).r(1);
        end
    end
        filepath_P_LCS = sprintf('%s%03.0f%s%02.0f%s%s', path_exp, info_site, ...
            '\', plot_t(p),'\mat\','P_LCS.mat');
        filepath_P_rad = sprintf('%s%03.0f%s%02.0f%s%s', path_exp, info_site, ...
            '\', plot_t(p),'\mat\','P_rad.mat');
        save(filepath_P_LCS, 'P_LCS');        
        save(filepath_P_rad, 'P_rad');
end
%}
%%
for s = 1:n_site;
    info_site = site_unique(s);
    if info_site==17||info_site==22||info_site==23;
        continue
    end
    is_site = (C{1}==info_site);
    plot_t = C{2}(is_site);
    n_plot = numel(plot_t);
    info_valid_plot = cell(n_plot,1);
    for p = 1:n_plot;
        info_valid_plot{p} = sprintf('%02.0f', p);
    end    
    aux.info_valid_plot = info_valid_plot;
    aux.info_site = info_site;
    fprintf('\nSite %d\n', info_site)
    %kelbe_registration_combine_dis3fun( aux ) % Pairwise registrations ++ 
    %kelbe_registration_saveplot( aux )
    %kelbe_registration_combine_posecons( aux ) % Replaced by Dijkstrapose
    %kelbe_registration_combine_graph( aux )
    % kelbe_registration_combine_dijkstraposecons( aux )
    %kelbe_registration_combine_dijkstraposewcs( aux ) % This one used ++
        filepath_P_LCS = sprintf('%s%03.0f%s%02.0f%s%s', path_exp, info_site, ...
            '\', 25,'\mat\','P_LCS.mat');
        filepath_P_rad = sprintf('%s%03.0f%s%02.0f%s%s', path_exp, info_site, ...
            '\', 25,'\mat\','P_rad.mat');
        load(filepath_P_LCS);
        load(filepath_P_rad);
        aux.P_LCSn = P_LCS;
        aux.P_radn = P_rad;
    aux.path_mat = sprintf('%s%03.0f%s%02.0f%s',...
            'D:\Users\djk2312\Documents\Harvard\reg\',...
            info_site, '\', 25,'\mat\');
    % Sensitivity analysis on exponential r 
        aux.sensitivity_r = true;
        kelbe_registration_combine_dijkstraposewcs_graphanalyses( aux ) % This one used
    
    % kelbe_registration_combine_dijkstraposewcs_graphanalyses( aux ) % This one used
end