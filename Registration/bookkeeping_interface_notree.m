function [ aux ] = bookkeeping_interface_notree( path_up, info_experiment, info_suffix, plot_register, info_sites, info_plot  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_site = numel(info_sites);
aux.info_exp = info_experiment;
aux.info_suffix = info_suffix;
path_exp = sprintf('%s%s%s%s%s', path_up, aux.info_exp,'/', aux.info_suffix, '/');
n_plot = numel(plot_register);
for s = 1:n_site;
    info_site_curr = info_sites(s);
    % Initialize variable to store registration input
    P_LCS = cell(n_plot,1);
    P_rad = cell(n_plot,1);
    for p = 1:n_plot;
        info_plot = plot_register(p);
        tree = [];
        path_plot = sprintf('%s%03.0f%s%03.0f%s', path_exp, info_site_curr, ...
            '/', plot_register(p),'/mat/');
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
    filepath_P_LCS = sprintf('%s%03.0f%s%03.0f%s%s', path_exp, info_site_curr, ...
        '/', plot_register(p),'/mat/','P_LCS.mat');
    filepath_P_rad = sprintf('%s%03.0f%s%03.0f%s%s', path_exp, info_site_curr, ...
        '/', plot_register(p),'/mat/','P_rad.mat');
    save(filepath_P_LCS, 'P_LCS');
    save(filepath_P_rad, 'P_rad');
end

end

