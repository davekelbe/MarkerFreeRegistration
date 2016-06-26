%% Registration of Forest Terrestrial Laser Scanning Data
% This is a wrapper script that runs Dave Kelbe's code for 
%   (1) stem detection,
%   (2) manual refinement, 
%   (3) feature extraction, 
%   (4) pairwise registration, and 
%   (5) multiview registration. 
% 
% Please refer to the README for more information.
% All output data (DEM, registered .ply point clouds, etc.) are saved to disk. 
% The Code is organized as follows:

%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   All the user modifications are made here in the first block of code

%%%%% Run stem detection (Paper 1)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For more information, refer to:
%       Kelbe, D., van Aardt, J., Romanczyk, P., Cawse-Nicholson, 
%       K., and van Leeuwen, M. (2015). Single-scan stem reconstruction 
%       using low-resolution terrestrial laser scanner data. IEEE Journal 
%       of Selected Topics in Applied Earth Observations and Remote Sensing 
%       8 (7). doi: 10.1109/JSTARS. 2015.2416001.

%%%%% Manual stem refinement (Hawaii addition)   %%%%%%%%%%%%%%%%%%%%%%%%%%
%   The user can choose to refine stem detection results for very difficult
%   forest types, e.g., Hawaii, but manually adding or removing tree
%   segments using a user interface 

%%%%% Extract tie points   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Extract tie points (stem-terrain intersection) from each detected 
%   segment to serve as inputs to the registration module 

%%%%% Pairwise registration (Paper 2)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Determine pairwise pose estimates (R,t) between all pairs of scans
%   For more information, refer to: 
%       Kelbe, D., van Aardt, J., Romanczyk, P., and van Leeuwen, M. (2016).
%       Marker-free registration of forest terrestrial laser scanner data 
%       pairs with embedded confidence metrics. IEEE Transactions on 
%       Geoscience and Remote Sensing 54 (7) DOI:10.1109/TGRS.2016.2539219. 

%%%%% Graph-based registration (Paper 3)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Analyze the ensemble of pairwise estimates to find the optimal global pose
%   parameters, which link each node (scan) to the reference node
%   For more information, refer to: 
%       Kelbe, D., van Aardt, J., Romanczyk, P., and van Leeuwen, M. 
%       Multi-view, marker-free registration of forest terrestrial laser 
%       scanner data with embedded confidence metrics. IEEE Transactions on 
%       Geonscience and Remote Sensing (In Review).


%% Setup 
% Set up variables 
if ispc();
    info_slash = '\';
else
    info_slash = '/';
end

% Read preferences file (it should be saved in the same directory 
filepath_preferences = sprintf('%s%spreferences.txt', pwd, info_slash);
if ~exist(filepath_preferences, 'file');
     filepath_preferences = uipickfiles('Type', {'*.txt', 'txt files'},'Prompt', 'Please choose the preferences.txt file if it exists or press cancel if it does not exist');
    if ~iscell(filepath_preferences)
     if filepath_preferences == 0;
        filepath_preferences = [];
     end
    end
end

[ path_up, path_source, info_experiment, info_sites, info_plots, info_terminate, info_suffix, options_skipseg ] = setup( filepath_preferences );
clear filepath_preferences
%% Run stem detection (Paper 1) 
n_plot = numel(info_plots);
n_site = numel(info_sites);

% Loop through each scan by site and then plot 
% 
for s = 1:n_site;
    for p = 1:n_plot;
        info_site = info_sites(s);
        info_plot = info_plots(p);
        lidar_JVA(info_site,info_plot, info_terminate, info_experiment, info_suffix, path_up, path_source, options_skipseg);
        % This has the following subfunctions:   
        %   load 
        %   dem 
        %   detect_trees
        %   tree* Note: this function was removed and added as a
        %   subsequent block to allow manual stem refinement 
    end
end

% Check terminate string to determine execution return/continue
if strcmp(info_terminate, 'load') || strcmp(info_terminate, 'dem') || strcmp(info_terminate, 'detect_trees');
    return
end

%% Manual stem refinement (Hawaii addition) 

%
% Loop through each scan by site and then plot 
for s = 1:n_site;
    for p = 1:n_plot;
        info_site = info_sites(s);
        info_plot = info_plots(p);
        manual_stem_selection(path_up, info_site,info_plot, info_experiment, info_suffix);
    end
end

% Check terminate string to determine execution return/continue
if strcmp(info_terminate, 'manual');
    return
end
%}

%% Extract tie points (to link output from (1) with input of (2,3))
info_plot_register = info_plots; 

% Restructure inputs for continuity with registration sequence 
% Note: registration is applied at the site level, collecting all plots
extract_tiepoints( path_up, info_experiment, info_suffix, ...
    info_plot_register, info_sites, info_plots  );

% Check terminate string to determine execution return/continue
if strcmp(info_terminate, 'tiepoints');
    return
end
%% Pairwise Registration (Paper 2)

% site-by-site basis 
for s = 1:n_site;
    % Preparation work 
    info_site = info_sites(s);
    [ aux ] = register_prep( path_up, info_experiment, info_suffix, info_site, ...
        info_plot_register );
    % Pairwise registration 
    kelbe_registration_combine_dis3fun_JVA( aux )
end

% Check terminate string to determine execution return/continue
if strcmp(info_terminate, 'pairwise');
    return
end
%%  Multiview Registration (Paper 3)

% site-by-site basis 
for s = 1:n_site;
    % Preparation work
    [ aux ] = register_prep( path_up, info_experiment, info_suffix, info_site, ...
        info_plot_register );
    % Graph-based registration 
    kelbe_registration_combine_dijkstraposewcs_JVA( aux )
end

% Check terminate string to determine execution return/continue
if strcmp(info_terminate, 'graphbased');
    return
end

%%  Transform PLY based on registration 

% For each scan
n_plot_register = numel(info_plot_register);
for s = 1:n_site;
    for p = 1:n_plot_register;
    % Transform PLY into registered coordinate system
    transform_ply( path_up, info_experiment, info_suffix, ...
    info_plot_register, info_sites(s),p);
    end
end

% Check terminate string to determine execution return/continue
if strcmp(info_terminate, 'transformdata');
    return
end
