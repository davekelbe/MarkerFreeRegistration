%% Registration of Forest Terrestrial Laser Scanning Data
% This is a wrapper script that runs Dave Kelbe's code for 
%   (1) stem detection,
%   (2) manual refinement, 
%   (3) feature extraction, 
%   (4) pairwise registration, and 
%   (5) multiview registration. 
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
%       Geoscience and Remote Sensing (In Press).

%%%%% Graph-based registration (Paper 3)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Analyze the ensemble of pairwise estimates to find the optimal global pose
%   parameters, which link each node (scan) to the reference node
%   For more information, refer to: 
%       Kelbe, D., van Aardt, J., Romanczyk, P., and van Leeuwen, M. (2016). 
%       Multi-view, marker-free registration of forest terrestrial laser 
%       scanner data with embedded confidence metrics. IEEE Transactions on 
%       Geonscience and Remote Sensing (In Review).


%% Setup 
% Set up variables 

path_up = '/Volumes/DJK/Asterix/Processed/';
% A string identifying the path to save output data into 
info_experiment = 'TEST2';
% A string identifying the experiment 
% Used to identify path of raw CBL data and to organize output data 
% Acceptable values include HAKT1, HAKT2, KIP3, KIP6, KIP9, KIP13, KIP18,
% KIP30, HAVO, MAT, harvard, michigan, soaproot, sanjoaquin
info_sites=31;
% An array of sites to process, e.g., =1:4; will process sites 1,2,3, and 
% then 4 sequentially in batch mode, or =1; will process just site 1. See
% README for acceptable sites for each experiment. 
info_plots = [8,13,14];
% An array of plots to process, e.g., =1:25; will process plots 1 through
% 25 in batch mode; while =1; will process only site one. Note that the
% registration steps require multiple plots (obviously), while the stem
% detection can be performed on a per-plot basis. See README for acceptable
% plots for each site/experiment. 
info_terminate = 'transformdata';
% A string identify when to terminate the algorithm. This script is a
% wrapper that contains subfunctions representing each building block/paper
% of my dissertation. You can choose to only do DEM extraction, for
% example, by setting the info_terminate string to 'dem'. Acceptable values
% include: 
%   load:  load lidar data and perform initial preprocessing 
%   dem:  extract digital elevation model (DEM) 
%   detect_trees:  voxel-based tree stem segment detection (most computationally intensive) 
%   manual:  user-interface to correct bad segments for challenging sites like Hawaii 
%   tiepoints:  interface to ensure continuity between (1) and (2,3) 
%   pairwise:  Generate pairwise registration correspondences  
%   graphbased:  Link all scans in a graph-based network 
%   transformdata:  Transform PLY back to WCS 

info_suffix = '02-14';
% A string specifying the date so that you can re-process the data (for 
% example, with code modifications without wiping out your previous results.
% In general, keep this constant for your full set of analyses.  
% Example acceptable value: '12-03'  
info_plot_register = [8,13,14];
% An array of plots to register, as all plots for a lare site may 
% not be able to be registered together. Example: =0:4; will register all 
% plots together for KIP3, but HAVO may need to be split up within
% contiguous regions, e.g., = [23 24 25 34 35 36 45 46 47];
options_skipseg = true;
% Set to true to skip the computationally-intensive segment detection and
% solely digitize manually. Set to false to initialize candidate segments
% automatically and then refine using manual digitization

%% Error catching
AinB = ismember(info_plot_register, info_plots);
if any(~AinB);
    error('info_plot_register must be subset of info_plot');
end
    
%% Run stem detection (Paper 1) 
n_plot = numel(info_plots);
n_site = numel(info_sites);

% Loop through each scan by site and then plot 
for s = 1:n_site;
    for p = 1:n_plot;
        info_site = info_sites(s);
        info_plot = info_plots(p);
        lidar_JVA(info_site,info_plot, info_terminate, info_experiment, info_suffix, path_up, options_skipseg);
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
