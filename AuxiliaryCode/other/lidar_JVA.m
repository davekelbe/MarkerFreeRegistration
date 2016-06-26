function [  ] = lidar_JVA(info_site, info_plot, info_terminate, info_exp, info_suffix, path_up, path_source, options_skipseg )
%LIDAR_PROCESSING Master code for lidar processing
%
% ********** Input ***************
%
% info_site             - site number
% info_plot             - plot number
% info_terminate        - string which defines the last cell evaluated
% info_site             - string for experiment (harvard, california, etc)

%
% ********** Output ***************
%
% There are no outputs to this function
% Data are procesed and stored into appropriate folders during execution
% Data are stored in according to the following folder structure:
%   <path_up><info_exp>/<info_suffix>/<info_site>/<info_plot>/<'png','mat','ply',...>
%
% ********** Cells ****************
%
% 'initialization' - Preparatory work
% info_exp              - Experiment name [string]
% info_suffix           - Clarifier, e.g., date or run [string]
%
% 'load' - Reads a lidar txt file and parses output to variables
% t_max_r               - Maximum xy-range [m]
% t_min_r               - Minimum xyz-range (due to sensor scattering) [m]
%
% 'subsample' - Subsamples a point cloud to give uniform density
% ** NOT USED RIGHT NOW
% t_boxr                - Box size for reducing point cloud
% t_n_max               - Number of points in each box
%
% 'dem' - Produces the Digital Elevation Model (DEM)
% t_azimuth_bins        - Azimuthal resampling for Delaunay Triangulation
% t_angle_deviation     - Maximum face normal (w.r.t. mean face normal)
% t_grid_size           - Gridded DEM spacing
%
% 'detect_trees'
% t_r = 2;              - Size of voxel [m]
% t_min_density = 0.1;  - percent of expected # points for density removal
% t_coverage = 0.5;     - Required coverage [%]
%
% p_theta = 20;         - Maximum cylinder lean angle [deg]
% p_r1_max = 1;         - Maximum radius [m]
% p_zmax = 15;          - Crown base height
% p_zmin = .5;          - Ground vegetation height
% p_taper_min = -5;     - Minimum taper angle
% p_taper_max = 5;      - Maximum taper angle
%
% *********** Example *************
%
% lidar(31,13, 'detect_trees');
%
% *********** Copyright ***********
%
% Copyright (C) David Kelbe, 2013 Rochester Institute of Technology

%% Initialization

close all

% Print current status
%{
fprintf('\n**************\n');
fprintf('info_site = %3.0f \n', info_site);
fprintf('info_plot = %3.0f \n', info_plot);
%}

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
options_ptsize = 5;
options_verbose = false;
options_verbose_progress = false;
options_verbose_output = false;
options_show_fig = false;
%options_save_eps = false;
%options_save_fig = false;
%options_save_png = false;
options_write_ply = false;

if ispc();
    info_slash = '\';
else
    info_slash = '/';
end

% Make directories
path_common = sprintf('%s%s%s%s%s',path_up,...
    info_exp, info_slash, 'Common',info_slash);
path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
    info_exp, info_slash, info_suffix,info_slash,info_site, info_slash,info_plot,info_slash);
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
path_fig = sprintf('%s%s%s',path_top,'fig',info_slash);
path_eps = sprintf('%s%s%s',path_top,'eps',info_slash);
path_ply = sprintf('%s%s%s',path_top,'ply',info_slash);
path_png = sprintf('%s%s%s',path_top,'png',info_slash);
% Execute UNIX command
if ~exist(path_common,'dir'); mkdir(path_common); end
if ~exist(path_top,'dir'); mkdir(path_top); end
if ~exist(path_mat,'dir'); mkdir(path_mat); end
if ~exist(path_fig,'dir'); mkdir(path_fig); end
if ~exist(path_eps,'dir'); mkdir(path_eps); end
if ~exist(path_ply,'dir'); mkdir(path_ply); end
if ~exist(path_png,'dir'); mkdir(path_png); end

% Outputs:
%       options               - user-defined run-type options
%       info                  - experiment-specific info
%% Load data

% Thresholds
t_max_r = 100; %11.34*sqrt(2); %11.34 is a standard forestry plot size
t_min_r = .25; % Scattering near sensor

% Print current status
if options_verbose;
    fprintf('\nLoad data\n');
    fprintf('t_maxr = %3.2f m\n', t_max_r);
    fprintf('t_origin = %3.2f m\n', t_min_r);
end

filename_source = sprintf('%s-%03.0f-%03.0f', info_exp, info_site, info_plot);

if numel(filename_source)==0;
    return
end
filepath_source = sprintf('%s%s.txt',path_source,filename_source);

% Filepaths for saved matlab variables
filepath_data_i  = sprintf('%s%s',path_mat,'data_i.mat');
filepath_data_ieq  = sprintf('%s%s',path_mat,'data_ieq.mat');
%filepath_data_index  = sprintf('%s%s',path_mat,'data_index.mat');
filepath_data_x  = sprintf('%s%s',path_mat,'data_x.mat');
filepath_data_y  = sprintf('%s%s',path_mat,'data_y.mat');
filepath_data_z  = sprintf('%s%s',path_mat,'data_z.mat');
filepath_data_n  = sprintf('%s%s',path_mat,'data_n.mat');
filepath_data_xy = sprintf('%s%s',path_mat,'data_xy.mat');
%filepath_data_r  = sprintf('%s%s',path_mat,'data_r.mat');
filepath_data_a  = sprintf('%s%s',path_mat,'data_a.mat');
filepath_data_e  = sprintf('%s%s',path_mat,'data_e.mat');
% If matlab variables exist, just load them
if exist(filepath_data_ieq,'file');
    if options_verbose_progress;
        fprintf('\n %% Loading data \n');
    end
    load(filepath_data_i);
    load(filepath_data_ieq);
    %load(filepath_data_index);
    load(filepath_data_x);
    load(filepath_data_y);
    load(filepath_data_z);
    load(filepath_data_n);
    load(filepath_data_xy);
    %load(filepath_data_r);
    load(filepath_data_a);
    load(filepath_data_e);
    % If matlab variables do not exist, produce them
elseif ~exist(filepath_data_ieq, 'file');
    % Read raw text file
    if options_verbose_progress;
        fprintf('\n %% Reading data \n');
    end
    [data_i, data_z, data_n, data_r, data_a, data_e]  = ...
        read_lidar(filepath_source);
   % if strcmp(options_las_processor, 'new');
       % data_r = data_r;
       % data_z = data_z;
   % end
    % If empty, return
    if numel(data_z)==0;
        return
    end
    %data_index = (1:numel(data_z))';
    % Convert scan angle & rotation stage to comply with Matlab
    [data_a, data_e] = az_to_360(data_a,data_e, data_z);
    data_a = -(data_a-pi/2);
    data_a(data_a<0) = data_a(data_a<0) + 360;
    % Compute x,y,z
    [data_x,data_y, data_z] = sph2cart(deg2rad(data_a),deg2rad(data_e),data_r); % Note coordinate flip
    data_xy = sqrt((data_x).^2+(data_y).^2);
    % Crop based on min and max range parameters
    index_far = (data_xy<t_max_r);
    index_near = (data_r>t_min_r);
    index = index_far&index_near;
    data_i = data_i(index);
    %data_index = data_index(index);
    data_x = data_x(index);
    data_y = data_y(index);
    data_z = data_z(index);
    data_n = data_n(index);
    data_xy = data_xy(index);
    %data_r = data_r(index);
    data_a = data_a(index);
    data_e = data_e(index);
    % Intensity Stretching
    data_ieq = data_i;
    Y = prctile(data_i,[2,98]);
    data_ieq(data_ieq<Y(1))=Y(1);
    data_ieq(data_ieq>Y(2))=Y(2);
    %{
     [I12ieq,~,~] = andrieu_image_from_points2(...
        data_a, data_e, data_n, data_ieq);
     figure; imagesc(I12ieq); colormap('gray'); axis image
    %}
    if options_verbose_progress;
        fprintf('\n %% Saving data \n');
    end
    save(filepath_data_i, 'data_i');
    save(filepath_data_ieq, 'data_ieq');
    %save(filepath_data_index, 'data_index');
    save(filepath_data_x, 'data_x');
    save(filepath_data_y, 'data_y');
    save(filepath_data_z, 'data_z');
    save(filepath_data_n, 'data_n');
    save(filepath_data_xy,'data_xy');
    %save(filepath_data_r, 'data_r');
    save(filepath_data_a, 'data_a');
    save(filepath_data_e, 'data_e');
    clear index_far index_near index Y
end

n_pts = numel(data_x);

% Index of n_plot # points for memory efficient visualization
%{
n_plot = 100000;
temp = randperm(n_pts)';
index_plot = temp(1:n_plot);
%}

% Plot the point cloud
if options_show_fig;
    figure('Color', 'w');
    handles_scatterdata = scatter3(data_x,data_y,data_z,options_ptsize,data_z, 'filled');
    axis equal;
    axis([-t_max_r t_max_r -t_max_r t_max_r min(data_z) max(data_z)]);
    grid on; hold on; %view(options_plotview);
    h1 = gcf;
end;

if true;%options_write_ply
    % Full point cloud
    filepath_ply_full  = sprintf('%s%s_%03.0f-%03.0f.ply',path_ply,'points_full',info_site,info_plot);
    if ~exist(filepath_ply_full,'file');
        color = plyintensity2color(data_ieq,'gray');
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply_full);
        end
        write2ply(filepath_ply_full,[data_x data_y data_z], color);
    end
end

if options_write_ply
    % Full dark point cloud
    filepath_ply_full_dark  = sprintf('%s%s_%03.0f-%02.0f.ply',path_ply,'points_full_dark',info_site,info_plot);
    if ~exist(filepath_ply_full_dark,'file');
        color = plyintensity2color(data_ieq,'gray');
        color = color./255;
        color = brighten(color,-.5);
        color = color.*255;
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply_full_dark);
        end
        write2ply(filepath_ply_full_dark,[data_x data_y data_z], color);
    end
end

if options_write_ply;
    % Just 1st returns
    filepath_ply_1return  = sprintf('%s%s_%03.0f-%02.0f.ply',path_ply,'points_1return',info_site,info_plot);
    if ~exist(filepath_ply_1return,'file');
        n1 = (data_n==1);
        color1 = plyintensity2color(data_ieq(n1),'gray');
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply_1return);
        end
        write2ply(filepath_ply_1return,[data_x(n1) data_y(n1) data_z(n1)], color1);
    end
end

if options_write_ply;
    % Just 2nd Returns
    filepath_ply_2return  = sprintf('%s%s_%03.0f-%02.0f.ply',path_ply,'points_2return',info_site,info_plot);
    if ~exist(filepath_ply_2return,'file');
        n2 = (data_n==2);
        color2 = plyintensity2color(data_ieq(n2),'gray');
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply_2return);
        end
        write2ply(filepath_ply_2return,[data_x(n2) data_y(n2) data_z(n2)], color2);
    end
end

if options_verbose_output;
    fprintf('\n');
    if options_write_ply;
        fprintf('Output PLY: %s \n', filepath_ply_full);
        fprintf('Output PLY: %s \n', filepath_ply_full_dark);
        fprintf('Output PLY: %s \n', filepath_ply_1return);
        fprintf('Output PLY: %s \n', filepath_ply_2return);
    end
end

%clear data_i
clear path_source filepath_source filename_source
clear temp Y
clear filepath_data_x filepath_data_y filepath_data_z filepath_data_xy
clear filepath_data_n filepath_data_i filepath_data_ieq
clear filepath_data_i filepath_data_a filepath_data_e
if strcmp(info_terminate, 'load');
    return
end
% Outputs:
%       data_x                - x values from lidar
%       data_y                - y values from lidar
%       data_z                - z values from lidar
%       data_n                - Return number from lidar
%       data_ieq              - histogram equalized intensity from lidar
%       data_a                - azimuth values from lidar
%       data_e                - elevation values from lidar
%% Andrieu Projection Images

% Define filenames
filepath_I12ieq = sprintf('%s%s_%03.0f-%03.0f%s',path_png,...
    'I12ieq',info_site,info_plot,'.png');
axis_a = sort(unique(data_a)); % 1:numel(axis_a);
axis_e = sort(unique(data_e),'descend'); %4*(axis_e + 45) %Changed from flipud

filepath_axis_a = sprintf('%s%s%s',path_mat,'axis_a','.mat');
filepath_axis_e = sprintf('%s%s%s',path_mat,'axis_e','.mat');
save(filepath_axis_a, 'axis_a');
save(filepath_axis_e, 'axis_e');


if ~exist(filepath_I12ieq,'file');
    [I12ieq,I1ieq,I2ieq] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_ieq);
    n_cmap = round(max(I12ieq(:))*1000);
    cmap = gray(n_cmap);
    if options_show_fig
        figure;
    else
        figure('visible','off');
    end
    h = imagesc(I12ieq);
    colormap(cmap);
    freezeColors;
    I12ieq = get(h,'CData');
    imwrite(I12ieq, filepath_I12ieq);
else
    I12ieq = imread(filepath_I12ieq);
end

% Create Andrieu Range image
filepath_I12r = sprintf('%s%s%s',path_mat,'I12r','.mat');
data_r = sqrt(data_x.^2 + data_y.^2 + data_z.^2);
axis_a = sort(unique(data_a)); % 1:numel(axis_a);
axis_e = sort(unique(data_e),'descend'); %4*(axis_e + 45) %Changed from flipud
if ~exist(filepath_I12r,'file');
    [I12r,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_r);
    save(filepath_I12r,'I12r')
else
    load(filepath_I12r);
end
% Create Andrieu xy image
filepath_I12xy = sprintf('%s%s%s',path_mat,'I12xy','.mat');
if ~exist(filepath_I12xy,'file');
    [I12xy,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_xy);
    save(filepath_I12xy,'I12xy')
else
    load(filepath_I12xy);
end

% Create Andrieu Range image
filepath_I12e = sprintf('%s%s%s',path_mat,'I12e','.mat');
if ~exist(filepath_I12e,'file');
    [I12e,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_e);
    save(filepath_I12e,'I12e')
else
    load(filepath_I12e);
end
filepath_I12x = sprintf('%s%s%s',path_mat,'I12x','.mat');
if ~exist(filepath_I12x,'file');
    [I12x,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_x);
    save(filepath_I12x,'I12x')
else
    load(filepath_I12x);
end
filepath_I12y = sprintf('%s%s%s',path_mat,'I12y','.mat');
if ~exist(filepath_I12y,'file');
    [I12y,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_y);
    save(filepath_I12y,'I12y')
else
    load(filepath_I12y);
end
filepath_I12z = sprintf('%s%s%s',path_mat,'I12z','.mat');
if ~exist(filepath_I12z,'file');
    [I12z,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_z);
    save(filepath_I12z,'I12z')
else
    load(filepath_I12z);
end

clear filepath_I12ieq n_cmap cmap h
% Output
% axis_a                - Vector of azimuth values for Andrieu Projection
% axis_e                - Vector of elevation values for Andrieu Projection
% I12ieq                - Andrieu intensity image
%% Subsampling
%{
t_boxr = 2; % Subsampling box size for reducing point cloud
t_n_max = 300; % Subsampling numbrer of points in each box

if options_verbose;
    fprintf('\nSubsampling\n');
    fprintf('t_boxr = %f m\n', t_boxr);
    fprintf('t_n_max = %f\n', t_n_max);
end

if exist(filepath_logical_sub,'file');
    load(filepath_logical_sub);
elseif ~exist(filepath_logical_sub, 'file');
    logical_sub = subsetLidarPoints(data_x, data_y, data_z, t_boxr, t_n_max);
    save(filepath_logical_sub, 'filepath_logical_sub');
end

if options_write_ply
    filepath_ply_subset  = sprintf('%s%s',path_ply,'points_subset.ply');
    write2ply(filepath_ply_subset,...
        [data_x(logical_sub) data_y(logical_sub) data_z(logical_sub)], ...
        color(logical_sub,:));
end

if strcmp(terminate, 'subsample');
    return
end
%}
%% DEM
% digital_elevation_model2
%   radial_subsample
%   point removal

t_azimuth_bins = 16; % Azimuthal resampling for Delaunay Triangulation
t_angle_deviation = 15; % Maximum face normal (w.r.t. mean face normal)
t_grid_size = 0.5; % Gridded DEM spacing

if options_verbose;
    fprintf('\nDEM\n');
    fprintf('t_azimuth_bins = %2.0f degrees\n', t_azimuth_bins);
    fprintf('t_angle_deviation = %3.1f x\n', t_angle_deviation);
    %    fprintf('t_grid_size = %3.2f m\n', t_grid_size);
end

filepath_a_blocked = sprintf('%s%s',path_common,'a_blocked.mat');
filepath_e_blocked = sprintf('%s%s',path_common,'e_blocked.mat');
filepath_xy_blocked = sprintf('%s%s',path_common,'xy_blocked.mat');
if exist(filepath_a_blocked,'file');
    load(filepath_a_blocked);
    load(filepath_e_blocked);
    load(filepath_xy_blocked);
elseif ~exist(filepath_a_blocked, 'file');
    [ a_blocked, e_blocked, xy_blocked ] = lidar_stage_shadow( data_a, data_e, data_xy );
    save(filepath_a_blocked, 'a_blocked');
    save(filepath_e_blocked, 'e_blocked');
    save(filepath_xy_blocked, 'xy_blocked', '-v7.3');
end

filepath_data_z0 = sprintf('%s%s',path_mat,'data_z0.mat'); % Height above ground
filepath_data_dem  = sprintf('%s%s',path_mat,'data_dem.mat'); % Ground height
filepath_dem_qx  = sprintf('%s%s',path_mat,'dem_qx.mat');
filepath_dem_qy  = sprintf('%s%s',path_mat,'dem_qy.mat');
filepath_dem_qz  = sprintf('%s%s',path_mat,'dem_qz.mat');
if exist(filepath_data_dem,'file');
    load(filepath_data_dem);
    load(filepath_dem_qx);
    load(filepath_dem_qy);
    load(filepath_dem_qz);
    load(filepath_data_z0);
elseif ~exist(filepath_data_dem, 'file');
    [data_dem, dem_qx, dem_qy, dem_qz] = ...
        digital_elevation_model(data_x,data_y,data_z,data_xy, data_a,data_e, ...
        t_max_r, t_azimuth_bins,...
        t_angle_deviation, t_grid_size, options_write_ply, path_top,...
        a_blocked,xy_blocked, info_site, info_plot);
    data_z0 = data_z - data_dem;
    save(filepath_data_dem, 'data_dem');
    save(filepath_dem_qx, 'dem_qx');
    save(filepath_dem_qy, 'dem_qy');
    save(filepath_dem_qz, 'dem_qz');
    save(filepath_data_z0, 'data_z0');
end

if options_show_fig;
    hold on;
    handles_dem = surf(dem_qx,dem_qy,dem_qz);
end

%
%Written in dem function
if true
    filepath_ply_dem = sprintf('%sply_dem_%03.0f-%03.0f.ply',path_ply,info_site,info_plot);
    %   filepath_ply_dem  = sprintf('%s%s',path_ply,'dem', );
    if ~exist(filepath_ply_dem,'file');
        [ vertices, indices ] = ply_qxqyqz_tri( dem_qx,dem_qy,dem_qz );
        plycolor = vec2cmap(vertices(:,3),'jet');
        write2plyfaces_2( filepath_ply_dem, vertices(1:4225,:), plycolor(1:4225,:), indices );
    end
end
%}

% Create Andrieu images
filepath_I12z0 = sprintf('%s%s%s',path_mat,'I12z0','.mat');
if ~exist(filepath_I12z0,'file');
    [I12z0,~,~] = andrieu_image_from_points(...
        data_a, data_e, data_n, data_z0);
    save(filepath_I12z0,'I12z0')
else
    load(filepath_I12z0);
end


if strcmp(info_terminate, 'dem');
    return
end

clear dem_qx dem_qy dem_qz
clear t_azimuth_bins t_angle_deviation t_grid_size
clear filepath_data_dem filepath_data_z0
%clear filepath_dem_qx filepath_dem_qy filepath_dem_qz
clear filepath_a_blocked filepath_e_blocked filepath_xy_blocked
clear a_blocked e_blocked xy_blocked
% Output
% axis_a                - Vector of azimuth values for Andrieu Projection
% axis_e                - Vector of elevation values for Andrieu Projection
% I12ieq                - Andrieu intensity image
% data_z0               - Height of each point above ground
%% Tree candidate points
% Prepare filepaths
filepath_seg_row = sprintf('%s%s',path_mat,'seg_row.mat');
filepath_seg_col = sprintf('%s%s',path_mat,'seg_col.mat');

filepath_seg_xb = sprintf('%s%s',path_mat,'seg_xb.mat');
filepath_seg_yb = sprintf('%s%s',path_mat,'seg_yb.mat');
filepath_seg_z = sprintf('%s%s',path_mat,'seg_z.mat');
filepath_seg_y = sprintf('%s%s',path_mat,'seg_y.mat');
filepath_seg_r = sprintf('%s%s',path_mat,'seg_r.mat');
filepath_seg_a = sprintf('%s%s',path_mat,'seg_a.mat');
filepath_seg_iter = sprintf('%s%s',path_mat,'seg_iter.mat');
filepath_seg_index = sprintf('%s%s',path_mat,'seg_index.mat');
filepath_seg_fill = sprintf('%s%s',path_mat,'seg_fill.mat');
filepath_seg_min_obj_size = sprintf('%s%s',path_mat,'seg_min_obj_size.mat');
filepath_seg_taper = sprintf('%s%s',path_mat,'seg_taper.mat');
filepath_seg_z_left = sprintf('%s%s',path_mat,'seg_z_left.mat');
filepath_seg_y_left = sprintf('%s%s',path_mat,'seg_y_left.mat');
filepath_seg_z_right = sprintf('%s%s',path_mat,'seg_z_right.mat');
filepath_seg_y_right = sprintf('%s%s',path_mat,'seg_y_right.mat');
filepath_seg_status = sprintf('%s%s',path_mat,'seg_status.mat');
filepath_seg_lean = sprintf('%s%s',path_mat,'seg_lean.mat');
filepath_seg_anorm = sprintf('%s%s',path_mat,'seg_anorm.mat');

% Initialize in case of skip
seg_row = cell(0,0);
seg_col = cell(0,0);
seg_iter = [];
seg_z = [];
seg_y = [];
seg_r = [];
seg_index = cell(0,0);
seg_fill = [];
seg_taper = [];
seg_lean = [];
seg_anorm = [];

% Save in case of skip
save(filepath_seg_row,'seg_row');
save(filepath_seg_col,'seg_col');
save(filepath_seg_iter,'seg_iter');
save(filepath_seg_z,'seg_z');
save(filepath_seg_y,'seg_y');
save(filepath_seg_r,'seg_r');
save(filepath_seg_index,'seg_index');
save(filepath_seg_fill,'seg_fill');
save(filepath_seg_taper,'seg_taper');
save(filepath_seg_lean,'seg_lean');
save(filepath_seg_anorm,'seg_anorm');

if options_skipseg ==true
    return
end

filepath_time = sprintf('%s%s',path_mat,'time.mat');

% Set system specifications
% D:\Users\djk2312\Documents\Paper01\SICK_specs.m

% Set thresholds
t_r = 2;                % Size of voxel [m]
t_min_density = .1; % percent of expected # points for density removal
t_coverage = 0.15 ;       % Required coverage [%]

% Set parameters
p_theta = 25;           % Maximum lean angle [deg]
p_r1_max = 1;           % Maximum radius [m]
p_zmax = 15;            % Crown base height
p_zmin = .5;            % Ground vegetation height
p_taper_min = -5;       % Minimum taper angle
p_taper_max = 5;        % Maximum taper angle

if options_verbose;
    fprintf('\nTree Candidate Points\n');
    fprintf('t_zmax = %2.0f m\n', p_zmax);
    fprintf('t_zmin = %2.0f m\n', p_zmin);
    fprintf('t_r = %2.0f m\n', t_r);
end

%Images of blocks
filepath_andrieu_block_xyz = sprintf('%s%s',path_mat,'andrieu_block_xyz.mat');
filepath_xyz_b = sprintf('%s%s',path_mat,'xyz_b.mat');
if exist(filepath_andrieu_block_xyz,'file');
    load(filepath_andrieu_block_xyz);
    load(filepath_xyz_b);
else
    [I_block_xyz, xyz_b] = andrieu_block_xyz(data_x,data_y,data_z,data_a, data_e, data_n,t_r );
    save(filepath_andrieu_block_xyz,'I_block_xyz');
    save(filepath_xyz_b,'xyz_b');
end
%{
    h = figure('color','w');
    imagesc(I_block_xyz)
    cmap = jet(size(xyz_b,1)+1);
    cmap(1,:) = [0 0 0];
    colormap(cmap);
    axis off; axis image
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
%}
tic

if exist(filepath_seg_z,'file');
    load(filepath_seg_xb);
    load(filepath_seg_yb);
    load(filepath_seg_z);
    load(filepath_seg_y);
    load(filepath_seg_r);
    load(filepath_seg_a);
    load(filepath_seg_iter);
    load(filepath_seg_index);
    load(filepath_seg_fill);
    load(filepath_seg_min_obj_size);
    load(filepath_seg_taper);
    load(filepath_seg_z_left);
    load(filepath_seg_y_left);
    load(filepath_seg_z_right);
    load(filepath_seg_y_right);
    load(filepath_seg_status);
    %load(filepath_seg_lean);
    %load(filepath_seg_anorm);
    % load(filepath_time);
elseif ~exist(filepath_seg_z, 'file');
    [seg_xb, seg_yb, seg_zb, seg_z, seg_y, seg_r, seg_iter, ...
        seg_index, seg_fill, seg_taper, seg_min_obj_size...
        seg_z_left,seg_y_left, seg_z_right, seg_y_right, seg_status] = ...
        detect_trees(data_x,data_y,data_z,data_z0,data_dem, ...
        I12ieq, axis_a,axis_e, t_r, p_zmin, p_zmax, path_top,...
        p_theta,p_r1_max,p_taper_min, p_taper_max,...
        t_coverage,t_min_density);
    time = toc;
    seg_a = seg_z - seg_y;
    seg_anorm = seg_a./repmat(sqrt(sum(seg_a.*seg_a,2)),[1,3]);
    seg_lean = acosd(seg_anorm(:,3));
    seg_min_obj_size = (2*mean(seg_r,2))./seg_min_obj_size;
    save(filepath_seg_xb,'seg_xb');
    save(filepath_seg_yb,'seg_yb');
    save(filepath_seg_z,'seg_z');
    save(filepath_seg_y,'seg_y');
    save(filepath_seg_r,'seg_r');
    save(filepath_seg_a,'seg_a');
    save(filepath_seg_iter, 'seg_iter');
    save(filepath_seg_index, 'seg_index');
    save(filepath_seg_fill, 'seg_fill');
    save(filepath_seg_min_obj_size, 'seg_min_obj_size');
    save(filepath_seg_taper, 'seg_taper');
    save(filepath_seg_z_left,'seg_z_left');
    save(filepath_seg_y_left,'seg_y_left');
    save(filepath_seg_z_right,'seg_z_right');
    save(filepath_seg_y_right,'seg_y_right');
    save(filepath_seg_status,'seg_status');
    save(filepath_seg_lean,'seg_lean');
    save(filepath_seg_anorm,'seg_anorm');
    save(filepath_time,'time');
end
n_seg = size(seg_r,1);

if size(seg_r,2) == 1;
    seg_r = [seg_r seg_r];
end
% Output
% seg_xb                - x-block number
% seg_yb                - y-block number
% seg_zb                - z_block number
% seg_z                 - xyz location of lower axis point
% seg_y                 - xyz location of upper axis point
% seg_z_left            - xyz location of upper left point
% seg_z_right           - xyz location of upper right point
% seg_y_left            - xyz location of lower left point
% seg_y_right           - xyz location of lower right point
% seg_r                 - radius vector for [lower upper] segment radius
% seg_iter              - iteration number
% seg_index             - Index to the points which are inliers
% seg_fill              - Fill percentage (may be superceded?)
% seg_taper             - Taper angle of each segment
% seg_min_obj_size      - Ratio of the radius to the minimum object size

%if strcmp(info_terminate, 'detect_trees');
%    return
%end

%%
seg_a = seg_z - seg_y;
seg_anorm = seg_a./repmat(sqrt(sum(seg_a.*seg_a,2)),[1,3]);
seg_lean = acosd(seg_anorm(:,3));

%{
logical_inliers = false(n_pts,1);
for s = 1:n_seg;
    logical_inliers(seg_index{s}) = true;
end
index_inliers = find(logical_inliers);
index_outliers = find(~logical_inliers);
%}

%% Adjust segments for wrap-around
filepath_seg_col = sprintf('%s%s%s',path_mat,'seg_col','.mat');
filepath_seg_row = sprintf('%s%s%s',path_mat,'seg_row','.mat');
filepath_seg_iter = sprintf('%s%s%s',path_mat,'seg_iter','.mat');
if exist(filepath_seg_col,'file');
    load(filepath_seg_col);
    load(filepath_seg_row);
    load(filepath_seg_iter);
else
    [seg_col, seg_row, seg_iter] = andrieu_fix_wraparound(...
        I12r,axis_a,axis_e,seg_iter,...
        seg_z_left, seg_y_left, seg_z_right, seg_y_right);
    is_finite = false(numel(seg_col),1);
    for s = 1:numel(seg_col);
        is_finite(s) = ~(any(~isfinite(seg_col{s})) || any(~isfinite(seg_row{s})));
    end
    seg_iter = seg_iter(is_finite);
    seg_row = seg_row(is_finite);
    seg_col = seg_col(is_finite);
    save(filepath_seg_col,'seg_col');
    save(filepath_seg_row,'seg_row');
    save(filepath_seg_iter,'seg_iter');
end

% Update seg arrays for split polygons
seg_xb = seg_xb(seg_iter);
seg_yb = seg_yb(seg_iter);
seg_z = seg_z(seg_iter,:);
seg_y = seg_y(seg_iter,:);
seg_z_left = seg_z_left(seg_iter,:);
seg_y_left = seg_y_left(seg_iter,:);
seg_z_right = seg_z_right(seg_iter,:);
seg_y_right = seg_y_right(seg_iter,:);
seg_r = seg_r(seg_iter,:);
seg_index = seg_index(seg_iter);
seg_fill = seg_fill(seg_iter);
seg_taper = seg_taper(seg_iter);
seg_min_obj_size = seg_min_obj_size(seg_iter);
seg_lean = seg_lean(seg_iter);
seg_anorm = seg_anorm(seg_iter,:);

% Show all segments outline image
%filepath_I12ieq_outline = sprintf('%s%s%s',path_png,'I12ieq_outline','.png');
filepath_I12ieq_outline = sprintf('%sI12ieq_outline_%03.0f-%03.0f.png',path_png,info_site,info_plot);
if ~exist(filepath_I12ieq_outline);
    figure;
    imagesc(I12ieq);
    hold on
    for s = 1:numel(seg_col);
        plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-r')
    end
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_I12ieq_outline,'png');
end
close all

% Show all segments
%{
figure;
hold on;
options_cyl_alpha = 0.5;
for s = 1:n_seg;
[Coneh,End1h,End2h] = Cone(seg_z(s,:),...
            seg_y(s,:),...
            [seg_r(s,1) seg_r(s,2)],...
            30,...
            'b',1,0);
        set(Coneh,'facealpha',options_cyl_alpha)
        set(End1h,'facealpha',options_cyl_alpha)
        set(End2h,'facealpha',options_cyl_alpha)
end
%}
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

save(filepath_seg_row, 'seg_row');
save(filepath_seg_col, 'seg_col');
save(filepath_seg_iter, 'seg_iter');
save(filepath_seg_z, 'seg_z');
save(filepath_seg_y, 'seg_y');
save(filepath_seg_r, 'seg_r');
save(filepath_seg_index, 'seg_index');
save(filepath_seg_fill, 'seg_fill');
save(filepath_seg_taper, 'seg_taper');
save(filepath_seg_lean, 'seg_lean');
save(filepath_seg_anorm, 'seg_anorm');

return 

%if strcmp(info_terminate, 'manual');
%    return
%end

%{



%% Filter segments
% Compute Nearer/ Farther Image
filepath_seg_farther = sprintf('%s%s%s',path_mat,'seg_farther','.mat');
filepath_seg_inlier = sprintf('%s%s%s',path_mat,'seg_inlier','.mat');
filepath_I_behind = sprintf('%s%s%s',path_mat,'I_behind','.mat');
if exist(filepath_seg_farther,'file');
    load(filepath_seg_farther);
    load(filepath_seg_inlier);
    load(filepath_I_behind);
else
    [seg_farther,seg_inlier,I_behind] = andrieu_check_showthrough2(...
        data_x, data_y,data_z,I12r,seg_fill,...
        seg_index,seg_iter,axis_a,axis_e,[0 2],...
        seg_col, seg_row, seg_z,seg_y);
    save(filepath_seg_farther,'seg_farther');
    save(filepath_seg_inlier,'seg_inlier');
    save(filepath_I_behind,'I_behind');
end


%filepath_I12_farther = sprintf('%s%s%s',path_png,'I12_farther','.png');
filepath_I12_farther = sprintf('%sI12_farther_%03.0f-%03.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_I12_farther);
    figure;
    imshow(I_behind)
    axis off; axis image
    hold on
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_I12_farther,'png');
end
%}

% Filter by # inliers
t_inlier = 0.5;
is_valid = (seg_inlier>t_inlier)&(seg_lean<p_theta)&(seg_taper<p_taper_max);


%% Additional visualization
%{
% Andrieu Inliers %
%
filepath_andrieu_inlier = sprintf('%sandrieu_inlier_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_inlier,'file');
    I12ieq_inlier = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,seg_inlier,...
        seg_index,axis_a,axis_e,'Inlier [%]',[0 2],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_inlier,'png');
    close all;
end

% Andrieu fill %
%
%filepath_andrieu_inlier = sprintf('%sandrieu_inlier_%03.0f-%02.0f.png',path_png,info_site,info_plot);
filepath_andrieu_fill = sprintf('%sandrieu_fill_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_fill,'file');
    I12ieq_fill = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,seg_fill,...
        seg_index,axis_a,axis_e,'Fill [%]',[0 2],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_fill,'png');
    close all;
end

% Andrieu taper angle
%
%filepath_andrieu_fill = sprintf('%sandrieu_fill_%03.0f-%02.0f.png',path_png,info_site,info_plot);
filepath_andrieu_fill = sprintf('%sandrieu_fill_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_fill,'file');
    I12ieq_taper = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,seg_taper,...
        seg_index,axis_a,axis_e,'Taper [degrees]',[],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_fill,'png');
    close all;
end


% Andrieu Lean angle
%
%filepath_andrieu_lean = sprintf('%sandrieu_lean_%03.0f-%02.0f.png',path_png,info_site,info_plot);
filepath_andrieu_lean = sprintf('%sandrieu_lean_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_lean,'file');
    I12ieq_lean = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,seg_lean,...
        seg_index,axis_a,axis_e,'Lean Angle [degrees]',[],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_lean,'png');
    close all;
end

% Andrieu obj_size
%
%filepath_andrieu_obj_size = sprintf('%sandrieu_obj_size_%03.0f-%02.0f.png',path_png,info_site,info_plot);
filepath_andrieu_obj_size = sprintf('%sandrieu_obj_size_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_obj_size,'file');
    I12ieq_obj_size = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,seg_min_obj_size,...
        seg_index,axis_a,axis_e,'Object Size',[],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_obj_size,'png');
    close all;
end

% Andrieu Diameter
%
%filepath_andrieu_diameter = sprintf('%sandrieu_fill_%03.0f-%02.0f.png',path_png,info_site,info_plot);
filepath_andrieu_diameter = sprintf('%sandrieu_diameter_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_diameter,'file');
    I12ieq_diameter = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,mean(2*seg_r,2),...
        seg_index,axis_a,axis_e,'Diameter',[],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_diameter,'png');
    close all;
end


% Andrieu Farther %
%
%filepath_andrieu_farther = sprintf('%sandrieu_farther_%03.0f-%02.0f.png',path_png,...
%    info_site,info_plot);
filepath_andrieu_farther= sprintf('%sandrieu_farther_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if false;%~exist(filepath_andrieu_farther,'file');
    I12ieq_farther = andrieu_colorLUT_axes4(data_x, data_y,data_z,I12ieq,seg_farther,...
        seg_index,axis_a,axis_e,'Farther [%]',[0 2],...
        seg_col,seg_row,is_valid,false);
    f = getframe(gcf);
    imwrite(f.cdata,filepath_andrieu_farther,'png');
    close all;
end

figure;
hold on;
options_cyl_alpha = 0.5;
for s = 1:n_seg;
    if is_valid(s); color = 'b'; end
    if ~is_valid(s); color = 'r'; end
[Coneh,End1h,End2h] = Cone(seg_z(s,:),...
            seg_y(s,:),...
            [seg_r(s,1) seg_r(s,2)],...
            30,...
            color,1,0);
        set(Coneh,'facealpha',options_cyl_alpha)
        set(End1h,'facealpha',options_cyl_alpha)
        set(End2h,'facealpha',options_cyl_alpha)
end
%}
%% Point cloud segment visualization
%{
% All cylinders
cyl(1,:) = seg_y(:,1);
cyl(2,:) = seg_y(:,2);
cyl(3,:) = seg_y(:,3);
cyl(4,:) = seg_r(:,1);
cyl(5,:) = seg_z(:,1);
cyl(6,:) = seg_z(:,2);
cyl(7,:) = seg_z(:,3);
cyl(8,:) = seg_r(:,2);
%cyl = cyl*1000;
filepath_segments = sprintf('%sseg_all_%03.0f-%02.0f.ply',path_ply,info_site,info_plot);
if ~exist(filepath_segments,'file');
cyl2ply(filepath_segments,cyl);
end
path_local = 'Z:\Desktop\Local\';
filepath_segments = sprintf('%sseg_all_%03.0f-%02.0f.ply',path_local,info_site,info_plot);
if ~exist(filepath_segments,'file');
cyl2ply(filepath_segments,cyl);
end
%}
%{
clear cyl
cyl(1,:) = seg_y(is_valid,1);
cyl(2,:) = seg_y(is_valid,2);
cyl(3,:) = seg_y(is_valid,3);
cyl(4,:) = seg_r(is_valid,1);
cyl(5,:) = seg_z(is_valid,1);
cyl(6,:) = seg_z(is_valid,2);
cyl(7,:) = seg_z(is_valid,3);
cyl(8,:) = seg_r(is_valid,2);
%cyl = cyl*1000;
filepath_segments_valid = sprintf('%sseg_valid_%03.0f-%02.0f.ply',path_ply,info_site,info_plot);
if ~exist(filepath_segments_valid,'file');
cyl2ply(filepath_segments_valid,cyl);
end
path_local = 'Z:\Desktop\Local\';
filepath_segments_valid = sprintf('%sseg_valid_%03.0f-%02.0f.ply',path_local,info_site,info_plot);
if ~exist(filepath_segments_valid,'file');
cyl2ply(filepath_segments_valid,cyl);
end
%}
filepath_tree = sprintf('%s%s%s',path_mat,'tree','.mat');
if exist(filepath_tree,'file');
    load(filepath_tree);
else
    t_rsearch = 0.15;
    tree  = seg2tree( seg_z, seg_y, seg_r, seg_anorm, seg_iter,...
        is_valid,t_rsearch, path_mat );
    save(filepath_tree,'tree');
end

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

if strcmp(info_terminate, 'tree');
    return
end

%% Pixel-based Validation

% Load truth ROI image
path_ROI = sprintf('%sROI%s',path_common,'\');
filepath_ROI = sprintf('%s%s%s',path_mat,'ROI','.mat');
if exist(filepath_ROI,'file');
    load(filepath_ROI);
else
    filepath_ROI_source = sprintf('%sROI_%03.0f-%02.0f.txt',path_ROI,info_site,info_plot);
    ROI = roiRead(filepath_ROI_source);
    save(filepath_ROI,'ROI');
end

filepath_pixel_results = sprintf('%s%s%s',path_mat,'pixel_results','.mat');
filepath_I_ROI_id = sprintf('%s%s%s',path_mat,'I_ROI_id','.mat');
filepath_I_ROI_range = sprintf('%s%s%s',path_mat,'I_ROI_range','.mat');
filepath_I_ROI_is_v = sprintf('%s%s%s',path_mat,'I_ROI_is_v','.mat');
filepath_I_SEG_id = sprintf('%s%s%s',path_mat,'I_SEG_id','.mat');
filepath_I_results_v = sprintf('%s%s%s',path_mat,'I_results_v','.mat');
filepath_I_results_vo2 = sprintf('%s%s%s',path_mat,'I_results_vo2','.mat');
if exist(filepath_pixel_results,'file');
    load(filepath_I_ROI_id);
    load(filepath_I_ROI_range);
    load(filepath_pixel_results);
    load(filepath_I_ROI_is_v);
    load(filepath_I_SEG_id);
    load(filepath_I_results_v);
    load(filepath_I_results_vo2);
elseif ~exist(filepath_pixel_results, 'file');
    t_angular_resolution = 0.25; %deg
    is_valid_tree = true(size(tree_row,1),1);
    tree_inlier = [1:size(tree_row,1)]';
    %[I_ROI_id,I_ROI_range,I_ROI_is_v, I_SEG_id, I_results_v,I_results_vo2, pixel_results] = ROI_analysis(ROI, ...
    %    I12r,I12xy,I12z0,I12e,seg_row,seg_col,seg_iter, seg_inlier,seg_r, is_valid,t_angular_resolution);
    [I_ROI_id,I_ROI_range,I_ROI_is_v, I_SEG_id, I_results_v,I_results_vo2, pixel_results] = ROI_analysis(ROI, ...
        I12r,I12xy,I12z0,I12e,tree_row,tree_col,tree_iter, tree_inlier,seg_r, is_valid_tree,t_angular_resolution);
    save(filepath_pixel_results,'pixel_results');
    save(filepath_I_ROI_id,'I_ROI_id');
    save(filepath_I_ROI_range,'I_ROI_range');
    save(filepath_I_ROI_is_v,'I_ROI_is_v');
    save(filepath_I_SEG_id,'I_SEG_id');
    save(filepath_I_results_v,'I_results_v');
    save(filepath_I_results_vo2,'I_results_vo2');
end

% t_angular_resolution = 0.25; %deg
% [I_ROI_id,I_ROI_range,I_ROI_is_v, I_SEG_id, I_results_v,I_results_vo2, pixel_results] = ROI_analysis(ROI, ...
%     I12r,I12xy,I12z0,I12e,seg_row,seg_col,seg_iter, seg_inlier,seg_r, is_valid,t_angular_resolution);


%filepath_Ipng_results_v = sprintf('%s%s%s',path_png,'Ipng_results_v','.png');
filepath_Ipng_results_v= sprintf('%sIpng_results_v_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if exist(filepath_Ipng_results_v, 'file');
elseif ~exist(filepath_Ipng_results_v,'file');
    cmap = zeros(4,3);
    cmap(1,:) = [ 1 0 0];
    cmap(2,:) = [ 0 0 1];
    cmap(3,:) = [ 0 0 0];
    cmap(4,:) = [ 0 1 0];
    figure; imagesc(I_results_v); axis image; colormap(cmap);
    Ipng_results_r = zeros(size(I_results_v));
    Ipng_results_g = zeros(size(I_results_v));
    Ipng_results_b = zeros(size(I_results_v));
    Ipng_results_r(I_results_v==1) = 255;
    Ipng_results_g(I_results_v==1) = 0;
    Ipng_results_b(I_results_v==1) = 0;
    Ipng_results_r(I_results_v==2) = 0;
    Ipng_results_g(I_results_v==2) = 0;
    Ipng_results_b(I_results_v==2) = 255;
    Ipng_results_r(I_results_v==4) = 0;
    Ipng_results_g(I_results_v==4) = 255;
    Ipng_results_b(I_results_v==4) = 0;
    Ipng_results_v(:,:,1) = Ipng_results_r;
    Ipng_results_v(:,:,2) = Ipng_results_g;
    Ipng_results_v(:,:,3) = Ipng_results_b;
    imwrite(Ipng_results_v,filepath_Ipng_results_v,'png');
end

%filepath_Ipng_results_vo2 = sprintf('%s%s%s',path_png,'Ipng_results_vo2','.png');
filepath_Ipng_results_vo2= sprintf('%sIpng_results_vo2_%03.0f-%02.0f.png',path_latex,info_site,info_plot);
if exist(filepath_Ipng_results_vo2, 'file');
elseif ~exist(filepath_Ipng_results_vo2,'file');
    cmap = zeros(6,3);
    cmap(1,:) = [1 0 0]; %red
    cmap(2,:) = [0 0 1]; %blue
    cmap(3,:) = [0 1 1]; % cyan
    cmap(4,:) = [0 0 0]; %black
    cmap(5,:) = [0 1 0]; %green
    cmap(6,:) = [1 1 0]; % yellow
    figure; imagesc(I_results_vo2); axis image; colormap(cmap);
    Ipng_results_r = zeros(size(I_results_v));
    Ipng_results_g = zeros(size(I_results_v));
    Ipng_results_b = zeros(size(I_results_v));
    Ipng_results_r(I_results_vo2==1) = 255;
    Ipng_results_g(I_results_vo2==1) = 0;
    Ipng_results_b(I_results_vo2==1) = 0;
    Ipng_results_r(I_results_vo2==2) = 0;
    Ipng_results_g(I_results_vo2==2) = 0;
    Ipng_results_b(I_results_vo2==2) = 255;
    Ipng_results_r(I_results_vo2==3) = 0;
    Ipng_results_g(I_results_vo2==3) = 255;
    Ipng_results_b(I_results_vo2==3) = 255;
    Ipng_results_r(I_results_vo2==5) = 0;
    Ipng_results_g(I_results_vo2==5) = 255;
    Ipng_results_b(I_results_vo2==5) = 0;
    Ipng_results_r(I_results_vo2==6) = 255;
    Ipng_results_g(I_results_vo2==6) = 255;
    Ipng_results_b(I_results_vo2==6) = 0;
    Ipng_results_vo2(:,:,1) = Ipng_results_r;
    Ipng_results_vo2(:,:,2) = Ipng_results_g;
    Ipng_results_vo2(:,:,3) = Ipng_results_b;
    imwrite(Ipng_results_vo2,filepath_Ipng_results_vo2,'png');
end

%% Forestry Evaluation

% Parse Field Data
filepath_NEON_DBH_10over_csv = sprintf('%sField_Data\\NEON_DBH_10over.csv',path_common);
filepath_NEON_DBH_10over = sprintf('%sNEON_DBH_10over.mat',path_mat);
if exist(filepath_NEON_DBH_10over,'file');
    load(filepath_NEON_DBH_10over);
elseif ~exist(filepath_NEON_DBH_10over);
    NEON_DBH_10over = parse_NEON_DBH_10over(filepath_NEON_DBH_10over_csv,info_site);
    save(filepath_NEON_DBH_10over,'NEON_DBH_10over');
end

filepath_NEON_DBH_10under_csv = sprintf('%sField_Data\\NEON_DBH_10under.csv',path_common);
filepath_NEON_DBH_10under = sprintf('%sNEON_DBH_10under.mat',path_mat);
if exist(filepath_NEON_DBH_10under,'file');
    load(filepath_NEON_DBH_10under);
elseif ~exist(filepath_NEON_DBH_10under);
    NEON_DBH_10under = parse_NEON_DBH_10under(filepath_NEON_DBH_10under_csv,info_site);
    save(filepath_NEON_DBH_10under,'NEON_DBH_10under');
end

% Get ROI parameter information;
%[roi_z, roi_y, roi_r] = get_roi_param(I12x, I12y, I12z, I_ROI_id);

filepath_n_stems = sprintf('%s%s%s',path_mat,'n_stems','.mat');
filepath_ba = sprintf('%s%s%s',path_mat,'ba','.mat');
filepath_match_SEG = sprintf('%s%s%s',path_mat,'match_SEG','.mat');
filepath_match_TREE = sprintf('%s%s%s',path_mat,'match_TREE','.mat');

if false;%exist(filepath_n_stems,'file');
    load(filepath_n_stems);
    load(filepath_ba);
    load(filepath_match_SEG);
    load(filepath_match_TREE);
else
    [ n_stems, ba, atdbh, match_SEG, match_TREE] ...
        = forestry_evaluation_tree( seg_z, seg_y,seg_r,tree, I12z0,I12x, I12y, I12z,I12r,I_SEG_id,I_ROI_id, I_ROI_range,...
        NEON_DBH_10over, NEON_DBH_10under, path_mat, info_site, info_plot);
    %   [ n_stems, ba, atdbh, match] ...
    %        = forestry_evaluation_tree_basic( seg_z, seg_y, seg_r, tree,...
    %        NEON_DBH_10over, NEON_DBH_10under);
    save(filepath_n_stems,'n_stems');
    save(filepath_ba,'ba');
    save(filepath_match_SEG,'match_SEG');
    save(filepath_match_TREE,'match_TREE');
end

% Plot DBH maps
%filepath_dbhmap = sprintf('%s%s%s',path_eps,'dbhmap','.eps');
filepath_dbhmap= sprintf('%sdbhmap_%03.0f-%02.0f.eps',path_latex,info_site,info_plot);
if ~exist(filepath_dbhmap,'file');
    figure('position', [587 95 1026 832]);
    hold on
    for s = 1:numel(atdbh.ROI_x);
        h = filledCircle([atdbh.ROI_x(s) atdbh.ROI_y(s)],atdbh.ROI_r(s),1000,'g');
        set(h,'FaceAlpha',.5);
    end
    % for s = 1:numel(atdbh.SEG_x);
    %     h = filledCircle([atdbh.SEG_x(s) atdbh.SEG_y(s)],atdbh.SEG_r(s),1000,'r');
    %     set(h,'FaceAlpha',.2);
    % end
    for s = 1:numel(atdbh.NEON_x)
        h = filledCircle([atdbh.NEON_x(s) atdbh.NEON_y(s)],atdbh.NEON_r(s),1000,'b');
        set(h,'FaceAlpha',.5);
    end
    for s = 1:numel(atdbh.TREE_x)
        h = filledCircle([atdbh.TREE_x(s) atdbh.TREE_y(s)],atdbh.TREE_r(s),1000,'r');
        set(h,'FaceAlpha',.5);
    end
    h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
    axis equal;
    axis([-16 16 -16 16]);
    xlabel('Easting');
    ylabel('Northing');
    grid on
    print(gcf,'-depsc','-opengl',filepath_dbhmap)
end

figure;
hold on
n_step = 100;
theta = linspace(0,2*pi,n_step);
% h_leg = scatter(0,0,10,[139,255,124]/255,'filled');
% set(h_leg, 'visible','off');
% legend_str{1} = 'Inventory';
% h_leg = scatter(0,0,10,[245,117,127]/255,'filled');
% set(h_leg, 'visible','off');
% legend_str{2} = 'Lidar';
% h_leg = scatter(0,0,10,[140,125,254]/255,'filled');
% set(h_leg, 'visible','off');
% legend_str{3} = 'ROI';
% h_leg = plot(0,0,'sk');
% set(h_leg, 'visible','off');
% legend_str{4} = 'Plot boundary';
for s = 1:numel(atdbh.ROI_x);
    xcirc = atdbh.ROI_x(s) + atdbh.ROI_r(s)*cos(theta);
    ycirc = atdbh.ROI_y(s) + atdbh.ROI_r(s)*sin(theta);
    fill(xcirc, ycirc, 'b', 'FaceAlpha', 0.5)
end
for s = 1:numel(atdbh.NEON_x);
    xcirc = atdbh.NEON_x(s) + atdbh.NEON_r(s)*cos(theta);
    ycirc = atdbh.NEON_y(s) + atdbh.NEON_r(s)*sin(theta);
    fill(xcirc, ycirc, 'g', 'FaceAlpha', 0.5)
end
for s = 1:numel(atdbh.TREE_x);
    xcirc = atdbh.TREE_x(s) + atdbh.TREE_r(s)*cos(theta);
    ycirc = atdbh.TREE_y(s) + atdbh.TREE_r(s)*sin(theta);
    fill(xcirc, ycirc, 'r', 'FaceAlpha', 0.5)
end
axis equal;
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
%hleg = legend(legend_str);
%set(hleg, 'position', [.8 .4 .2 .2]);
axis([-16 16 -16 16]);
xlabel('Easting from plot center [m]');
ylabel('Northing from plot center [m]');
grid on
path_latex = 'Z:\Desktop\Paper\';
filename = sprintf('%s%s', path_latex,'dbhmap.tex');
matlab2tikz(filename);

%%
%{
% Plot all tree candidate points
if options_show_fig
    h = figure;
elseif any([options_save_eps, options_save_png, options_save_fig]);
    h = figure('visible','off');
end
if any([options_show_fig, options_save_eps, options_save_png, options_save_fig]);
    hold on
    scatter3(data_x(index_inliers),data_y(index_inliers),data_z(index_inliers),...
        10,'r','filled')
    scatter3(data_x(index_outliers),data_y(index_outliers),data_z(index_outliers),...
        10,'b','filled')
    %             [Coneh,End1h,End2h] = Cone(z1,y1,[r1 r1],30,'g',1,0);
    %             set(Coneh,'facealpha',.2)
    %             set(End1h,'facealpha',.2)
    %             set(End2h,'facealpha',.2)
    %axis(options_plotlim);
    %view(options.view);
    view(0,0);
    axis equal
    grid on; hold on;
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    % title(info_case);
end
if options_save_eps
    if ~exist(path_eps, 'dir');
        mkdir(path_eps);
    end
    name = sprintf('Sub_%03.0f',iteration);
    filename = sprintf('%s%s.eps',path_eps,name);
    print(h, '-depsc', filename);
end
if options_save_png
    if ~exist(path_png, 'dir');
        mkdir(path_png);
    end
    name = sprintf('Sub_%03.0f',iteration);
    filename = sprintf('%s%s.png',path_png,name);
    print(h, '-dpng', filename);
end
if options_save_fig
    if ~exist(path_fig, 'dir');
        mkdir(path_fig);
    end
    name = sprintf('Sub_%03.0f',iteration);
    filename = sprintf('%s%s.fig',path_fig,name);
    saveas(h,filename,'fig');
end
if  any([options_save_eps, options_save_png, options_save_fig])
    delete(h);
end

if options_write_ply
    filepath_ply_tree_cand_pts  = sprintf('%s%s%03.0f-%02.0f.ply',path_ply,'tree_cand_pts_',info_site, info_plot);
    filepath_ply_tree_cand_pts_old  = sprintf('%s%s.ply',path_ply,'tree_cand_pts');
    if exist(filepath_ply_tree_cand_pts_old,'file');
        command = sprintf('move %s %s',filepath_ply_tree_cand_pts_old, filepath_ply_tree_cand_pts);
        [~,~] = system(command);
    end
    if ~exist(filepath_ply_tree_cand_pts,'file');
        color = zeros(numel(data_x),3);
        color(index_inliers,:) = repmat([255 0 0],[numel(index_inliers),1]);
        color(index_outliers,:) = repmat([0 0 255],[numel(index_outliers),1]);
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply_tree_cand_pts);
        end
        write2ply(filepath_ply_tree_cand_pts,[data_x data_y data_z], color);
    end
end

if options_write_ply
    filepath_ply_tree_cand_cyl  = sprintf('%s%s_%03.0f-%02.0f.ply',path_ply,'tree_cand_cyl', info_site,info_plot);
    filepath_ply_tree_cand_cyl_old  = sprintf('%s%s.ply',path_ply,'tree_cand_cyl');
    if exist(filepath_ply_tree_cand_cyl_old,'file');
        command = sprintf('move %s %s',filepath_ply_tree_cand_cyl_old, filepath_ply_tree_cand_cyl);
        [~,~] = system(command);
    end
    if ~exist(filepath_ply_tree_cand_cyl,'file')
        cyl = [seg_z'; seg_r'; seg_y'; seg_r'];
        h = sqrt((cyl(1,:)-cyl(5,:)).^2 + (cyl(2,:)-cyl(6,:)).^2 + (cyl(3,:)-cyl(7,:)).^2);
        degenerate = (h==0);
        color = repmat([255; 0; 0;], [1,n_seg]);
        if options_verbose_progress;
            fprintf('Writing %s \n',filepath_ply_tree_cand_cyl);
        end
        cyl2ply(filepath_ply_tree_cand_cyl, cyl(:,~degenerate), color(:,~degenerate));
    end
end

if options_verbose_output;
    fprintf('\n');
    if options_write_ply;
        fprintf('Output PLY: %s \n', filepath_ply_tree_cand_pts);
        fprintf('Output PLY: %s \n', filepath_ply_tree_cand_cyl);
    end
end
%}
%}

%}
