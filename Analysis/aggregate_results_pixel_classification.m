% Aggregate Results
%% Initialization
% Set Internal Matlab Parameters
set(0,'DefaultFigureRenderer','OpenGL')
set(0,'defaultfigureposition',[966  322  557 472]);
set(0,'DefaultFigureColor', 'White');
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontName', 'Calibri')
set(0,'DefaultAxesFontSize', 14)

% Info
info_exp = 'Harvard';
info_suffix = '03-01';
info_slash = '\';

site_arr = 1:31;
plot_arr = 31;
n_site = numel(site_arr);
n_plot = numel(plot_arr);

%pixel_results_all = cell(n_site*n_plot,1);
site_all = zeros(n_site*n_plot,1);
plot_all = zeros(n_site*n_plot,1);
pixel_results_all = cell(n_site*n_plot,1);
i = 1;

% Load data
for info_site = 1:31;
    info_plot = 13;
    
    
    % Make directories
    path_common = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, 'Common',info_slash);
    path_top = sprintf('%s%s%s%s%s%03.0f%s%02.0f%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, info_suffix,info_slash,info_site,info_slash,info_plot,info_slash);
    path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
    path_fig = sprintf('%s%s%s',path_top,'fig',info_slash);
    path_eps = sprintf('%s%s%s',path_top,'eps',info_slash);
    path_ply = sprintf('%s%s%s',path_top,'ply',info_slash);
    path_png = sprintf('%s%s%s',path_top,'png',info_slash);
    
    if ~exist(path_mat,'dir');
        i = i + 1;
        continue
    end
    
    site_all(i) = info_site;
    plot_all(i) = info_plot;
    
    filepath_pixel_results = sprintf('%s%s%s',path_mat,'pixel_results','.mat');
    %filepath_I_ROI_id = sprintf('%s%s%s',path_mat,'I_ROI_id','.mat');
    %filepath_I_ROI_range = sprintf('%s%s%s',path_mat,'I_ROI_range','.mat');
    %filepath_I_ROI_is_v = sprintf('%s%s%s',path_mat,'I_ROI_is_v','.mat');
    %filepath_I_SEG_id = sprintf('%s%s%s',path_mat,'I_SEG_id','.mat');
    %filepath_I_results_v = sprintf('%s%s%s',path_mat,'I_results_v','.mat');
    %filepath_I_results_vo2 = sprintf('%s%s%s',path_mat,'I_results_vo2','.mat');
    if exist(filepath_pixel_results,'file');
        load(filepath_pixel_results);
        pixel_results_all{i} = pixel_results;
    end
     
    %{
    filepath_Ipng_results_v = sprintf('%s%s%s',path_png,'Ipng_results_v','.png');
    filepath_Ipng_results_vo2 = sprintf('%s%s%s',path_png,'Ipng_results_vo2','.png');
    
    filepath_NEON_DBH_10over_csv = sprintf('%sField_Data\\NEON_DBH_10over.csv',path_common);
    filepath_NEON_DBH_10over = sprintf('%sNEON_DBH_10over.mat',path_mat);
    
    filepath_NEON_DBH_10under_csv = sprintf('%sField_Data\\NEON_DBH_10under.csv',path_common);
    filepath_NEON_DBH_10under = sprintf('%sNEON_DBH_10under.mat',path_mat);
    
    filepath_n_stems = sprintf('%s%s%s',path_mat,'n_stems','.mat');
    filepath_ba = sprintf('%s%s%s',path_mat,'ba','.mat');
    filepath_match = sprintf('%s%s%s',path_mat,'match','.mat');
    
    filepath_dbhmap = sprintf('%s%s%s',path_eps,'dbhmap','.eps');
    %}
    
    i = i + 1;
end

is_valid = ~cellfun(@isempty,pixel_results_all);
site_all = site_all(is_valid);
plot_all = plot_all(is_valid);
pixel_results_all = pixel_results_all(is_valid);

fig_classification_by_site_range( site_all, plot_all, pixel_results_all, 'Range [m]' )
fig_classification_by_site_z0( site_all, plot_all, pixel_results_all, 'Height above ground [m]' )
fig_classification_by_site_e( site_all, plot_all, pixel_results_all, 'Elevation angle [deg]' )
fig_classification_by_site_xy( site_all, plot_all, pixel_results_all, ' xy-Range [m]' )
fig_classification_all_r( site_all, plot_all, pixel_results_all, 'Range [m]');
fig_classification_all_z0( site_all, plot_all, pixel_results_all, 'Height above ground [m]');
fig_classification_all_e( site_all, plot_all, pixel_results_all, 'Elevation angle [deg]');
fig_classification_all_xy( site_all, plot_all, pixel_results_all, 'xy-Range [m]');

foo = 1;


