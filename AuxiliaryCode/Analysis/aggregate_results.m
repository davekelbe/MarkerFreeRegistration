% Aggregate Results
%% Initialization
% Set Internal Matlab Parameters
set(0,'DefaultFigureRenderer','painters')
set(0,'defaultfigureposition',[966  322  557 472]);
set(0,'DefaultFigureColor', 'White');
%set(0,'DefaultAxesFontName', 'Times New Roman')
%set(0,'DefaultAxesFontName', 'Calibri')
set(0,'DefaultAxesFontSize', 14)
filepath_latex = 'Z:\Desktop\Paper\Plots\';
% Info
info_exp = 'Harvard';
info_suffix = 'stem';
info_slash = '\';

site_arr = 1:31;
plot_arr = 31;
n_site = numel(site_arr);
n_plot = numel(plot_arr);

%pixel_results_all = cell(n_site*n_plot,1);
site_all = zeros(n_site*n_plot,1);
plot_all = zeros(n_site*n_plot,1);
pixel_results_all = cell(n_site*n_plot,1);
n_stems_all = cell(n_site*n_plot,1);
ba_all = cell(n_site*n_plot,1);
match_TREE_all = cell(n_site*n_plot,1);
match_SEG_all = cell(n_site*n_plot,1);

i = 1;

% Load data
for info_site = 1:31;
    info_plot = 13;
    if info_site==9;
        continue
    end
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
    %{
    %filepath_I_ROI_id = sprintf('%s%s%s',path_mat,'I_ROI_id','.mat');
    %filepath_I_ROI_range = sprintf('%s%s%s',path_mat,'I_ROI_range','.mat');
    %filepath_I_ROI_is_v = sprintf('%s%s%s',path_mat,'I_ROI_is_v','.mat');
    %filepath_I_SEG_id = sprintf('%s%s%s',path_mat,'I_SEG_id','.mat');
    %filepath_I_results_v = sprintf('%s%s%s',path_mat,'I_results_v','.mat');
    %filepath_I_results_vo2 = sprintf('%s%s%s',path_mat,'I_results_vo2','.mat');
    %}
    if exist(filepath_pixel_results,'file');
        load(filepath_pixel_results);
        pixel_results_all{i} = pixel_results;
    end
    
    %{
    filepath_Ipng_results_v = sprintf('%s%s%s',path_png,'Ipng_results_v','.png');
    filepath_Ipng_results_vo2 = sprintf('%s%s%s',path_png,'Ipng_results_vo2','.png');
    %}
    
    filepath_n_stems = sprintf('%s%s%s',path_mat,'n_stems','.mat');
    if exist(filepath_n_stems,'file');
        load(filepath_n_stems);
        n_stems_all{i} = n_stems;
    end
    
    filepath_ba = sprintf('%s%s%s',path_mat,'ba','.mat');
    if exist(filepath_ba,'file');
        load(filepath_ba);
        ba_all{i} = ba;
    end
    
    filepath_match_TREE = sprintf('%s%s%s',path_mat,'match_TREE','.mat');
    if exist(filepath_match_TREE,'file');
        load(filepath_match_TREE);
        match_TREE_all{i} = match_TREE;
    end
    
    filepath_match_SEG = sprintf('%s%s%s',path_mat,'match_SEG','.mat');
    if exist(filepath_match_SEG,'file');
        load(filepath_match_SEG);
        match_SEG_all{i} = match_SEG;
    end
    
    filepath_dbhmap = sprintf('%s%s%s',path_eps,'dbhmap','.eps');
    
    
    i = i + 1;
end
clear n_stems ba match 
is_valid = ~cellfun(@isempty,pixel_results_all);
site_all = site_all(is_valid);
plot_all = plot_all(is_valid);
pixel_results_all = pixel_results_all(is_valid);
n_stems_all = n_stems_all(is_valid);
ba_all = ba_all(is_valid);
match_TREE_all = match_TREE_all(is_valid);
match_SEG_all = match_SEG_all(is_valid);


%% Pixel classification
%
fig_classification_by_site_all( site_all, plot_all, pixel_results_all )
table_classification_by_site_all( site_all, plot_all, pixel_results_all);
%print(gcf,'-depsc2','-painters',sprintf('%s%s.eps',filepath_latex,'classification_by_site_all'));

fig_classification_by_site_range( site_all, plot_all, pixel_results_all, 'Range [m]' )
%fig_classification_by_site_z0( site_all, plot_all, pixel_results_all, 'Height above ground [m]' )
%fig_classification_by_site_e( site_all, plot_all, pixel_results_all, 'Elevation angle [deg]' )
%fig_classification_by_site_xy( site_all, plot_all, pixel_results_all, ' xy-Range [m]' )
%}
%
%
fig_classification_all_r( site_all, plot_all, pixel_results_all, 'Range [m]');
%{
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'classification_all_r'));
fig_classification_all_z0( site_all, plot_all, pixel_results_all, 'Height above ground [m]');
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'classification_all_z0'));
fig_classification_all_e( site_all, plot_all, pixel_results_all, 'Elevation angle [deg]');
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'classification_all_e'));
fig_classification_all_xy( site_all, plot_all, pixel_results_all, 'xy-Range [m]');
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'classification_all_xy'));
%}
%% N stems
%
fig_n_stems_tree_roi(site_all,plot_all,n_stems_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'n_stems_tree_roi'));
fig_n_stems_tree_neon(site_all,plot_all,n_stems_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'n_stems_tree_neon'));

% fig_n_stems_roi(site_all,plot_all,n_stems_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'n_stems_seg_roi'));
%fig_n_stems_seg_neon(site_all,plot_all,n_stems_all);
%saveas(gcf,sprintf('%s%s.eps',filepath_latex,'n_stems_seg_neon'),'psc2')
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'n_stems_seg_neon'));
%}
%% Basal area
fig_ba_tree_roi(site_all,plot_all,ba_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'ba_tree_roi'));
fig_ba_tree_neon(site_all,plot_all,ba_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'ba_tree_neon'));
foo = 1; 
%fig_ba_seg_roi(site_all,plot_all,ba_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'ba_seg_roi'));
%fig_ba_seg_neon(site_all,plot_all,ba_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'ba_seg_neon'));
%}
%% Range
% 
fig_range_tree_neon(site_all, plot_all,match_TREE_all)
%%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'tree_range'));
%}
%% Diameter
%

fig_match_tree_neon_color_xy(site_all, plot_all,match_TREE_all)
fig_match_tree_neon_color_xy_witherror(site_all, plot_all,match_TREE_all)

%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'DBH_vs_xy'));

%{
%fig_match_seg_neon(site_all, plot_all,match_SEG_all)
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'match_SEG_neon'));
%fig_match_tree_neon(site_all, plot_all,match_TREE_all)
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'match_TREE_neon'));
%fig_match_seg_roi(site_all, plot_all,match_SEG_all)
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'match_SEG_roi'));
%fig_match_tree_roi(site_all, plot_all,match_TREE_all)
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'match_TREE_roi'));
%fig_match_tree_neon_color_dbh(site_all, plot_all,match_TREE_all)
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'DBH_vs_DBH'));
%}
%%

%
fig_dbh_seg_neon_xy(site_all,plot_all,match_SEG_all);
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_SEG_neon_xy'));
fig_dbh_seg_neon_dbh(site_all,plot_all,match_SEG_all);
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_SEG_neon_dbh'));
%fig_dbh_roi_xy(site_all,plot_all,match_SEG_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_SEG_roi_xy'));
%fig_dbh_roi_dbh(site_all,plot_all,match_SEG_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_SEG_roi_dbh'));

fig_dbh_tree_neon_xy(site_all,plot_all,match_TREE_all);
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_TREE_neon_xy'));
fig_dbh_tree_neon_dbh(site_all,plot_all,match_TREE_all);
print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_TREE_neon_dbh'));
%fig_dbh_roi_xy(site_all,plot_all,match_TREE_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_TREE_roi_xy'));
%fig_dbh_roi_dbh(site_all,plot_all,match_TREE_all);
%print(gcf,'-depsc','-painters',sprintf('%s%s.eps',filepath_latex,'dbh_TREE_roi_dbh'));
%}
foo = 1;


