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
tic;
site_unique = unique(C{1});
n_site = numel(site_unique);
aux.n_trial = 1;
for s = 1:n_site;
    info_site = site_unique(s);
    is_site = (C{1}==info_site);
    plot = C{2}(is_site);
    n_plot2 = numel(plot);
    n_plot = 2;
    info_valid_plot = cell(n_plot,1);
    %pix = randperm(numel(plot),n_plot);
    %pvalid = plot(pix);
    for pi = 1:n_plot2;
        %p = pvalid(pi);
        p = plot(pi);
        fprintf('\nSite %d - Plot %d \n', info_site, p);
        info_valid_plot{1} = sprintf('%02.0f', p);
        info_valid_plot{2} = sprintf('%02.0fm', p);
        aux.info_valid_plot = info_valid_plot;
        aux.info_site = info_site;
        %kelbe_analysis_pose_fun( aux );
        %kelbe_analysis_subset_fun( aux );
        %kelbe_analysis_noise_fun( aux );
        kelbe_analysis_rmse_fun( aux );
        %kelbe_analysis_dbh_fun( aux );
    end
end
time = toc; 
foo = 1;



foo = 1;
