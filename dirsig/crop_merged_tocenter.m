%% Crop merged point clouds to site 

path_exp = 'D:\Users\djk2312\Documents\Harvard\reg\';
info_site = 31;
path_ply = sprintf('%s%03.0f%s',path_exp,info_site, '\13\ply\paper\');
path_crop = sprintf('%s%03.0f%s',path_exp,info_site, '\13\ply\crop\');
if ~exist(path_crop,'dir');
    mkdir(path_crop);
end
cd(path_ply);
D = dir('*Ptx*');
D = remove_hiddenfiles(D);

xmin = -20;
xmax = 20;
ymin = -20;
ymax = 20;

n_plot = numel(D);
for p = 1:n_plot;
    filepath_ply = sprintf('%s%s',path_ply, D{p});        
    fid = fopen(filepath_ply);
    data = textscan(fid, '%f%f%f%u%u%u', 'headerlines',10 );
    x = data{1};
    y = data{2};
    isvalid = (x>xmin)&(x<xmax)&(y>ymin)&(y<ymax);
    data{1} = data{1}(isvalid);
    data{2} = data{2}(isvalid);
    data{3} = data{3}(isvalid);
    data{4} = data{4}(isvalid);
    data{5} = data{5}(isvalid);
    data{6} = data{6}(isvalid);
    n_pts = numel(data{1});
    
    
    filepath_out = sprintf('%s%s', path_crop, D{p});
    pointcloud = [data{1} data{2} data{3}];
    color = [data{4} data{5} data{6}];
    
    write2ply(filepath_out, pointcloud, color);
    fclose all
end