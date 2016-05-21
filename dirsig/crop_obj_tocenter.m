%% Crop martin obj to site 
% maybe still in process of being written 
path_obj = 'D:\Users\djk2312\Documents\Harvard\dirsig\martin_may26\';
cd(path_obj);
D = dir('*.obj');
D = remove_hiddenfiles(D);

xmin = -20;
xmax = 20;
ymin = -20;
ymax = 20;

n_plot = numel(D);
for p = 1:n_plot;
    filepath_obj = sprintf('%s%s',path_obj, D{p});
    fid = fopen(filepath_obj);
    data = textscan(fid, '%s%f%f%f', 'headerlines', 2);
    fclose all
    xmean = mean(data{2});
    ymean = mean(data{3});
    
    if xmean > xmin && xmean < xmax && ymean > ymin && ymean < ymax;
        fprintf('Tree %s: %3.2f, %3.2f\n', D{p}, xmean, ymean)
    else
        filepath_move = sprintf('%s%s%s',path_obj, 'outside\',D{p});
        command = sprintf('mv %s %s',filepath_obj, filepath_move);
        [a,b] = system(command);
    end   
end