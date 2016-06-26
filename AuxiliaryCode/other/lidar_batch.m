
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

for s = 306:n_scans;
    site = C{1}(s);
    plot = C{2}(s);
    fprintf('\nSite %d Plot %d\n', site, plot)
    lidar(site,plot,'tree', 'harvard');
end





