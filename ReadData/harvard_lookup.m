function [ filename_lidar ] = harvard_lookup( site,plot )
%HARVARD_LOOKUP Get LiDAR base filename (excl. extension) from site & plot
%   Reads 2012-08-Harvard.csv to get LiDAR filename associated with give
%   site and plot number. 
%
%   Example: filename_lidar = harvard_lookup(31,13)
%   ans = 'SICK_NONE_2012-08-15_111819'
%
%   (c) David Kelbe, Rochester Institute of Technology 


%filename_csv = '/dirs/wasp/ground_lidar_collects/2012-08-Harvard/2012-08-Harvard.csv';
filename_csv = 'D:\Users\djk2312\Documents\thisshouldbecyclone\2012-08-Harvard.csv';
fid = fopen(filename_csv);
C = textscan(fid, '%u8%u8%s%s%s%u16%u16%u16%u16%u16%u16%u16%u16%u16',...
    'delimiter', ',',...
    'headerlines',8);
% 1 Site
% 2 Plot 
% 3 Lidar Filename
% 4 Date 
% 5 Lidar Orientation
i = find((C{1}==site)&(C{2}==plot));

if numel(i)>0
    filename_lidar = char(C{3}(i));
else
    filename_lidar = '';
end

end

