function [ data_i, data_z, data_n, data_r, data_a, data_e ] = read_lidar_newlasprocessor( filename )
% READ_LIDAR Reads lidar txt files to produce arrays
%   filename should include full path and .txt extension
% 
%   Example:
%   filepath_source = 'D:\Users\djk2312\Documents\Harvard\2012-08-Harvard-Lidar\SICK_NONE_2012-08-15_111923.txt'
%   [data_i, data_z, data_n, data_r, data_a, data_e]  = read_lidar2(filepath_source);
%   % data_x, data_y, and data_z can then be calculated from the angular
%   % measurements 
% 
% (c) David Kelbe, Rochester Institute of Technology



xyzPath = [  filename(1:end-4) '_xyz.txt' ];
polarPath = [  filename(1:end-4) '_polar.txt' ];

fid = fopen(xyzPath);
rawdata = textscan(fid, '%f%f%f%f%f');
data_i = rawdata{4};
%data_x = rawdata{1}/1000; 
%data_y = rawdata{2}/1000; 
data_z = rawdata{3}/1000;
data_n = rawdata{5};

%data_xy = sqrt((data_x).^2+(data_y).^2);

fclose all;
fid = fopen(polarPath);
rawdata = textscan(fid, '%f%f%f%f%f');
data_r = rawdata{1}/1000;
data_a = rawdata{2};
data_e = rawdata{3};
end

