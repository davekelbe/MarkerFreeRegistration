function [ datadem, qx, qy, qz ] = digital_elevation_model( ...
    data_x,data_y,data_z,data_xy, data_a, data_e,...
    t_maxr, t_azimuth_bins, ...
    t_angle_deviation, t_grid_size, ...
    write_ply, path_top,...
    a_blocked, xy_blocked, info_site, info_plot)
%DIGITAL_ELEVATION_MODEL Generates the DEM for a lidar scan "filename"
%   filename = full path of lidar file without the extension, e.g.
%              SICK_graciesforest_2012-01-27_172248
%   output is a gridded DEM (Matlab variable) at 0.5m resolution
%         and a list of ground heights for each point in the input file,
%         This file is a matlab variable called 'datadem'
%
% Questions contact Dave Kelbe, djk2312@cis.rit.edu
% 4: Changed NaN problem

verbose =false; %options.verbose; % Set to 1 if you want to see images
saveim =false ;%options.saveim;  % Set to 1 if you want to save images
savedata= false;%options.savedata; % Set to 1 if you want to save DEM to computer
pngPath = 'foo';%[options.path,'fig/dem/'];
options_name = 'foo';
% Change this path to 
                                                 % save Images in a diff
                                                 % directory
%datadirout = [options.path,'out/dem/']; % Change this path to 
                                                 % directory
%name = graciesplot(plot1);% Change this to refert back to old lookup
%system
% Define a new variable 'name' = 'fullpath filename of lidar file (no
% extensions)

% Visualization point subsampling
nplotp = 50000; 
ixp = randi([1 size(data_x,1)],nplotp,1);
%ixp = 1*data_x;

%if  verbose; lidarplot2(data,ixp); end;
%if saveim;  saveas( gcf(), [ pngPath 'scatter' num2str( options.name ),'.png' ], 'png' ); end

%%
[minindex, minindex_bad] = radial_subsample( data_xy, data_a, data_e, data_z,...
    t_azimuth_bins,...
    a_blocked, xy_blocked);

temp = 1:1:numel(minindex);
minindexl = temp(logical(minindex));

%if verbose; lidarplot(data,minindex); end;
%if saveim;  saveas( gcf(), [ pngPath 'minindex' num2str( options.name ),'.png' ], 'png' ); end

 %%
[~,ix_unique,~] = unique([data_x(minindexl),data_y(minindexl)],'rows');
is_unique = false(size(minindexl));
is_unique(ix_unique) = true;
minindexl = minindexl(is_unique);
dt = DelaunayTri(data_x(minindexl),data_y(minindexl));
minindex = false(size(minindex));
minindex(minindexl)= true;

if verbose;
figure;
hidden on
trimesh(dt.Triangulation,data_x(minindexl),...
    data_y(minindexl),...
    data_z(minindexl))
xlabel('x, [m]'),ylabel('y, [m]'),zlabel('z, [m]');   
view(0,90);
axis equal 
end;
if saveim;  saveas( gcf(), [ pngPath 'delaunay' num2str( options.name ),'.png' ], 'png' ); end
   
if write_ply;
    plycolor = vec2cmap(data_z(minindexl),'jet');
    filepath_ply_delaunay = sprintf('%sply\\delaunay_%03.0f-%02.0f.ply',path_top,info_site, info_plot);
    write2plyfaces( filepath_ply_delaunay,  ...
        [data_x(minindexl) data_y(minindexl) data_z(minindexl)], ...
        plycolor, ...
        dt.Triangulation )
end
    
%%


[tri_x, tri_y, tri_z] = point_removal(data_x, data_y, data_z, minindex, dt, t_angle_deviation, pngPath, options_name, 0, saveim);

if verbose;
figure;
hidden on
trimesh(dt.Triangulation,tri_x,...
    tri_y,...
    tri_z);
xlabel('x, [m]'),ylabel('y, [m]'),zlabel('z, [m]');   
view(0,90);
end;
if saveim;  saveas( gcf(), [ pngPath 'delaunaynn' num2str( options.name ),'.png' ], 'png' ); end

plycolor = vec2cmap(tri_z,'summer');
if write_ply;
    filepath_ply_delaunay = sprintf('%sply\\delaunay_smoothed_%03.0f-%03.0fdemcmap.ply',path_top, info_site, info_plot);
    write2plyfaces( filepath_ply_delaunay,  ...
        [tri_x tri_y tri_z], ...
        plycolor, ...
        dt.Triangulation )
end
%%
%tx = min(data_x(nnindexl)):0.5:max(data_x(nnindexl));
%ty = min(data_y(nnindexl)):0.5:max(data_y(nnindexl));

tx =  -t_maxr:t_grid_size:t_maxr; 
ty = -t_maxr:t_grid_size:t_maxr;
F = TriScatteredInterp(tri_x,tri_y,tri_z);
[qx,qy] = meshgrid(tx,ty);
qz = F(qx,qy);

qz2 = inpaint_nans(qz,1);
if verbose
figure('color', [1 1 1 ]);
mesh(qx,qy,qz);
    axis('equal');
   % xlabel('x, [meters]', 'fontsize', 14);
   % ylabel('y, [meters]', 'fontsize', 14);
   % zlabel('z, [meters]', 'fontsize', 14);
    colormap('jet');
    view(45,30);
    grid off;
   % title('DEM', 'fontsize', 14)
   % ylabel('y range, [meters]');
   % colorbar
end   
    
if saveim; saveas( gcf(), [ pngPath 'fullplot' num2str( options.name ),'.png' ], 'png' ); end;
if verbose;
hold on;
axis auto;
scatter3(data_x(ixp),data_y(ixp),data_z(ixp),3,[0.25, 0.25 0.25],'filled')
axis equal;
hold off
end
%if saveim; saveas( gcf(), [ pngPath 'dempointsplot' num2str( options.name ),'.png' ], 'png' ); end;

F2 = TriScatteredInterp(qx(:), qy(:), qz2(:));
datadem = F2(data_x,data_y);


if savedata;
    filename = [ datadirout 'griddem' num2str( options.name )];
    save(filename, 'qx','qy','qz');
    filename = [ datadirout 'pointdem' num2str( options.name )];
    save(filename, 'datadem');
end
end

