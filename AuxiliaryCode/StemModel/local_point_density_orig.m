function [ npts_2d, spacing_u, spacing_v ] = local_point_density_orig( all_u,all_v,all_w, M1, info_point_spacing )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Centers in uvw space 
%
npts_2d = [];


spacing_u = (min(all_u):5*info_point_spacing:max(all_u)+info_point_spacing)';
spacing_v = (min(all_v):5*info_point_spacing:max(all_v)+info_point_spacing)';
spacing_w = (min(all_w):5*info_point_spacing:max(all_w)+info_point_spacing)';

nu = numel(spacing_u);
nv = numel(spacing_v);
nw = numel(spacing_w);

if nu<2 || nv<2 || nw<2;
    return
end

[mesh_u, mesh_v, mesh_w] = meshgrid(spacing_u, spacing_v, spacing_w);
mesh1d_u = reshape(mesh_u, [numel(mesh_u),1]);
mesh1d_v = reshape(mesh_v, [numel(mesh_v),1]);
mesh1d_w = reshape(mesh_w, [numel(mesh_w),1]);
%{
figure;
scatter3(mesh1d_u, mesh1d_v, mesh1d_w,10,mesh1d_w,'filled')
xlabel('u'); ylabel('v'); zlabel('w'); view([0,90]); daspect([1 1 1]);
%}
% spacing_xyz = (M1\[ all_w all_v all_u]')';
% Convert to xyz space 
spacing_xyz = (M1\[ mesh1d_w mesh1d_v mesh1d_u]')';
mesh1d_x = spacing_xyz(:,1);
mesh1d_y = spacing_xyz(:,2);
mesh1d_z = spacing_xyz(:,3);
%}

%{
spacing_x = -10:.1:10;
spacing_y = -10:.1:10;
spacing_z = -10:.1:10;
[mesh_x, mesh_y, mesh_z] = meshgrid(spacing_x, spacing_y, spacing_z);
mesh1d_x = reshape(mesh_x, [numel(mesh_x),1]);
mesh1d_y = reshape(mesh_y, [numel(mesh_y),1]);
mesh1d_z = reshape(mesh_z, [numel(mesh_z),1]);
nu = numel(spacing_x);
nv = numel(spacing_y);
nw = numel(spacing_z);
%}

%{
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh1d_z,'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]);
%}
% Convert to spherical space
[mesh1d_a, mesh1d_e, ~] = cart2sph(mesh1d_x, mesh1d_y, mesh1d_z);
mesh1d_a = rad2deg(mesh1d_a);
mesh1d_e = rad2deg(mesh1d_e);
%{
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh1d_a,'filled')
xlabel('x'); ylabel('y'); zlabel('z'); campos([0 0 0]); daspect([1 1 1]); colorbar
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh1d_e,'filled')
xlabel('x'); ylabel('y'); zlabel('z'); campos([0 0 0]); daspect([1 1 1]); colorbar
%}

mesh_a = reshape(mesh1d_a, [nv,nu,nw]);
mesh_e = reshape(mesh1d_e, [nv,nu,nw]);
%{
figure; 
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_a(:),'filled')
figure; 
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_e(:),'filled')
%}
% Find differences in each direction 
mesh_a_diff_x = abs(mesh_a - circshift(mesh_a,[-1,0,0]));
mesh_a_diff_y = abs(mesh_a - circshift(mesh_a,[0,-1,0]));
mesh_a_diff_z = abs(mesh_a - circshift(mesh_a,[0,0,-1]));
mesh_e_diff_x = abs(mesh_e - circshift(mesh_e,[-1,0,0]));
mesh_e_diff_y = abs(mesh_e - circshift(mesh_e,[0,-1,0]));
mesh_e_diff_z = abs(mesh_e - circshift(mesh_e,[0,0,-1]));
% Address edges 
mesh_a_diff_x(1,:,:)= mesh_a_diff_x(2,:,:);
mesh_a_diff_x(end,:,:)= mesh_a_diff_x(end-1,:,:);
mesh_a_diff_x(:,1,:)= mesh_a_diff_x(:,2,:);
mesh_a_diff_x(:,end,:)= mesh_a_diff_x(:,end-1,:);
mesh_a_diff_x(:,:,1)= mesh_a_diff_x(:,:,2);
mesh_a_diff_x(:,:,end)= mesh_a_diff_x(:,:,end-1);
mesh_a_diff_x(mesh_a_diff_x>300)= 360-mesh_a_diff_x(mesh_a_diff_x>300);
mesh_a_diff_x(mesh_a_diff_x<-300)= 360+mesh_a_diff_x(mesh_a_diff_x<-300);


mesh_a_diff_y(1,:,:)= mesh_a_diff_y(2,:,:);
mesh_a_diff_y(end,:,:)= mesh_a_diff_y(end-1,:,:);
mesh_a_diff_y(:,1,:)= mesh_a_diff_y(:,2,:);
mesh_a_diff_y(:,end,:)= mesh_a_diff_y(:,end-1,:);
mesh_a_diff_y(:,:,1)= mesh_a_diff_y(:,:,2);
mesh_a_diff_y(:,:,end)= mesh_a_diff_y(:,:,end-1);
mesh_a_diff_y(mesh_a_diff_y>300)= 360-mesh_a_diff_y(mesh_a_diff_y>300);
mesh_a_diff_y(mesh_a_diff_y<-300)= 360+mesh_a_diff_y(mesh_a_diff_y<-300);

mesh_a_diff_z(1,:,:)= mesh_a_diff_z(2,:,:);
mesh_a_diff_z(end,:,:)= mesh_a_diff_z(end-1,:,:);
mesh_a_diff_z(:,1,:)= mesh_a_diff_z(:,2,:);
mesh_a_diff_z(:,end,:)= mesh_a_diff_z(:,end-1,:);
mesh_a_diff_z(:,:,1)= mesh_a_diff_z(:,:,2);
mesh_a_diff_z(:,:,end)= mesh_a_diff_z(:,:,end-1);
mesh_a_diff_z(mesh_a_diff_z>300)= 360-mesh_a_diff_z(mesh_a_diff_z>300);
mesh_a_diff_z(mesh_a_diff_z<-300)= 360+mesh_a_diff_z(mesh_a_diff_z<-300);

mesh_e_diff_x(1,:,:)= mesh_e_diff_x(2,:,:);
mesh_e_diff_x(end,:,:)= mesh_e_diff_x(end-1,:,:);
mesh_e_diff_x(:,1,:)= mesh_e_diff_x(:,2,:);
mesh_e_diff_x(:,end,:)= mesh_e_diff_x(:,end-1,:);
mesh_e_diff_x(:,:,1)= mesh_e_diff_x(:,:,2);
mesh_e_diff_x(:,:,end)= mesh_e_diff_x(:,:,end-1);

mesh_e_diff_y(1,:,:)= mesh_e_diff_y(2,:,:);
mesh_e_diff_y(end,:,:)= mesh_e_diff_y(end-1,:,:);
mesh_e_diff_y(:,1,:)= mesh_e_diff_y(:,2,:);
mesh_e_diff_y(:,end,:)= mesh_e_diff_y(:,end-1,:);
mesh_e_diff_y(:,:,1)= mesh_e_diff_y(:,:,2);
mesh_e_diff_y(:,:,end)= mesh_e_diff_y(:,:,end-1);

mesh_e_diff_z(1,:,:)= mesh_e_diff_z(2,:,:);
mesh_e_diff_z(end,:,:)= mesh_e_diff_z(end-1,:,:);
mesh_e_diff_z(:,1,:)= mesh_e_diff_z(:,2,:);
mesh_e_diff_z(:,end,:)= mesh_e_diff_z(:,end-1,:);
mesh_e_diff_z(:,:,1)= mesh_e_diff_z(:,:,2);
mesh_e_diff_z(:,:,end)= mesh_e_diff_z(:,:,end-1);

%{
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_a_diff_x(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Azimuth difference x')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_a_diff_y(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Azimuth difference y')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_a_diff_z(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Azimuth difference z')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_e_diff_x(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Elevation difference x')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_e_diff_y(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Elevation difference y')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_e_diff_z(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Elevation difference z')
%}

% Find maxmimum 
mesh_a_diff = sqrt(mesh_a_diff_x(:).^2 +  mesh_a_diff_y(:).^2 + mesh_a_diff_z(:).^2);
mesh_e_diff = sqrt(mesh_e_diff_x(:).^2 +  mesh_e_diff_y(:).^2 + mesh_e_diff_z(:).^2);
%mesh_a_diff_max = max(max(abs(mesh_a_diff_x(:)), abs(mesh_a_diff_y(:))),abs(mesh_a_diff_z(:)));
mesh_a_diff = reshape(mesh_a_diff, size(mesh_a_diff_x));
%mesh_e_diff_max = max(max(abs(mesh_e_diff_x(:)), abs(mesh_e_diff_y(:))),abs(mesh_e_diff_z(:)));
mesh_e_diff = reshape(mesh_e_diff, size(mesh_e_diff_x));
%{
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_a_diff(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Azimuth difference')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_e_diff(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Elevation difference')
%}

% Find number of points 
mesh_npts_a = mesh_a_diff/.25; 
mesh_npts_e = mesh_e_diff/.25;
mesh_npts = mesh_npts_a.* mesh_npts_e; 
%{
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_npts_a(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Number of points Expected given a')
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_npts_e(:),'filled')
xlabel('x'); ylabel('y'); zlabel('z'); view([0,90]); daspect([1 1 1]); colorbar
title('Number of points Expected given e')
%}
%{
figure;
scatter3(mesh1d_x, mesh1d_y, mesh1d_z,10,mesh_npts(:),'filled')
xlabel('u'); ylabel('v'); zlabel('w'); view([0,90]); daspect([1 1 1]); colorbar
title('Number of points Expected')
%}
npts_2d = sum(mesh_npts,3);

%Vq = interp3(spacing_u, spacing_v,spacing_w, mesh_npts, all_u,all_v, all_w);
%scatter3(all_u, all_v,all_w, 10, Vq,'filled')

end

