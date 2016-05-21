function [  ] = kelbe_product_errfield(  )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;
space =2;
xv = xmin:space:xmax;
yv = ymin:space:ymax;
zv = zmin:space:zmax;
[X,Y,Z] = meshgrid(xv, yv, zv);
[nx,ny,nz] = size(X);
test_rx = 0;
test_ry = 3;
test_rz = 1;
test_tx = .1;
test_ty = 0.5;
test_tz = 0.0;
test_R = compose_rotation(deg2rad(test_rx),deg2rad(test_ry),deg2rad(test_rz));
test_t = [test_tx test_ty test_tz]';
XYZ = [X(:),Y(:),Z(:)]';
n_pts = size(XYZ,2);
XYZhat = test_R*XYZ + repmat(test_t,[1, n_pts]);
error_XYZ = XYZ - XYZhat;
error = abs(mean((XYZ - XYZhat),1));

figure;
hold on
scatter3(XYZ(1,:),XYZ(2,:),XYZ(3,:),10,'r','filled');
scatter3(XYZhat(1,:),XYZhat(2,:),XYZhat(3,:),10,'b','filled');
axis equal
legend('True','Estimated');

figure;
hold on
scatter3(XYZ(1,:),XYZ(2,:),XYZ(3,:),10,error,'filled');
axis equal
colorbar
emin = min(error);
emax = max(error);
cmin = emin -.1;
cmax = emax + .1;
caxis([cmin, cmax]);

figure;
xbin = linspace(min(error),max(error),100);
count = hist(error, xbin);
plot(xbin,count,'-x');

figure;
scale = 2;
quiver3(XYZ(1,:),XYZ(2,:),XYZ(3,:),error_XYZ(1,:), error_XYZ(2,:), error_XYZ(3,:),...
    scale);



end

