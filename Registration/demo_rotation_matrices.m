% rotation testing 

x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

% YAW ONLY 
rz_min = -10;
rz_int = 2;
rz_max = 10;
rz = rz_min:rz_int:rz_max;
ry = zeros(size(rz));
rx = zeros(size(rz));
rx = deg2rad(rx);
ry = deg2rad(ry);
rz = deg2rad(rz);

n_r = numel(rz);

figure; 
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);
for r = 1:n_r;
    R = compose_rotation(rx(r),ry(r),rz(r));
    x_axis_R = (R*x_axis')';
    y_axis_R = (R*y_axis')';
    z_axis_R = (R*z_axis')';
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
end
title('Only YAW');

%% 

x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

% only ROLL/PITCH
rz_min = 0;
rz_int = 10;
rz_max = 50;
rx_min = -10;
rx_int = 2;
rx_max = 10;
ry_min = -10;
ry_int = 2;
ry_max = 10;
rx = rx_min:rx_int:rx_max;
ry = ry_min:ry_int:ry_max;
[rX, rY] = meshgrid(rx, ry);
rx = rX(:);
ry = rY(:);
rz = zeros(size(rx));
rx = deg2rad(rx);
ry = deg2rad(ry);
rz = deg2rad(rz);
n_r = numel(rz);

figure; 
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);
for r = 1:n_r;
    R = compose_rotation(rx(r),ry(r),rz(r));
    x_axis_R = (R*x_axis')';
    y_axis_R = (R*y_axis')';
    z_axis_R = (R*z_axis')';
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
end
title('Only ROLL/PITCH');
%% 

x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

% small YAW/ROLL/PITCH

rx_min = -10;
rx_int = 3;
rx_max = 10;
ry_min = -10;
ry_int = 3;
ry_max = 10;
rz_min = -10;
rz_int = 3;
rz_max = 10;
rx = rx_min:rx_int:rx_max;
ry = ry_min:ry_int:ry_max;
rz = rz_min:rz_int:rz_max;
[rX, rY, rZ] = meshgrid(rx, ry, rz);
rx = rX(:);
ry = rY(:);
rz = rZ(:);
rx = deg2rad(rx);
ry = deg2rad(ry);
rz = deg2rad(rz);
n_r = numel(rz);

figure; 
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);
for r = 1:n_r;
    R = compose_rotation(rx(r),ry(r),rz(r));
    x_axis_R = (R*x_axis')';
    y_axis_R = (R*y_axis')';
    z_axis_R = (R*z_axis')';
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
end
title('Small YAW/ROLL/PITCH');
%% 

x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

% small ROLL/PITCH + YAW

rx_min = -10;
rx_int = 3;
rx_max = 10;
ry_min = -10;
ry_int = 3;
ry_max = 10;
rz_min = -180;
rz_int = 45;
rz_max = 180;
rx = rx_min:rx_int:rx_max;
ry = ry_min:ry_int:ry_max;
rz = rz_min:rz_int:rz_max;
[rX, rY, rZ] = meshgrid(rx, ry, rz);
rx = rX(:);
ry = rY(:);
rz = rZ(:);
rx = deg2rad(rx);
ry = deg2rad(ry);
rz = deg2rad(rz);
n_r = numel(rz);

figure; 
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);
for r = 1:n_r;
    R = compose_rotation(rx(r),ry(r),rz(r));
    x_axis_R = (R*x_axis')';
    y_axis_R = (R*y_axis')';
    z_axis_R = (R*z_axis')';
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
end
title('Small ROLL/PITCH + YAW');

%% 

x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

% big YAW/ROLL/PITCH

rx_min = -180;
rx_int = 45;
rx_max = 180;
ry_min = -180;
ry_int = 45;
ry_max = 180;
rz_min = -180;
rz_int = 45;
rz_max = 180;
rx = rx_min:rx_int:rx_max;
ry = ry_min:ry_int:ry_max;
rz = rz_min:rz_int:rz_max;
[rX, rY, rZ] = meshgrid(rx, ry, rz);
rx = rX(:);
ry = rY(:);
rz = rZ(:);
rx = deg2rad(rx);
ry = deg2rad(ry);
rz = deg2rad(rz);
n_r = numel(rz);

figure; 
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);
for r = 1:n_r;
    R = compose_rotation(rx(r),ry(r),rz(r));
    x_axis_R = (R*x_axis')';
    y_axis_R = (R*y_axis')';
    z_axis_R = (R*z_axis')';
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
end
title('Big YAW/ROLL/PITCH');

%% 

x_axis = [0 0 0; 1 0 0];
y_axis = [0 0 0; 0 1 0];
z_axis = [0 0 0; 0 0 1];

% small YAW/ROLL/PITCH with TRANSLATION

rx_min = -10;
rx_int = 3;
rx_max = 10;
ry_min = -10;
ry_int = 3;
ry_max = 10;
rz_min = -10;
rz_int = 3;
rz_max = 10;
rx = rx_min:rx_int:rx_max;
ry = ry_min:ry_int:ry_max;
rz = rz_min:rz_int:rz_max;
[rX, rY, rZ] = meshgrid(rx, ry, rz);
rx = rX(:);
ry = rY(:);
rz = rZ(:);
rx = deg2rad(rx);
ry = deg2rad(ry);
rz = deg2rad(rz);
n_r = numel(rz);

t = [ 2 1 0];

figure; 
hold on
plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'-r', 'linewidth',2)
plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'-g', 'linewidth',2)
plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'-b', 'linewidth',2)
xlabel('x');
ylabel('y');
zlabel('z');
view(65,20);
for r = 1:n_r;
    R = compose_rotation(rx(r),ry(r),rz(r));
    x_axis_R = (R*x_axis')' + repmat(t, [2,1]);
    y_axis_R = (R*y_axis')' + repmat(t, [2,1]);
    z_axis_R = (R*z_axis')' + repmat(t, [2,1]);
    plot3(x_axis_R(:,1),x_axis_R(:,2),x_axis_R(:,3),'-r', 'linewidth',2)
    plot3(y_axis_R(:,1),y_axis_R(:,2),y_axis_R(:,3),'-g', 'linewidth',2)
    plot3(z_axis_R(:,1),z_axis_R(:,2),z_axis_R(:,3),'-b', 'linewidth',2)
end
title('Small YAW/ROLL/PITCH + TRANSLATION');
