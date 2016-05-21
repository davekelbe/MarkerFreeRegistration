

x = [0 1; 0 0; 0 0];
y = [0 0; 0 1; 0 0];
z = [0 0; 0 0; 0 1];

xyz = [x y z];

rx = 20;
ry = 0;
rz = 0;
R = compose_rotation(deg2rad(rx),deg2rad(ry),deg2rad(rz));
xyzR = R*xyz;

lw = 3;
figure;
hold on
xlabel('x');
ylabel('y');
zlabel('z');
plot3(xyz(1,1:2),xyz(2,1:2),xyz(3,1:2),'-r', 'linewidth', lw);
plot3(xyz(1,3:4),xyz(2,3:4),xyz(3,3:4),'-g', 'linewidth', lw);
plot3(xyz(1,5:6),xyz(2,5:6),xyz(3,5:6),'-b', 'linewidth', lw);
view(60,20);

rx = 20; % around x axis 
ry = 0;
rz = 0;
R = compose_rotation(deg2rad(rx),deg2rad(ry),deg2rad(rz));
xyzR = R*xyz;
plot3(xyzR(1,1:2),xyzR(2,1:2),xyzR(3,1:2),'--r', 'linewidth', lw);
plot3(xyzR(1,3:4),xyzR(2,3:4),xyzR(3,3:4),'--g', 'linewidth', lw);
plot3(xyzR(1,5:6),xyzR(2,5:6),xyzR(3,5:6),'--b', 'linewidth', lw);

rx = 0;
ry = 20; % around y axis 
rz = 0;
R = compose_rotation(deg2rad(rx),deg2rad(ry),deg2rad(rz));
xyzR = R*xyz;
plot3(xyzR(1,1:2),xyzR(2,1:2),xyzR(3,1:2),':r', 'linewidth', lw);
plot3(xyzR(1,3:4),xyzR(2,3:4),xyzR(3,3:4),':g', 'linewidth', lw);
plot3(xyzR(1,5:6),xyzR(2,5:6),xyzR(3,5:6),':b', 'linewidth', lw);

rx = 0;
ry = 0;
rz = 50; % around z axis (azimuth pointing)
R = compose_rotation(deg2rad(rx),deg2rad(ry),deg2rad(rz));
xyzR = R*xyz;
plot3(xyzR(1,1:2),xyzR(2,1:2),xyzR(3,1:2),'-.r', 'linewidth', lw);
plot3(xyzR(1,3:4),xyzR(2,3:4),xyzR(3,3:4),'-.g', 'linewidth', lw);
plot3(xyzR(1,5:6),xyzR(2,5:6),xyzR(3,5:6),'-.b', 'linewidth', lw);

rx = 5;
ry = 5;
rz = 70; % around z axis (azimuth pointing)
R = compose_rotation(deg2rad(rx),deg2rad(ry),deg2rad(rz));
xyzR = R*xyz;
plot3(xyzR(1,1:2),xyzR(2,1:2),xyzR(3,1:2),'-.r', 'linewidth', lw);
plot3(xyzR(1,3:4),xyzR(2,3:4),xyzR(3,3:4),'-.g', 'linewidth', lw);
plot3(xyzR(1,5:6),xyzR(2,5:6),xyzR(3,5:6),'-.b', 'linewidth', lw);


