%% Toy Example
% Based on Harvey's suggestion
% Only allows matches for points imaged by all cameras
% Requires LM_toy_nest2.m
%clear all; close all, clc;
set(0,'defaultfigureposition', [895   169   760   651]')
%% Create points
fprintf('\nCreate points\n');
% Point parameters
n_P = 15;
P_xmin = -10;
P_xmax = 10;
P_ymin = -10;
P_ymax = 10;
P_zmin = -1;
P_zmax = 1;
P_rmin = .04;
P_rmax = .5;

% Bookkeeping
P_index = 1:n_P;

% Random locations on uniform distribution
Px = unifrnd(P_xmin, P_xmax, [n_P,1]);
Py = unifrnd(P_ymin, P_ymax, [n_P,1]);
Pz = unifrnd(P_zmin, P_zmax, [n_P,1]);
Pr = unifrnd(P_rmin, P_rmax, [n_P,1]);

%
figure;
plot3(Px,Py,Pz,'+k','markersize',10,...
    'markerfacecolor',[0 0 0], 'linewidth',1.2);axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
title('Truth tree locations');
grid on
%legend('truth location','location','bestoutside');
%}

% Outputs
% Px            x locations of truth points
% Py            y locations of truth points
% Pz            z locations of truth points
% Pr            r of truth points
clear P_rmin P_rmax
%% Create sensors
fprintf('\nCreate sensors\n');

% Number of sensors
n_S = 3;

% Sensor parameters
T_xmin = -7;
T_xmax = 7;
T_ymin = -7;
T_ymax = 7;
T_zmin = -.5;
T_zmax = .5;
% If rotation angles are large, change matching function
R_xmin = deg2rad(-10);
R_xmax = deg2rad(10);
R_ymin = deg2rad(-10);
R_ymax = deg2rad(10);
R_zmin = deg2rad(-10);
R_zmax = deg2rad(10);

% Colormap for sensors
color = jet(n_S);

% Rotation and translation truth
rt = zeros(n_S,3);
tt = zeros(n_S,3);
for s = 1:n_S;
    rt(s,1) = unifrnd(R_xmin, R_xmax,1);
    rt(s,2) = unifrnd(R_ymin, R_ymax,1);
    rt(s,3) = unifrnd(R_zmin, R_zmax,1);
    tt(s,1) = unifrnd(T_xmin, T_xmax,1);
    tt(s,2) = unifrnd(T_ymin, T_ymax,1);
    tt(s,3) = unifrnd(T_zmin, T_zmax,1);
end

% Set S1 to origin
rt(1,:) = 0;
tt(1,:) = 0;

%
figure;
plot3(Px,Py,Pz,'+k','markersize',10,...
    'markerfacecolor',[0 0 0], 'linewidth',1.2);axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
hold on
legend_str{1} = 'Truth tree locations';
for s = 1:n_S;
    legend_str{s+1} = sprintf('Camera %g',s);
    plot3(tt(s,1),tt(s,2),tt(s,3),'^k','markersize',10,...
        'markerfacecolor',color(s,:));
end
title('Points and camera positions');
legend(legend_str,'location','best');
grid on
%}

% Outputs
% S             structure of sensors
% tr            truth rotation
% tt            truth translation
% S.c           color
clear R_xmax R_xmin R_ymax R_ymin R_zmax R_zmin
clear T_xmax T_xmin T_ymax T_ymin T_zmax T_zmin
clear map legend_str
clear rx ry rz s tx ty tz
%% Image points with sensors
fprintf('\nImage points with sensors\n');

% Number of points visible to each sensor
V_min = n_P;
V_max = n_P;

% Noise in xyz location
N_min = -.2;
N_max = .2;

N_rmin = -0;
N_rmax = 0;
N_tmin = -.0;
N_tmax = .0;

P_truthix = cell(n_S,1);
for s = 1:n_S;
    % Find visible
    n_visible = randi([V_min V_max],1);
    ix_rand = randperm(n_P)';
    is_visible = ix_rand(1:n_visible);
    P_truthix{s} = is_visible;
    % Add noise
    Nx = unifrnd(N_min, N_max,[n_visible,1]);
    Ny = unifrnd(N_min, N_max,[n_visible,1]);
    Nz = unifrnd(N_min, N_max,[n_visible,1]);
    P_WCS = [Px(is_visible)'; Py(is_visible)'; Pz(is_visible)'];
    % Transform back to LCS
    R = compose_rotation(rt(s,1)+deg2rad(unifrnd(N_rmin, N_rmax,1)),...
        rt(s,2)+deg2rad(unifrnd(N_rmin, N_rmax,1)),...
        rt(s,3)+deg2rad(unifrnd(N_rmin, N_rmax,1)));
    T = repmat(tt(s,:)' + unifrnd(N_tmin, N_tmax, [3,1]),[1,n_visible]);
    P_LCS{s} = (R'*(P_WCS - T)) - [Nx Ny Nz]';
    Rt = compose_rotation(rt(s,1),rt(s,2),rt(s,3));
    Tt = repmat(tt(s,:)',[1,n_visible]);
    P_LCSt{s} = (Rt'*(P_WCS - Tt));
    P_rad{s} = Pr(is_visible);
    P_n(s) = n_visible;
end

% Individual camera views
% color_P_index = jet(n_P);
for s = 1:n_S;
    figure;
    hold on
    plot3(P_LCSt{s}(1,:),P_LCSt{s}(2,:),P_LCSt{s}(3,:),'+k','markersize',10,...
        'markerfacecolor',color(s,:));
    plot3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),'ok','markersize',5,...
        'markerfacecolor',color(s,:));
    %scatter3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),30,...
    %    color_P_index(truth_P_index{s},:),'filled');
    axis(1.5*[P_xmin P_xmax P_ymin P_ymax -10 10]);
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    grid on
    titlestr = sprintf('Camera %g',s);
    title(titlestr);
    plot3(0,0,0,'^k','markersize',10,...
        'markerfacecolor',color(s,:));
    legend_str{1} = 'Truth points';
    legend_str{2} = 'Points in LCS';
    legend_str{3} = 'Camera';
    legend(legend_str,'location','best');
end
%}
% Outputs
% P_LCS                 points in local coordinate system
% P_rad                 radius of points
% P_n                   number of points
clear V_min V_max N_min N_max Nx Ny Nz P_WCS
clear s n_visible ix_rand is_visible
%% Run registration code
fprintf('\nRun registration code\n');

match_i = cell(n_S); %base
match_j = cell(n_S); %mobile
match_R = cell(n_S);
match_t = cell(n_S);

% Find rotation and tranlation from j to i along with point matching pairs
for i = 1:1;%n_S;
    for j = i:n_S;
        fprintf('\n\tMatching %g to %g\n',j,i);
        [ match_R{i,j}, match_t{i,j}, match_i{i,j}, match_j{i,j} ] = toy_registrationfunction(P_LCS{i}',P_LCS{j}',P_rad{i},P_rad{j});
    end
end

% Find points j transformed into new coordinate system i
match_Pi = cell(n_S,n_S);           % Points that are matched
match_Pi_all = cell(n_S,n_S);       % All points
for i = 1:1;%n_S;
    for j = i:n_S;
        match_Pi{i,j} = match_R{i,j}*P_LCS{j}(:,match_j{i,j}) + repmat(match_t{i,j},[1,numel(match_j{i,j})]);
        match_Pi_all{i,j} = match_R{i,j}*P_LCS{j} + repmat(match_t{i,j},[1,P_n(j)]);
    end
end

%
% Transformed to camera 1
figure;
legend_str{1} = 'Truth points';
plot3(Px,Py,Pz,'+k','markersize',10,...
    'markerfacecolor',[0 0 0], 'linewidth',1.2);xlabel('x');ylabel('y');zlabel('z');
view(0,90);
hold on
for s = 1:n_S;
    legend_str{s+1} = sprintf('Camera %g',s);
    plot3(tt(s,1),tt(s,2),tt(s,3),'^k','markersize',10,...
        'markerfacecolor',color(s,:));
end
for s = 1:n_S;
    legend_str{n_S+1+s} = sprintf('Points %g in WCS',s);
    plot3(match_Pi{1,s}(1,:),match_Pi{1,s}(2,:),match_Pi{1,s}(3,:),'ok','markersize',5,...
        'markerfacecolor',color(s,:));
end
legend(legend_str,'location','best');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax -10 10]);
title('Points in WCS - Initial Transformation');
grid on
%}

% Outputs
% match_i               index to points in i (world)
% match_j               index to points in j (local)
% match_R               (i,j) gives rotation matrix of j into i
% match_t               (i,j) gives translation matrix of j into i
% match_Pi              matched points of j transformed into CS of i
% match_Pi_all          all points of j transformed into CS of i
%% Find unique points
fprintf('\nFind unique points\n');

% Gather all points in WCS (i = 1) from each camera
all_Pi_x = [];
all_Pi_y = [];
all_Pi_z = [];
all_Pi_r = [];
for j = 1:n_S;
    all_Pi_x = [all_Pi_x match_Pi_all{1,j}(1,:)]; %#ok<*AGROW>
    all_Pi_y = [all_Pi_y match_Pi_all{1,j}(2,:)];
    all_Pi_z = [all_Pi_z match_Pi_all{1,j}(3,:)];
    all_Pi_r = [all_Pi_r P_rad{j}'];
end
%
% Find error between points
n_all = numel(all_Pi_x);
A_Pi_x = repmat(all_Pi_x,[n_all,1]);
A_Pi_xt = repmat(all_Pi_x',[1,n_all]);
A_Pi_y = repmat(all_Pi_y,[n_all,1]);
A_Pi_yt = repmat(all_Pi_y',[1,n_all]);
A_Pi_z = repmat(all_Pi_z,[n_all,1]);
A_Pi_zt = repmat(all_Pi_z',[1,n_all]);
A_Pi_r = repmat(all_Pi_r,[n_all,1]);
A_Pi_rt = repmat(all_Pi_r',[1,n_all]);

% Points are the matches if within t_xyz and t_r
t_xyz = 1.^2;
t_r = 0.1;
D2_xyz = (A_Pi_x - A_Pi_xt).^2 + (A_Pi_y - A_Pi_yt).^2 + (A_Pi_z - A_Pi_zt).^2;
D2_r = (A_Pi_r - A_Pi_rt).^2;
D2_isxyz = (D2_xyz<t_xyz);
D2_isr = (D2_r<t_r);
D2_is = D2_isxyz&D2_isr&triu(D2_isr,1);

% Block out ii matches
D2_ii = false(size(D2_is));
ix = 1;
for s = 1:n_S;
    D2_ii(ix:ix+P_n(s)-1,ix:ix+P_n(s)-1) = true;
    ix = ix + P_n(s);
end
% Final valid matches
D2_is = D2_is&~D2_ii;

%{
figure;
imagesc(D2_is);
axis image
%}

% Find unique points
is_unique = (sum(D2_is,1)==0);
unique_x = all_Pi_x(is_unique);
unique_y = all_Pi_y(is_unique);
unique_z = all_Pi_z(is_unique);
unique_r = all_Pi_r(is_unique);
n_unique = numel(unique_x);

%
figure;
hold on
plot3(Px,Py,Pz,'+k','markersize',10,...
    'markerfacecolor',[0 0 0], 'linewidth',1.2);
scatter3(all_Pi_x, all_Pi_y, all_Pi_z,10,all_Pi_r,'filled');
scatter3(unique_x,unique_y,unique_z,60,'k');%,'markersize',20,...
%  'markerfacecolor',[1 1 1],'alpha',1, 'linewidth',.5);
legend_str{1} = 'Truth points';
legend_str{2} = 'All points';
legend_str{3} = 'Unique points';
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
legend(legend_str,'location','best');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax 2*P_zmin 2*P_zmax]);
title('Unique points ');
%}

% Output
% unique_x                  x values of unique points
% unique_y                  y values of unique points
% unique_z                  z values of unique points
% unique_r                  r values of unique points
% n_unique                  number of unique points
clear all_Pi_x all_Pi_y all_Pi_z all_Pi_r
clear A_Pi_xt A_Pi_yt A_Pi_zt A_Pi_rt
clear A_Pi_x A_Pi_y A_Pi_z A_Pi_r
clear D2_xyz D2_r D2_isxyz D2_isr D2_is
clear is_unique
%% Build up data array
fprintf('\nBuild up data array\n');

data_LCSi = nan(n_unique,n_S,3);
data_truthix = nan(n_unique,n_S);

% Points in WCS
for s = 1:n_S;
    P_WCS{s} = match_Pi_all{1,s};
end

for s = 1:n_S;
    % Matrices for distance calculations
    A_unique_x = repmat(unique_x',[1,P_n(s)]);
    A_unique_y = repmat(unique_y',[1,P_n(s)]);
    A_unique_z = repmat(unique_z',[1,P_n(s)]);
    A_unique_r = repmat(unique_r',[1,P_n(s)]);
    A_Pi_x = repmat(P_WCS{s}(1,:),[n_unique,1]);
    A_Pi_y = repmat(P_WCS{s}(2,:),[n_unique,1]);
    A_Pi_z = repmat(P_WCS{s}(3,:),[n_unique,1]);
    A_Pi_r = repmat(P_rad{s}',[n_unique,1]);
    % Match if within thresholds t_xyz t_r
    t_xyz = 1.^2;
    t_r = 0.1;
    D2_xyz = (A_Pi_x - A_unique_x).^2 + (A_Pi_y - A_unique_y).^2 + (A_Pi_z - A_unique_z).^2;
    D2_r = (A_Pi_r - A_unique_r).^2;
    D2_isxyz = (D2_xyz<t_xyz);
    D2_isr = (D2_r<t_r);
    D2_is = D2_isxyz&D2_isr;
    
    % Remove duplicates
    ix_dup = find(sum(D2_is,2)>1);
    for d = 1:numel(ix_dup);
        i = ix_dup(d);
        j = find(D2_is(i,:));
        [~,minix] = min(D2_xyz(i,j)+D2_r(i,j));
        D2_is(i,:) = false;
        D2_is(i,j(minix)) = true;
    end
    
    %{
    figure;
    imagesc(D2_is);
    axis image
    xlabel(sprintf('Camera %g points',s));
    ylabel('Unique points');
    %}
    [u,c] = find(D2_is);
    
    data_LCSi(u,s,1) = P_LCS{s}(1,c);
    data_LCSi(u,s,2) = P_LCS{s}(2,c);
    data_LCSi(u,s,3) = P_LCS{s}(3,c);
    data_truthix(u,s) = P_truthix{s}(c); % Should be identical for all sensors
    
end

data_WCSt = nan(n_unique,n_S,3);
for s = 1:n_S;
    data_WCSt(:,s,:) = [Px(data_truthix(:,s)) Py(data_truthix(:,s)) Pz(data_truthix(:,s))];
end

data_LCSt = nan(n_unique,n_S,3);
for s = 1:n_S;
    R = compose_rotation(rt(s,1),rt(s,2),rt(s,3));
    T = tt(s,:)';
    data_LCSt(:,s,:) = (R'*(squeeze(data_WCSt(:,s,:))' - repmat(T,[1,n_unique])))';
end

% Outputs
% S.P_WCS                   points in WCS (camera 1)
% data                      n_unique x n_S x 3 of LCS points
clear all_Pi_x all_Pi_y all_Pi_z all_Pi_r
clear A_Pi_x A_Pi_y A_Pi_z A_Pi_r
clear A_unique_x A_unique_y A_unique_z A_unique_r
clear t_xyz t_r
clear D2_xyz D2_r D2_isxyz D2_isr D2_is
clear ix_dup i j minix u c
clear is_unique
%% Levenberg Marquardt
fprintf('\nLevenberg Marquardt\n');

% Add transformation parameters to parameter vector
P0 = [];
ri = zeros(n_S,3);
ti = zeros(n_S,3);
for s = 2:n_S;
    [rx,ry,rz] = decompose_rotation(match_R{1,s});
    ri(s,1) = rx;
    ri(s,2) = ry;
    ri(s,3) = rz;
    t = match_t{1,s};
    ti(s,:) = t;
    P0 = [P0 rx ry rz t(1) t(2) t(3)];
end

% Add position 1 data (exclude entries with NaN's)
is_nan = any(isnan(data_LCSi(:,:,1)),2);
n_isnnan = sum(~is_nan);
P0 = [P0 data_LCSi(~is_nan,1,1)' data_LCSi(~is_nan,1,2)' data_LCSi(~is_nan,1,3)'];

% Levenberg-Marquardt (nested)
options = optimset('lsqnonlin');
options.Display = 'iter-detailed';
options.TolFun = 1e-10;
options.TolX = 1e-10;
[P, exitflag, output] = LM_toy_nest2(P0,data_LCSi(~is_nan,:,:),options);

clear rx ry rz t s
%% Reconstruct

% Find registration and translation matrices
ix = 1;
rf = zeros(n_S,3);
tf = zeros(n_S,3);
for s = 2:n_S;
    rx = P(ix);
    ry = P(ix+1);
    rz = P(ix+2);
    rf(s,1) = rx;
    rf(s,2) = ry;
    rf(s,3) = rz;
    tf(s,:) = P(ix+3:ix+5)';
    ix = ix + 6;
end

% Points in WCS before LM optimization
data_WCSi = zeros(n_isnnan,n_S,3);
data_WCSi(:,1,:) = reshape(P0(ix:end),[(numel(P0)-ix+1)/3,3]);
for s = 2:n_S;
    R = compose_rotation(ri(s,1),ri(s,2),ri(s,3));
    T = ti(s,:)';
    data_WCSi(:,s,:) = (R*squeeze(data_LCSi(~is_nan,s,:))'+repmat(T,[1,n_isnnan]))';
end

% Points in WCS after LM optimization
data_WCSf = zeros(n_isnnan,n_S,3);
data_WCSf(:,1,:) = reshape(P(ix:end),[(numel(P)-ix+1)/3,3]);
for s = 2:n_S;
    R = compose_rotation(rf(s,1),rf(s,2),rf(s,3));
    T = tf(s,:)';
    data_WCSf(:,s,:) = (R*squeeze(data_LCSi(~is_nan,s,:))'+repmat(T,[1,n_isnnan]))';
end

%
% Points in LCS after LM optimization
data_LCSf = zeros(n_isnnan,n_S,3);
data_LCSf(:,1,:) = reshape(P(ix:end),[(numel(P)-ix+1)/3,3]);
for s = 2:n_S;
    R = compose_rotation(rf(s,1),rf(s,2),rf(s,3));
    T = tf(s,:)';
    data_LCSf(:,s,:) = (R'*(squeeze(data_WCSf(:,1,:))' - repmat(T,[1,n_isnnan])))';
end
%}

%{
% Points in LCS after LM optimization
data_LCSf = zeros(n_isnnan,n_S,3);
data_LCSf(:,1,:) = reshape(P(ix:end),[(numel(P)-ix+1)/3,3]);
for s = 2:n_S;
    R = compose_rotation(rf(s,1),rf(s,2),rf(s,3));
    T = tf(s,:)';
    data_LCSf(:,s,:) = (R'*(squeeze(data_WCSf(:,s,:))' - repmat(T,[1,n_isnnan])))';
end
%}

clear legend_str
figure;
legend_str{1} = 'Truth points';
plot3(Px,Py,Pz,'+k','markersize',10,...
    'markerfacecolor',[0 0 0], 'linewidth',1.2);xlabel('x');ylabel('y');zlabel('z');
view(0,90);
hold on
for s = 1:n_S
    plot3(data_WCSi(:,s,1),data_WCSi(:,s,2),data_WCSi(:,s,3),'ok','markersize',5,...
        'markerfacecolor','b');%color(s,:));
    plot3(data_WCSf(:,s,1),data_WCSf(:,s,2),data_WCSf(:,s,3),'ok','markersize',5,...
        'markerfacecolor','r');%color(s,:));
    if s == 1;
        legend_str{2} = 'Initial Transformation';
        legend_str{3} = 'Nonlinear Optimization';
    end
end
legend(legend_str,'location','best');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax 2*P_zmin 2*P_zmax]);
title('Points in WCS');
grid on


% Outputs
% S.Ri                  initial rotation
% S.ti                  initial translation
% S.Rf                  final rotation after LM
% S.tf                  final translation after LM
% data_i                initial WCS
% data_f                final WCS
%   data                  LCS
clear rx ry rz ix
%% Plot LCS points
%
for s = 1:n_S;
    figure;
    hold on
    plot3(0,0,0,'^k','markersize',10,...
        'markerfacecolor',color(s,:));
    legend_str{1} = 'Camera';
    legend_str{2} = 'Truth points';
    plot3(data_LCSt(:,s,1),data_LCSt(:,s,2),data_LCSt(:,s,3),'+k','markersize',15,...
        'markerfacecolor','k');%color(s,:));xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    plot3(data_LCSi(:,s,1),data_LCSi(:,s,2),data_LCSi(:,s,3),'ok','markersize',5,...
        'markerfacecolor','b');%color(s,:));
    plot3(data_LCSf(:,s,1),data_LCSf(:,s,2),data_LCSf(:,s,3),'ok','markersize',5,...
        'markerfacecolor','r');%color(s,:));
    legend_str{3} = 'Initial Transformation';
    legend_str{4} = 'Nonlinear Optimization';
    legend(legend_str,'location','best');
    axis equal
    axis(1.5*[P_xmin P_xmax P_ymin P_ymax 2*P_zmin 2*P_zmax]);
    title(sprintf('Points in LCS - Camera %g', s));
    grid on
end
%}
%% Error metrics
sqerrori_p = sum((data_WCSi - data_WCSt).^2,3);
sqerrorf_p = sum((data_WCSf - data_WCSt).^2,3);

rmsei_p_point = sqrt(mean(sqerrori_p,2));
rmsef_p_point = sqrt(mean(sqerrorf_p,2));

rmsei_p_sensor = sqrt(mean(sqerrori_p,1));
rmsef_p_sensor = sqrt(mean(sqerrorf_p,1));

rmsei_p_all = sqrt(mean(sqerrori_p(:)));
rmsef_p_all = sqrt(mean(sqerrorf_p(:)));

fprintf('\nAll errors\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
fprintf('%3.3f \t%3.3f \tRMSE between points [meters]\n', rmsei_p_all, rmsef_p_all);

fprintf('\nErrors by point\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for p = 1:n_isnnan;
    fprintf('%3.3f \t%3.3f \tRMSE between points %g\n', rmsei_p_point(p), rmsef_p_point(p), p);
end

fprintf('\nErrors by sensor\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for s = 1:n_S;
    fprintf('%3.3f \t%3.3f \tRMSE for sensor %g\n', rmsei_p_sensor(s), rmsef_p_sensor(s), s);
end

rmsei_t_sensor = sqrt(mean((ti - tt).^2,2));
rmsef_t_sensor = sqrt(mean((tf - tt).^2,2));

rmsei_t_all = mean(rmsei_t_sensor);
rmsef_t_all = mean(rmsef_t_sensor);

rmsei_r_sensor = rad2deg(sqrt(mean((ri - rt).^2,2)));
rmsef_r_sensor = rad2deg(sqrt(mean((rf - rt).^2,2)));

rmsei_r_all = sum(rmsei_r_sensor);
rmsef_r_all = sum(rmsef_r_sensor);

fprintf('\nAll errors\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
fprintf('%3.3f \t%3.3f \tRMSE of rotation comp. to truth [degrees]\n', rmsei_r_all, rmsef_r_all);
fprintf('%3.3f \t%3.3f \tRMSE of translation comp. to truth [m]\n', rmsei_t_all, rmsef_t_all);

fprintf('\nErrors in Rotation [degrees]\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for s = 1:n_S;
    fprintf('%3.3f \t%3.3f \tSensor %g\n', rmsei_r_sensor(s), rmsef_r_sensor(s), s);
end

fprintf('\nErrors in Translation [meters]\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for s = 1:n_S;
    fprintf('%3.3f \t%3.3f \tSensor %g\n', rmsei_t_sensor(s), rmsef_t_sensor(s),s);
end

%{
% No a priori knowledge
errori_p = sqrt(sum((data_WCSi - repmat(mean(data_WCSi,2),[1,n_S,1])).^2,3));
errorf_p = sqrt(sum((data_WCSf - repmat(mean(data_WCSf,2),[1,n_S,1])).^2,3));

errori_p_point = sum(errori_p,2);
errorf_p_point = sum(errorf_p,2);

errori_p_sensor = sum(errori_p,1);
errorf_p_sensor = sum(errorf_p,1);

errori_p_all = sum(errori_p(:));
errorf_p_all = sum(errorf_p(:));

fprintf('\nAll errors\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
fprintf('%3.3f \t%3.3f \tRMSE between points [meters]\n', errori_p_all, errorf_p_all);

fprintf('\nErrors by point\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for p = 1:n_isnnan;
fprintf('%3.3f \t%3.3f \tRMSE between points %g\n', errori_p_point(p), errorf_p_point(p), p);
end

fprintf('\nErrors by sensor\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for s = 1:n_S;
fprintf('%3.3f \t%3.3f \tRMSE for sensor %g\n', errori_p_sensor(s), errorf_p_sensor(s), s);
end
%}


