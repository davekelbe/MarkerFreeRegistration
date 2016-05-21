% Limits of full tie point set
xmin = -30;
xmax = 30;
ymin = -30;
ymax = 30;
zmin = -1;
zmax = 1;
n_T = 250;

% Radii
rmin = .01;
rmax = .5;

% Tie point locations and radi
T_x = xmin + (xmax-xmin)*rand(n_T,1);
T_y = ymin + (ymax-ymin)*rand(n_T,1);
T_z = zmin + (zmax-zmin)*rand(n_T,1);
T_r = rmin + (rmax-rmin)*rand(n_T,1);

% Remove stems close to each other
t_dist_min = 2;
T_xrep1 = repmat(T_x,[1, n_T]);
T_yrep1 = repmat(T_y,[1, n_T]);
T_xrep2 = repmat(T_x',[n_T,1]);
T_yrep2 = repmat(T_y',[n_T,1]);
T_dist = sqrt((T_xrep1 - T_xrep2).^2 + (T_yrep1 - T_yrep2).^2);
is_close = (T_dist<t_dist_min);
is_close(logical(tril(is_close,0))) = false;
[close1, close2] = find(is_close);
T_x(close1) = nan;
T_y(close1) = nan;
T_z(close1) = nan;
T_r(close1) = nan;
T_x = T_x(~isnan(T_x));
T_y = T_y(~isnan(T_y));
T_z = T_z(~isnan(T_z));
T_r = T_r(~isnan(T_r));

% Sensors
n_Sg = 5;
n_S = round(n_Sg.^2);

% Nominal position
S_tx = repmat((-20:10:20),[n_Sg,1]);
S_ty = repmat((-20:10:20)',[1,n_Sg]);
S_tx = S_tx(:);
S_ty = S_ty(:);
S_tz = zeros(n_S,1);

% Deviation around nominal position
dxmin = -2;
dxmax = 2;
dymin = -2;
dymax = 2;
dzmin = -.5;
dzmax = .5;
S_dx = dxmin + (dxmax-dxmin)*rand(n_S,1);
S_dy = dymin + (dymax-dymin)*rand(n_S,1);
S_dz = dzmin + (dzmax-dzmin)*rand(n_S,1);
S_tx = S_tx + S_dx;
S_ty = S_ty + S_dy;
S_tz = S_tz + S_dz;

% Nominal Rotation angles
S_rx = zeros(n_S,1);
S_ry = zeros(n_S,1);
S_rz = zeros(n_S,1);

% Deviation of rotation angles
drxmin = -5;
drxmax = 5;
drymin = -5;
drymax = 5;
drzmin = 0;
drzmax = 360;
S_drx = drxmin + (drxmax-drxmin)*rand(n_S,1);
S_dry = drymin + (drymax-drymin)*rand(n_S,1);
S_drz = drzmin + (drzmax-drzmin)*rand(n_S,1);
S_rx = S_rx + S_drx;
S_ry = S_ry + S_dry;
S_rz = S_rz + S_drz;

% Instantiate sensor cell arrays
S_R = cell(n_S,1);
S_t = cell(n_S,1);
for i = 1:n_S;
    S_R{i} = compose_rotation(pi*S_rx(i)./180,...
        pi*S_ry(i)./180,pi*S_rz(i)./180);
    S_t{i} = [S_tx(i) S_ty(i) S_tz(i)]';
end

x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
% Plot WCS points and sensors
%{
figure;
hold on
%plot3(0,0,0,'^k','markersize',10,...
%    'markerfacecolor',P_color(s,:));
%hdummy = plot3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),'ok','markersize',5,...
%    'markerfacecolor',P_color(s,:));
%set(hdummy, 'visible', 'off');
for t = 1:numel(T_x);
    h = filledCircle([T_x(t) T_y(t)],T_r(t),1000,[.5 .5 .5]);
end
rectangle('Position', [xmin,ymin, xmax-xmin, ymax-ymin], 'linestyle', '--');
axis equal
for i =1:n_S;
        x_axist = S_R{i}*x_axis + repmat(S_t{i},[1,2]);
        y_axist = S_R{i}*y_axis + repmat(S_t{i},[1,2]);
        z_axist = S_R{i}*z_axis + repmat(S_t{i},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'r', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
           'b', 'linewidth',2)
end
%}

%% Create LCS tie points

% Range Cutoff
cmin = 10;
cmax = 35;
S_c = cmin + (cmax-cmin)*rand(n_S,1);
smin = 15;
smax = 30;

% Random Subset
S_r = round(smin + (smax-smin)*rand(n_S,1));
P_WCS = cell(n_S,1);
P_rad = cell(n_S,1);
for i = 1:n_S;
    sensor_dist = sqrt((T_x - S_t{i}(1)).^2 + (T_y - S_t{i}(2)).^2 + (T_z - S_t{i}(3)).^2);
    is_valid = sensor_dist < S_c(i);
    P_WCS{i} = [T_x(is_valid) T_y(is_valid) T_z(is_valid)]';
    P_rad{i} = T_r(is_valid);
    is_valid = randperm(sum(is_valid),min([S_r(i) sum(is_valid)]));
    P_WCS{i} = P_WCS{i}(:,is_valid);
    P_rad{i} = P_rad{i}(is_valid);
end

% Transform by true rotation and translation
P_LCS = cell(n_S,1);
for i = 1:n_S;
    P_LCS{i} = S_R{i}'*P_WCS{i} - S_R{i}'*repmat(S_t{i}, [1, size(P_WCS{i},2)]);
end

% Add noise
nmin = .2;
nmax = .2;
S_n = nmin + (nmax-nmin)*rand(n_S,1);
P_LCSn = cell(n_S,1);
P_WCSn = cell(n_S,1);
P_RMSEin = zeros(n_S,1);
for i = 1:n_S;
    P_noise = S_n(i).*randn(3,size(P_LCS{i},2));
    P_LCSn{i} = P_LCS{i} + P_noise;
    P_WCSn{i} = P_WCS{i} + P_noise;
    P_RMSEin(i) = sqrt(nanmean(sum((P_LCSn{i} - P_LCS{i}).^2)));
end

% Modify diameters
ddmin = .2;
ddmax = .2;
S_n = ddmin + (ddmax-ddmin)*rand(n_S,1);
P_radn = cell(n_S,1);
for i = 1:n_S;
    P_noise = S_n(i).*randn(1,size(P_LCS{i},2))';
    P_radn{i} = P_rad{i}.*(1+ P_noise);
end

%% Plot noise-added tie points in WCS
%{
P_color = uint8(jet(n_S).*255);
P_color = P_color(randperm(n_S,n_S),:);
figure;
hold on
for t = 1:numel(T_x);
    h = filledCircle([T_x(t) T_y(t)],T_r(t),1000,[.5 .5 .5]);
end
rectangle('Position', [xmin,ymin, xmax-xmin, ymax-ymin], 'linestyle', '--');
axis equal
for i =1:n_S;
    for t = 1:size(P_LCS{i},2);
        h = filledCircle([P_WCSn{i}(1,t) P_WCSn{i}(2,t)],P_radn{i}(t),500,P_color(i,:));
        set(h, 'facealpha',0.5);
    end
    x_axist = S_R{i}*x_axis + repmat(S_t{i},[1,2]);
    y_axist = S_R{i}*y_axis + repmat(S_t{i},[1,2]);
    z_axist = S_R{i}*z_axis + repmat(S_t{i},[1,2]);
    plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
        P_color(i,:), 'linewidth',2);
    plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
        P_color(i,:), 'linewidth',2)
    plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
        P_color(i,:), 'linewidth',2)
    % rectangle('Position', [S_t{i}(1)-S_c(i),S_t{i}(2)-S_c(i), 2*S_c(i), 2*S_c(i)], ...
    %  'curvature', [1 1], 'linestyle', '-', 'edgecolor', P_color(i,:));
end
axis([xmin xmax ymin ymax]);
%}
%%
% Set up aux variables
aux.info_site = 1;
aux.info_valid_plot = cell(n_S,1);
for p = 1:n_S;
    aux.info_valid_plot{p} = sprintf('%02.0f', p);
end
aux.info_exp = 'GraphAnal';
aux.info_suffix = 'test';
aux.P_LCS = P_LCSn;
aux.P_rad = P_radn;
aux.path_mat = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix,'\');

upperdir = sprintf('%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\');
if ~exist(upperdir, 'dir');
    mkdir(upperdir)
end
trialdir = sprintf('%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix,'\');
if ~exist(trialdir, 'dir');
    mkdir(trialdir)
end
% Pairwise registrations
%kelbe_registration_combine_dis3fun_graphanalyses( aux )
kelbe_registration_combine_dijkstraposewcs_graphanalyses( aux )



