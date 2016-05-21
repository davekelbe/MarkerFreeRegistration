function [  ] = graph_analysis_rmsefun_generic_param_analyses( aux )
% Limits of full tie point set
aux.path_mat = sprintf('%s%s%s%s%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix,'\',aux.trial,'\');
filepath_S_R = sprintf('%s%s',aux.path_mat, 'S_R.mat');

    % Sensors
    n_Sg = 5;
    n_S = round(n_Sg.^2);
    
if ~exist(filepath_S_R, 'file');
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
    [close1, ~] = find(is_close);
    T_x(close1) = nan;
    T_y(close1) = nan;
    T_z(close1) = nan;
    T_r(close1) = nan;
    T_x = T_x(~isnan(T_x));
    T_y = T_y(~isnan(T_y));
    T_z = T_z(~isnan(T_z));
    T_r = T_r(~isnan(T_r));
    n_T = numel(T_x);
    
    % Nominal position
    S_tx = repmat((-20:10:20),[n_Sg,1]);
    S_ty = repmat((-20:10:20)',[1,n_Sg]);
    S_tx = S_tx(:);
    S_ty = S_ty(:);
    S_tz = zeros(n_S,1);
    
    % Deviation around nominal position
    % was 1 and .2
    dxmin = -3; 
    dxmax = 3;
    dymin = -3;
    dymax = 3;
    dzmin = -.2;
    dzmax = .2;
    %}
    %{
    dxmin = 0;
    dxmax = 0;
    dymin = 0;
    dymax = 0;
    dzmin = 0;
    dzmax = 0;
    %}
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
    %
    drxmin = -10;
    drxmax = 10;
    drymin = -10;
    drymax = 10;
    drzmin = -10;
    drzmax = 10;
    %}
    %{
    drxmin = 0;
    drxmax = 0;
    drymin = 0;
    drymax = 0;
    drzmin = 0;
    drzmax = 0;
    %}
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
    cmin = 15;
    cmax = 30;
    S_c = cmin + (cmax-cmin)*rand(n_S,1);
    
    smin = 20;
    smax = 30;
    % Random Subset
    S_r = round(smin + (smax-smin)*rand(n_S,1));
    P_WCS = cell(n_S,1);
    P_rad = cell(n_S,1);
    P_ix = cell(n_S,1);
    for i = 1:n_S;
        P_ix{i} = 1:n_T;
        sensor_dist = sqrt((T_x - S_t{i}(1)).^2 + (T_y - S_t{i}(2)).^2 + (T_z - S_t{i}(3)).^2);
        is_valid = sensor_dist < S_c(i);
        P_WCS{i} = [T_x(is_valid) T_y(is_valid) T_z(is_valid)]';
        P_rad{i} = T_r(is_valid);
        P_ix{i} = P_ix{i}(is_valid);
        is_valid = randperm(sum(is_valid),min([S_r(i) sum(is_valid)]));
        P_WCS{i} = P_WCS{i}(:,is_valid);
        P_rad{i} = P_rad{i}(is_valid);
        P_ix{i} = P_ix{i}(is_valid);
    end
    
    % Transform by true rotation and translation
    P_LCS = cell(n_S,1);
    for i = 1:n_S;
        P_LCS{i} = S_R{i}'*P_WCS{i} - S_R{i}'*repmat(S_t{i}, [1, size(P_WCS{i},2)]);
    end
    
    % Add noise
    %n_mean = aux.n_mean;
    %n_std = aux.n_std;
    S_n = aux.S_n;%normrnd(n_mean, n_std, [n_S,1]);
    %S_n = nmin + (nmax-nmin)*rand(n_S,1);
    % Uniform
    %{
nmin = aux.nmin;
nmax = aux.nmax;
S_n = nmin + (nmax-nmin)*rand(n_S,1);
    %}
    P_LCSn = cell(n_S,1);
    P_WCSn = cell(n_S,1);
    P_RMSEin = zeros(n_S,1);
    for i = 1:n_S;
        P_noise = S_n.*randn(3,size(P_LCS{i},2));
        P_LCSn{i} = P_LCS{i} + P_noise;
        P_WCSn{i} = P_WCS{i} + P_noise;
        P_RMSEin(i) = sqrt(nanmean(sum((P_LCSn{i} - P_LCS{i}).^2)));
    end
    
    % Modify diameters
    %
    ddmin = .05;
    ddmax = .15;
    %}
    %{
    ddmin = .0;
    ddmax = .0;
    %}
    %ddmean = .06;
    %ddsigma =
    %S_n = normrnd(
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
    %aux.info_exp = 'GraphAnal';
    %aux.info_suffix = 'test';
    %aux.trial = 1;
    aux.P_LCSn = P_LCSn;
    aux.P_radn = P_radn;
    upperdir = sprintf('%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\');
    if ~exist(upperdir, 'dir');
        mkdir(upperdir)
    end
    suffixdir = sprintf('%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix,'\');
    if ~exist(suffixdir, 'dir');
        mkdir(suffixdir)
    end
    trialdir = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix,'\', aux.trial, '\');
    if ~exist(trialdir, 'dir');
        mkdir(trialdir)
    end
    
    %{
filepath_P_LCSn = sprintf('%s%s',aux.path_mat, 'P_LCSn.mat');
filepath_P_radn = sprintf('%s%s',aux.path_mat, 'P_radn.mat');
if exist(filepath_P_LCSn, 'file');
    load(filepath_P_LCSn);
    load(filepath_P_radn)
    aux.P_LCS = P_LCSn;
    aux.P_rad = P_radn;
end
    %}
end
% Set up aux variables
aux.info_site = 1;
aux.info_valid_plot = cell(n_S,1);
for p = 1:n_S;
    aux.info_valid_plot{p} = sprintf('%02.0f', p);
end

%{
filepath_P_LCSn = sprintf('%s%s',aux.path_mat, 'P_LCSn.mat');
load(filepath_P_LCSn);
aux.P_LCS = P_LCSn;
filepath_P_radn = sprintf('%s%s',aux.path_mat, 'P_radn.mat');
load(filepath_P_radn);
aux.P_rad = P_radn;
%}

% Pairwise registrations
%kelbe_registration_combine_dis3fun_graphanalyses( aux )
kelbe_registration_combine_dijkstraposewcs_graphanalyses( aux )

% Save sensor poses
filepath_S_R = sprintf('%s%s',aux.path_mat, 'S_R.mat');
if ~exist(filepath_S_R, 'file');
    save(filepath_S_R, 'S_R');
    filepath_S_t = sprintf('%s%s',aux.path_mat, 'S_t.mat');
    save(filepath_S_t, 'S_t');
    % Save truth locations and diameters
    filepath_T_x = sprintf('%s%s',aux.path_mat, 'T_x.mat');
    save(filepath_T_x, 'T_x');
    filepath_T_y = sprintf('%s%s',aux.path_mat, 'T_y.mat');
    save(filepath_T_y, 'T_y');
    filepath_T_z = sprintf('%s%s',aux.path_mat, 'T_z.mat');
    save(filepath_T_z, 'T_z');
    filepath_T_r = sprintf('%s%s',aux.path_mat, 'T_r.mat');
    save(filepath_T_r, 'T_r');
    % Save tie points
    filepath_P_LCSn = sprintf('%s%s',aux.path_mat, 'P_LCSn.mat');
    save(filepath_P_LCSn, 'P_LCSn');
    filepath_P_LCS = sprintf('%s%s',aux.path_mat, 'P_LCS.mat');
    save(filepath_P_LCS, 'P_LCS');
    filepath_P_radn = sprintf('%s%s',aux.path_mat, 'P_radn.mat');
    save(filepath_P_radn, 'P_radn');
    filepath_P_ix = sprintf('%s%s',aux.path_mat, 'P_ix.mat');
    save(filepath_P_ix, 'P_ix');
    filepath_P_WCSn = sprintf('%s%s',aux.path_mat, 'P_WCSn.mat');
    save(filepath_P_WCSn, 'P_WCSn');
    filepath_P_WCS = sprintf('%s%s',aux.path_mat, 'P_WCS.mat');
    save(filepath_P_WCS, 'P_WCS');
    % Save RMSE_in
    filepath_P_RMSEin = sprintf('%s%s',aux.path_mat, 'P_RMSEin.mat');
    save(filepath_P_RMSEin, 'P_RMSEin');
end
end
