%% Toy Example
% Based on Harvey's suggestion

%% Create points
% Point parameters 
n_P = 20;   
P_xmin = -10;
P_xmax = 10;
P_ymin = -10;
P_ymax = 10;
P_zmin = -1;
P_zmax = 1;
P_rmin = .04;
P_rmax = .5;

% Random locations on uniform distribution 
Px = unifrnd(P_xmin, P_xmax, [n_P,1]);
Py = unifrnd(P_ymin, P_ymax, [n_P,1]);
Pz = unifrnd(P_zmin, P_zmax, [n_P,1]);
Pr = unifrnd(P_rmin, P_rmax, [n_P,1]);

%
figure;
scatter3(Px,Py,Pz,30,'k','filled');
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
title('Truth tree locations');
%legend('truth location','location','bestoutside');
%}

% Outputs 
% Px            x locations of truth points 
% Py            y locations of truth points 
% Pz            z locations of truth points 
% Pr            r of truth points 
clear P_rmin P_rmax
%% Create sensors

% Number of sensors 
n_S = 4;

% Sensor parameters 
T_xmin = -10;
T_xmax = 10;
T_ymin = -10;
T_ymax = 10;
T_zmin = -1;
T_zmax = 1;
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
scatter3(Px,Py,Pz,30,'k','filled');
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
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

% Number of points visible to each sensor 
V_min = 15;
V_max = 20;

% Noise in xyz location 
N_min = -.2;
N_max = .2;

for s = 1:n_S;
    % Find visible 
    n_visible = randi([V_min V_max],1);
    ix_rand = randperm(n_P)';
    is_visible = ix_rand(1:n_visible);
    % Add noise 
    Nx = unifrnd(N_min, N_max,[n_visible,1]);
    Ny = unifrnd(N_min, N_max,[n_visible,1]);
    Nz = unifrnd(N_min, N_max,[n_visible,1]);
    P_WCS = [Px(is_visible)'; Py(is_visible)'; Pz(is_visible)'];
    % Transform back to LCS 
    S{s}.P_LCS = (S{s}.R)'*(P_WCS - repmat(S{s}.t,[1,n_visible])) - [Nx Ny Nz]';
    S{s}.r = Pr(is_visible);
    S{s}.nP = n_visible;
end

%
% Individual camera views
for s = 1:n_S;
    figure;
    hold on
    plot3(S{s}.P_LCS(1,:),S{s}.P_LCS(2,:),S{s}.P_LCS(3,:),'ok','markersize',5,...
        'markerfacecolor',S{s}.color);
    axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    grid on
    titlestr = sprintf('Camera %g',s);
    title(titlestr);
    plot3(0,0,0,'^k','markersize',10,...
        'markerfacecolor',S{s}.color);
    legend_str{1} = 'Points in LCS';
    legend_str{2} = 'Camera';
end
%}
% Outputs 
% S.P_LCS               points in local coordinate system 
% S.r                   radius of points 
% S.nP                  number of points 
clear V_min V_max N_min N_max Nx Ny Nz P_WCS 
clear s n_visible ix_rand is_visible
%% Run registration code
match_i = cell(n_S); %base
match_j = cell(n_S); %mobile
match_R = cell(n_S);
match_t = cell(n_S);

for i = 1:n_S;
    for j = i:n_S;
        [ match_R{i,j}, match_t{i,j}, match_i{i,j}, match_j{i,j} ] = toy_registrationfunction(S{i}.P_LCS',S{j}.P_LCS',S{i}.r,S{j}.r);
    end
end

match_Pi = cell(n_S,n_S);
match_Pi_all = cell(n_S,n_S);
for i = 1:n_S;
    for j = i:n_S;
        match_Pi{i,j} = match_R{i,j}*S{j}.P_LCS(:,match_j{i,j}) + repmat(match_t{i,j},[1,numel(match_j{i,j})]);
        match_Pi_all{i,j} = match_R{i,j}*S{j}.P_LCS + repmat(match_t{i,j},[1,S{j}.nP]);
    end
end

%
% Transformed to camera 1
figure;
legend_str{1} = 'Truth points';
scatter3(Px,Py,Pz,30,'k','filled');
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
hold on
for s = 1:n_S;
    legend_str{s+1} = sprintf('Camera %g',s);
    plot3(S{s}.t(1),S{s}.t(2),S{s}.t(3),'^k','markersize',10,...
        'markerfacecolor',S{s}.color);
end
for s = 1:n_S;
    legend_str{n_S+1+s} = sprintf('Points %g in WCS',s);
    plot3(match_Pi{1,s}(1,:),match_Pi{1,s}(2,:),match_Pi{1,s}(3,:),'ok','markersize',5,...
        'markerfacecolor',S{s}.color);
end
legend(legend_str,'location','bestoutside');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
title('Points in WCS');
%}

% Outputs 
% match_i               index to points in i (world) 
% match_j               index to points in j (local) 
% match_R               (i,j) gives rotation matrix of j into i 
% match_t               (i,j) gives translation matrix of j into i 
% match_Pi              matched points of j transformed into CS of i 
% match_Pi_all          all points of j transformed into CS of i 
%% Find unique points

% Gather all points in WCS (i = 1) from each camera
all_Pi_x = [];
all_Pi_y = [];
all_Pi_z = [];
all_Pi_r = [];
for j = 1:n_S;
    all_Pi_x = [all_Pi_x match_Pi_all{1,j}(1,:)]; %#ok<*AGROW>
    all_Pi_y = [all_Pi_y match_Pi_all{1,j}(2,:)];
    all_Pi_z = [all_Pi_z match_Pi_all{1,j}(3,:)];
    all_Pi_r = [all_Pi_r S{j}.r'];
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

% Only interested in first set for now
D2_ii = false(size(D2_is));
ix = 1;
for s = 1:n_S;
    D2_ii(ix:ix+S{s}.nP-1,ix:ix+S{s}.nP-1) = true;
    ix = ix + S{s}.nP;
end
D2_is = D2_is&~D2_ii;

%{
figure;
imagesc(D2_is);
axis image
%}

is_unique = (sum(D2_is,1)==0);
unique_x = all_Pi_x(is_unique);
unique_y = all_Pi_y(is_unique);
unique_z = all_Pi_z(is_unique);
unique_r = all_Pi_r(is_unique);
n_unique = numel(unique_x);

%
figure;
hold on
scatter3(Px,Py,Pz,60,'k');
scatter3(all_Pi_x, all_Pi_y, all_Pi_z,10,all_Pi_r,'filled');
plot3(unique_x,unique_y,unique_z,'+k','markersize',10,...
        'markerfacecolor',[0 0 0], 'linewidth',1.5);
legend_str{1} = 'Truth points';
legend_str{2} = 'All points';
legend_str{3} = 'Unique points';
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
legend(legend_str,'location','bestoutside');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
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

data = nan(n_unique,n_S,3);

% Points in WCS 
for s = 1:n_S;
    S{s}.P_WCS = match_Pi_all{1,s};
end

for s = 1:n_S;
    % Matrices for distance calculations 
    A_unique_x = repmat(unique_x',[1,S{s}.nP]);
    A_unique_y = repmat(unique_y',[1,S{s}.nP]);
    A_unique_z = repmat(unique_z',[1,S{s}.nP]);
    A_unique_r = repmat(unique_r',[1,S{s}.nP]);
    A_Pi_x = repmat(S{s}.P_WCS(1,:),[n_unique,1]);
    A_Pi_y = repmat(S{s}.P_WCS(2,:),[n_unique,1]);
    A_Pi_z = repmat(S{s}.P_WCS(3,:),[n_unique,1]);
    A_Pi_r = repmat(S{s}.r',[n_unique,1]);
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
    
    data(u,s,1) = S{s}.P_LCS(1,c);
    data(u,s,2) = S{s}.P_LCS(2,c);
    data(u,s,3) = S{s}.P_LCS(3,c);
    
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
is_nan = any(isnan(data(:,:,1)),2);
n_isnnan = sum(~is_nan);
P0 = [P0 data(~is_nan,1,1)' data(~is_nan,1,2)' data(~is_nan,1,3)'];

% Levenberg-Marquardt (nested) 
[P,F] = LM_toy_nest(P0,data(~is_nan,:,:));

clear rx ry rz t s 
%% Reconstruct

for s = 1:n_S;
    S{s}.Ri = match_R{1,s};
    [S{s}.rxi,S{s}.ryi,S{s}.rzi] = decompose_rotation(S{s}.Ri);
    S{s}.ti = match_t{1,s};
end

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
    S{s}.Rf = compose_rotation(rx, ry, rz);
    S{s}.tf = P(ix+3:ix+5)';
    tf(s,:) = P(ix+3:ix+5)';
    ix = ix + 6;
end

datai = zeros(n_isnnan,n_S,3);
datai(:,1,:) = reshape(P0(ix:end),[(numel(P0)-ix+1)/3,3]);
for s = 2:n_S;
    datai(:,s,:) = (S{s}.Ri*squeeze(data(~is_nan,s,:))'+repmat(S{s}.ti,[1,n_isnnan]))';
end


dataf = zeros(n_isnnan,n_S,3);
dataf(:,1,:) = reshape(P(ix:end),[(numel(P)-ix+1)/3,3]);
for s = 2:n_S;
    dataf(:,s,:) = (S{s}.Rf*squeeze(data(~is_nan,s,:))'+repmat(S{s}.tf,[1,n_isnnan]))';
end

clear legend_str
figure;
legend_str{1} = 'Truth points';
scatter3(Px,Py,Pz,30,'k','filled');
xlabel('x');ylabel('y');zlabel('z');
view(0,90);
hold on
for s = 1:n_S
    plot3(datai(:,s,1),datai(:,s,2),datai(:,s,3),'ok','markersize',5,...
        'markerfacecolor','b');%S{s}.color);
    plot3(dataf(:,s,1),dataf(:,s,2),dataf(:,s,3),'ok','markersize',5,...
        'markerfacecolor','r');%S{s}.color);
    if s == 1;
        legend_str{2} = 'Initial Transformation';
        legend_str{3} = 'Nonlinear Optimization';
    end
end
legend(legend_str,'location','bestoutside');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
title('Points in WCS');
  

% Outputs 
% S.Ri                  initial rotation  
% S.ti                  initial translation 
% S.Rf                  final rotation after LM 
% S.tf                  final translation after LM 
% data_i                initial WCS 
% data_f                final WCS 
%   data                  LCS 
clear rx ry rz ix 

%% Error metrics 

% No a priori knowledge 
errori_p = sqrt(sum((datai - repmat(mean(datai,2),[1,n_S,1])).^2,3));
errorf_p = sqrt(sum((dataf - repmat(mean(dataf,2),[1,n_S,1])).^2,3));

errori_p_point = sum(errori_p,2);
errorf_p_point = sum(errorf_p,2);

errori_p_sensor = sum(errori_p,1);
errorf_p_sensor = sum(errorf_p,1);

erorri_p_all = sum(errori_p(:));
errorf_p_all = sum(errorf_p(:));

% A priori knowledge 

errori_t_sensor = sqrt(sum((ti - tt).^2,2));
errorf_t_sensor = sqrt(sum((tf - tt).^2,2));

errori_t_all = sum(errori_t_sensor);
errorf_t_all = sum(errorf_t_sensor);

errori_r_sensor = sqrt(sum((ri - rt).^2,2));
errorf_r_sensor = sqrt(sum((rf - rt).^2,2));

errori_r_all = sum(errori_r_sensor);
errorf_r_all = sum(errorf_r_sensor);





