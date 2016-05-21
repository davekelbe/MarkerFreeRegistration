%% Toy Example3
% Allows matches between 2,3 2,4 etc.
% Corresponds to LM_toy_nest3
% Based on Harvey's suggestion
%clear all; close all, clc;
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_points = false;
options_sensors = false;
options_imagepoints = false;
options_initialmatch = true;
options_unique = false;
%% Create points
fprintf('\nCreate points\n');
% Point parameters
n_P = 8;
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

% Redefine points if too close
P_minspacing = 2.5;
continuesearch = true;
while continuesearch
    dPx = repmat(Px,[1,n_P]) - repmat(Px',[n_P,1]);
    dPy = repmat(Py,[1,n_P]) - repmat(Py',[n_P,1]);
    dPz = repmat(Pz,[1,n_P]) - repmat(Pz',[n_P,1]);
    dP = sqrt(dPx.^2 + dPy.^2 + dPz.^2);
    is_reject = (dP<P_minspacing & triu(dP,1));
    [~,col] = find(is_reject);
    if isempty(col);
        continuesearch = false;
    end
    Px(col) =  unifrnd(P_xmin, P_xmax, [numel(col),1]);
    Py(col) =  unifrnd(P_ymin, P_ymax,  [numel(col),1]);
    Pz(col) =  unifrnd(P_zmin, P_zmax,  [numel(col),1]);
end

if options_verbose && options_points;
    figure;
    plot3(Px,Py,Pz,'+k','markersize',10,...
        'markerfacecolor',[0 0 0], 'linewidth',1.2);axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    title('Truth tree locations');
    grid on
end
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
n_S = 10;

% Sensor parameters
T_xmin = -9;
T_xmax = 9;
T_ymin = -9;
T_ymax = 9;
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

% Grid sensor pattern 
n_S = 9;
color = jet(n_S);
rt = zeros(n_S,3);
xv = -5:5:5;
yv = -5:5:5;
[x,y] = meshgrid(xv,yv);
tt = [x(:),y(:),zeros(n_S,1)];
tt(ceil(end/2),:) = tt(1,:);
tt(1,:) = [0 0 0];
%}

if options_verbose && options_sensors;
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
end
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
V_min = round(n_P*.5);
V_max = round(n_P*.8);

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
    P_rad{s} = Pr(is_visible);
    P_n(s) = n_visible;
end

% Individual camera views
if options_verbose && options_imagepoints;
    for s = 1:n_S;
        figure;
        hold on
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
        legend_str{1} = 'Points in LCS';
        legend_str{2} = 'Camera';
        legend(legend_str,'location','best');
    end
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
for i = 1:n_S;
    for j = i:n_S;
        fprintf('\n\tMatching %g to %g\n',j,i);
        [ a,b, c, d ] = toy_registrationfunction(P_LCS{i}',P_LCS{j}',P_rad{i},P_rad{j});
        if ~isempty(a)&&~isempty(b)&&~isempty(c)&&~isempty(d);
            match_R{i,j} = a;
            match_t{i,j} = b;
            match_i{i,j} = c;
            match_j{i,j} = d;
        end
    end
end

% Make R,t non-directed
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(match_R{i,j});
            match_R{j,i} = match_R{i,j}';
            match_t{j,i} = -(match_R{i,j}')*match_t{i,j};
        end
    end
end

%{
% Examine a match
i = 2;
j = 6;
figure;
hold on
% Source points
h = filledCircle([P_LCS{i}(1,1); P_LCS{i}(2,1)]',P_rad{i}(1),1000,'r');
h = filledCircle([P_LCS{j}(1,1); P_LCS{j}(2,1)]',P_rad{j}(1),1000,'b');
legend_str{1} = 'Local points i';
legend_str{2} = 'Local points j';
for t = 2:numel(P_rad{i});
h = filledCircle([P_LCS{i}(1,t); P_LCS{i}(2,t)]',P_rad{i}(t),1000,'r');
end
for t = 2:numel(P_rad{j});
h = filledCircle([P_LCS{j}(1,t); P_LCS{j}(2,t)]',P_rad{j}(t),1000,'b');
end
axis equal
grid on
legend(legend_str,'location','best');
title('Points input to matching');
for m = 1:numel(match_j{i,j});
    plot3([P_LCS{i}(1,match_i{i,j}(m)) P_LCS{j}(1,match_j{i,j}(m))],...
        [P_LCS{i}(2,match_i{i,j}(m)) P_LCS{j}(2,match_j{i,j}(m))],...
        [P_LCS{i}(3,match_i{i,j}(m)) P_LCS{j}(3,match_j{i,j}(m))],...
        '-k','linewidth',2);
end
%}

%{
% Examine a match
i = 8;
j = 6;
figure;
hold on
% Source points
h1 = scatter3(0,0,0,50,'k','filled');
h2 = scatter3(0,0,0,50,'k','^','filled');
set(h1,'visible','off');
set(h2,'visible','off');
legend_str{1} = 'Local points i';
legend_str{2} = 'Local points j';
for c = 1:size(color,1);
    h3 = scatter3(0,0,0,50,color(c,:),'s','filled');
    set(h3,'visible','off');
    legend_str{2+c} = sprintf('Unique point %g',c);
end
scatter3(P_LCS{i}(1,:),P_LCS{i}(2,:),P_LCS{i}(3,:),50,color(P_truthix{i},:),'filled');
scatter3(P_LCS{j}(1,:),P_LCS{j}(2,:),P_LCS{j}(3,:),50,color(P_truthix{j},:),'fill','^');
axis equal
grid on
legend(legend_str,'location','best');
title('Points input to matching');
for m = 1:numel(match_j{i,j});
    plot3([P_LCS{i}(1,match_i{i,j}(m)) P_LCS{j}(1,match_j{i,j}(m))],...
        [P_LCS{i}(2,match_i{i,j}(m)) P_LCS{j}(2,match_j{i,j}(m))],...
        [P_LCS{i}(3,match_i{i,j}(m)) P_LCS{j}(3,match_j{i,j}(m))],...
        '-k','linewidth',2);
end
%}
%% Effective R and t

% Shortest path
% Weighted adjacency matrix
G = cellfun(@(x) sqrt(sum(x.^2)), match_t);
% Make undirected
% G = G + G' - tdiag(diag(G));
G(logical(eye(size(G)))) = 0;
Gsparse = sparse(G);
G_path = cell(1,n_S);
for j = 2:n_S;
    [~,G_path{1,j},~] = graphshortestpath(Gsparse,j,1);
end

% Effective R and t
match_Reff = cell(n_S,n_S);
match_teff = cell(n_S,n_S);
for i = 1:n_S;
    for j = i+1:n_S;
        if ~isempty(match_R{i,j});
            path = G_path{1,j};
            Rtemp = eye(3);
            ttemp = zeros(3,1);
            for k = 1:numel(path)-1;
                if any(size(match_R{i,k})~=size(Rtemp));
                    foo = 1;
                end
                Rtemp = match_R{path(k+1),path(k)}*Rtemp;
                ttemp = match_R{path(k+1),path(k)}*(ttemp)+ match_t{path(k+1),path(k)};
            end
            match_Reff{i,j} = Rtemp;
            match_teff{i,j} = ttemp;
        end
    end
end

%{
% Check incremental R,t
i = 2;
j = 6;
path = G_path{1,j};
Rtemp = eye(3);
ttemp = zeros(3,1);
k = 1;
Rtemp = match_R{path(k+1),path(k)}*Rtemp;
ttemp = match_R{path(k+1),path(k)}*(ttemp)+ match_t{path(k+1),path(k)};
%}

% Not necessary
% Make Reff,teff non-directed
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(match_Reff{i,j});
            match_Reff{j,i} = match_Reff{i,j}';
            match_teff{j,i} = -(match_Reff{i,j}')*match_teff{i,j};
        end
    end
end
match_Reff{1,1} = eye(3);
match_teff{1,1} = zeros(3,1);


% Find points j transformed into new coordinate system i
match_Pi = cell(n_S,n_S);           % Points that are matched
match_Pi_all = cell(n_S,n_S);       % All points
match_rad_all = cell(n_S,n_S);
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(match_Reff{i,j});
            match_Pi{i,j} = match_Reff{i,j}*P_LCS{j}(:,match_j{i,j}) + repmat(match_teff{i,j},[1,numel(match_j{i,j})]);
            match_Pi_all{i,j} = match_Reff{i,j}*P_LCS{j} + repmat(match_teff{i,j},[1,P_n(j)]);
            match_rad_all{i,j} = P_rad{j};%(match_j{i,j});
        end
    end
end

% Examine a pair
%{
   figure;
    legend_str{1} = 'Truth points';
    plot3(Px,Py,Pz,'+k','markersize',10,...
        'markerfacecolor',[0 0 0], 'linewidth',1.2);xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    hold on
    ix = 1;
    for i  =1:n_S;
        for j = 1:n_S;
            if ~isempty(match_Pi{i,j});
                % legend_str{n_S+1+ix} = sprintf('Points %g in WCS',j);
                hij = plot3(match_Pi{i,j}(1,:),match_Pi{i,j}(2,:),match_Pi{i,j}(3,:),'ok','markersize',5,...
                    'markerfacecolor',color(i,:));
                ix = ix + 1;
                set(hij,'visible','off');
            end
        end
    end
    legend(legend_str,'location','best');
    axis equal
    axis(1.5*[P_xmin P_xmax P_ymin P_ymax -10 10]);
    title('Points in WCS - Initial Transformation');
    grid on
%}

%
% Transformed to camera 1
if options_verbose && options_initialmatch;
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
    ix = 1;
    for i  =1:n_S;
        for j = 1:n_S;
            if ~isempty(match_Pi{i,j});
                % legend_str{n_S+1+ix} = sprintf('Points %g in WCS',j);
                hij = plot3(match_Pi{i,j}(1,:),match_Pi{i,j}(2,:),match_Pi{i,j}(3,:),'ok','markersize',5,...
                    'markerfacecolor',color(i,:));
                ix = ix + 1;
                %set(hij,'visible','off');
            end
        end
    end
    legend(legend_str,'location','best');
    axis equal
    axis(1.5*[P_xmin P_xmax P_ymin P_ymax -10 10]);
    title('Points in WCS - Initial Transformation');
    grid on
end
%}

% Outputs
% match_i               index to points in i (world)
% match_j               index to points in j (local)
% match_R               (i,j) gives rotation matrix of j into i
% match_t               (i,j) gives translation matrix of j into i
% match_Pi              matched points of j transformed into CS of 1 *NEW
% match_Pi_all          all points of j transformed into CS of 1 *NEW
%% Find unique points
fprintf('\nFind unique points\n');

% Gather all points in WCS (i = 1) from each camera
all_Pi_x = [];
all_Pi_y = [];
all_Pi_z = [];
all_Pi_r = [];
for i = 1:n_S;
    for j = 1:n_S;
        if ~isempty(match_Pi_all{i,j});
            all_Pi_x = [all_Pi_x match_Pi_all{i,j}(1,:)]; %#ok<*AGROW>
            all_Pi_y = [all_Pi_y match_Pi_all{i,j}(2,:)];
            all_Pi_z = [all_Pi_z match_Pi_all{i,j}(3,:)];
            all_Pi_r = [all_Pi_r match_rad_all{i,j}']; %match_R_all or P_rad is the same
        end
    end
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
t_xyz = 2.^2;
t_r = 0.2;
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
if any(size(D2_is)~= size(D2_ii));
    foo = 1;
end
D2_is = D2_is&~D2_ii;

%{
figure;
imagesc(D2_ii);
axis image
%}

% Find unique points
is_unique = (sum(D2_is,1)==0);
unique_x = all_Pi_x(is_unique);
unique_y = all_Pi_y(is_unique);
unique_z = all_Pi_z(is_unique);
unique_r = all_Pi_r(is_unique);
unique_ix = find(is_unique);
n_unique = numel(unique_x);
all_unique = nan(n_all,1);
D2_is = D2_is | eye(n_all,n_all);
for i = 1:n_unique;
    all_unique(D2_is(unique_ix(i),:)) = i;
end
    

%
if options_verbose && options_unique;
    figure;
    hold on
    plot3(Px,Py,Pz,'+k','markersize',15,...
        'markerfacecolor',[0 0 0], 'linewidth',1.4);
    scatter3(all_Pi_x, all_Pi_y, all_Pi_z,30,all_unique,'filled','markeredgecolor','k');
    scatter3(unique_x,unique_y,unique_z,180,'k');%,'markersize',20,...
    %  'markerfacecolor',[1 1 1],'alpha',1, 'linewidth',.5);
    legend_str{1} = 'Truth points';
    legend_str{2} = 'All points';
    legend_str{3} = 'Unique points';
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    legend(legend_str,'location','best');
    axis equal
    axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
    title('Unique points ');
end
%}

% Output
% unique_x                  x values of unique points
% unique_y                  y values of unique points
% unique_z                  z values of unique points
% unique_r                  r values of unique points
% n_unique                  number of unique points
%clear all_Pi_x all_Pi_y all_Pi_z all_Pi_r
clear A_Pi_xt A_Pi_yt A_Pi_zt A_Pi_rt
clear A_Pi_x A_Pi_y A_Pi_z A_Pi_r
clear D2_xyz D2_r D2_isxyz D2_isr D2_is
clear is_unique
%% Build up data array
fprintf('\nBuild up data array\n');

data_LCSi = nan(n_unique,n_S,3);
data_truthix = nan(n_unique,n_S);

% Points in WCS
is_match = ~cellfun(@isempty,match_Pi_all);
P_WCS = cell(n_S,1);
for s = 1:n_S;
    ix_match = find(is_match(:,s));
    if isempty(ix_match);
        continue
    end
    ix1 = ix_match(1);
    P_WCS{s} = match_Pi_all{ix1,s};
end

for s = 1:n_S;
    if isempty(P_WCS{s});
        continue
    end
    % Matrices for distance calculations
    % Unique xyz values
    A_unique_x = repmat(unique_x',[1,P_n(s)]);
    A_unique_y = repmat(unique_y',[1,P_n(s)]);
    A_unique_z = repmat(unique_z',[1,P_n(s)]);
    A_unique_r = repmat(unique_r',[1,P_n(s)]);
    %
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
    % data_LCSi is n_unique x n_S x 3
    data_LCSi(u,s,1) = P_LCS{s}(1,c);
    data_LCSi(u,s,2) = P_LCS{s}(2,c);
    data_LCSi(u,s,3) = P_LCS{s}(3,c);
    data_truthix(u,s) = P_truthix{s}(c); % Should be identical for all sensors
end

% True WCS
data_WCSt = nan(n_unique,n_S,3);
for s = 1:n_S;
    temp_ix =  data_truthix(:,s);
    temp_isnan = isnan(temp_ix);
    data_WCSt(temp_isnan,s,:) = NaN;
    data_WCSt(~temp_isnan,s,:) = [Px(temp_ix(~temp_isnan)) Py(temp_ix(~temp_isnan)) Pz(temp_ix(~temp_isnan))];
end

% True LCS
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
%% Find graph paths

% Selects paths used for effective R,t
t_max_cxn_to1 = 5;
t_max_cxndist_to1 = 20;

% Find combinations and permutations
G_comb = cell(n_S,1);
G_perm = cell(n_S,1);
n_comb = zeros(n_S,1);
n_perm = zeros(n_S,1);
for s = 2:n_S;
    G_comb{s} = combnk(1:n_S,s);
    n_comb(s) = size(G_comb{s},1);
end
for s = 2:t_max_cxn_to1;
    all_perm = perms(1:s);
    G_perm{s} = all_perm(1:end/2,:);
    n_perm(s) = size(G_perm{s},1);
end
n_cp = n_comb.*n_perm;

G_cp = cell(n_S,1);
G_t = cell(n_S,1);
for s = 2:t_max_cxn_to1;
    G_cp{s} = nan(n_cp(s),s);
    G_t{s} = nan(n_cp(s),s,3);
    ix = 1;
    for c = 1:n_comb(s);
        for p = 1:n_perm(s);
            G_cp{s}(ix,:) = G_comb{s}(c,G_perm{s}(p,:));
            G_t{s}(ix,:,:) = tt(G_comb{s}(c,G_perm{s}(p,:)),:);
            ix = ix + 1;
        end
    end
end

% Find distance
G_dist = cell(n_S,1);
for s = 2:t_max_cxn_to1;
    dist = circshift(G_t{s},[0,0,0]) - circshift(G_t{s},[0,-1,0]);
    dist = dist(:,1:end-1,:);
    dist = sqrt(sum(sum(dist.^2,3),2));
    G_dist{s} = dist;
end

% Remove paths which are too far
G_valid = cell(n_S,1);
for s = 1:n_S;
    G_valid{s} = (G_dist{s}<t_max_cxndist_to1);
end
G_path = cell(n_S,1);
n_path = zeros(n_S,1);
for s = 1:n_S;
    G_path{s} = G_cp{s}(G_valid{s},:);
    G_t{s} = G_t{s}(G_valid{s},:,:);
    G_dist{s} = G_dist{s}(G_valid{s});
    n_path(s) = size(G_path{s},1);
end

% Remove paths which are not connected
G_valid = cell(n_S,1);
% Make non-directed
match_empty = triu(cellfun(@isempty,match_i),1);
match_empty = or(match_empty, match_empty');
for g = 2:t_max_cxn_to1;
    is_reject = false(n_path(g),1);
    for p = 1:n_path(g);
        path = G_path{g}(p,:);
        for k = 1:g-1;
            if match_empty(path(k+1),path(k));
                is_reject(p) = true;
            end
        end
    end
    G_valid{g} = ~is_reject;
end
for s = 1:n_S;
    G_path{s} = G_path{s}(G_valid{s},:);
    G_t{s} = G_t{s}(G_valid{s},:,:);
    G_dist{s} = G_dist{s}(G_valid{s});
    n_path(s) = size(G_path{s},1);
end

all_dist = [];
all_path = nan(sum(n_path),t_max_cxn_to1);
all_endstart = nan(sum(n_path),1);
for g = 2:t_max_cxn_to1;
    all_dist = [all_dist; G_dist{g}];
    all_path(sum(n_path(2:g-1))+1:sum(n_path(2:g)),1:g) = G_path{g};
    all_endstart(sum(n_path(2:g-1))+1:sum(n_path(2:g)),1) = G_path{g}(:,1);
    all_endstart(sum(n_path(2:g-1))+1:sum(n_path(2:g)),2) = G_path{g}(:,end);
end

% best path to WCS structure (i==1)
best_path = cell(n_S,1);
for j = 1:n_S;
    ixs = (all_endstart(:,2) == 1);
    ixe = (all_endstart(:,1) == j);
    ix_valid = find(ixs & ixe);
    if isempty(ix_valid);
        continue
    end
    [~,ix_min] = min(all_dist(ix_valid));
    is_nnan = ~isnan(all_path(ix_valid(ix_min),:));
    best_path{j} = all_path(ix_valid(ix_min),is_nnan);
end

% R,t from best path
best_Rieff = cell(n_S,1);
best_tieff = cell(n_S,1);
for j = 1:n_S;
    path = best_path{j};
    if isempty(path);
        continue
    end
    Reff = eye(3);
    teff = zeros(3,1);
    for k = 1:numel(path)-1;
        Reff = match_R{path(k+1),path(k)}*Reff;
        teff = match_R{path(k+1),path(k)}*(teff)+ match_t{path(k+1),path(k)};
    end
    best_Rieff{j} = Reff;
    best_tieff{j} = teff;
end
best_Rieff{1} = eye(3);
best_tieff{1} = zeros(3,1);

%{
% make non-directed
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(data_Rieff{i,j});
            data_Rieff{j,i} = data_Rieff{i,j}';
            data_tieff{j,i} = -(data_Rieff{i,j}')*data_tieff{i,j};
        end
    end
end
foo = 1;
%}


%{
for s = 2:t_max_cxn;
    clear legend_str
    color_cxn = jet(n_path(s));
    figure;
    plot3(Px,Py,Pz,'+k','markersize',10,...
        'markerfacecolor',[0 0 0], 'linewidth',1.2);axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    hold on
    legend_str{1} = 'Truth tree locations';
    for c = 1:n_S;
        legend_str{c+1} = sprintf('Camera %g',c);
        plot3(tt(c,1),tt(c,2),tt(c,3),'^k','markersize',10,...
            'markerfacecolor',color(c,:));
    end
    for p = 1:n_path(s);
         plot3(G_t{s}(p,:,1),G_t{s}(p,:,2),G_t{s}(p,:,3),'-.','color',color_cxn(p,:), 'linewidth',2);
    end
    legend_str{c+2} = sprintf('Paths of %g',s);
    title('Points and camera positions');
    legend(legend_str,'location','best');
    grid on
end
%}


%% Levenberg Marquardt
fprintf('\nLevenberg Marquardt\n');

% Selects paths used in LM
t_max_cxn = 3;
t_max_cxndist = 11;
%
% Remove paths which are too far
G_valid = cell(n_S,1);
for s = 1:t_max_cxn;
    G_valid{s} = (G_dist{s}<t_max_cxndist);
end
%G_path = cell(n_S,1);
n_path = zeros(n_S,1);
for s = 1:n_S;
    G_path{s} = G_path{s}(G_valid{s},:);
    G_t{s} = G_t{s}(G_valid{s},:,:);
    G_dist{s} = G_dist{s}(G_valid{s});
    n_path(s) = size(G_path{s},1);
end
is_valid = (n_path>0);
for s = 1:n_S;
    if ~is_valid(s);
        G_path{s} = [];
        G_t{s} = [];
        G_dist{s} = [];
    end
end
%}

% Add transformation parameters to parameter vector
P0 = [];
ri = nan(n_S,3);
ti = nan(n_S,3);
%rri = zeros(n_S,n_S,3);
%tti = zeros(n_S,n_S,3);
%foo = zeros(n_S,n_S);
ctr = 1;
R_empty = cellfun(@isempty,match_R);
for i = 1:n_S-1;
    for j = i+1:n_S;
        if ~R_empty(i,j);
            %foo(i,j) = ctr;
            ctr = ctr  + 1;
            [rx,ry,rz] = decompose_rotation(match_R{i,j});
            t = match_t{i,j};
            if i == 1;
                ri(j,1) = rx;
                ri(j,2) = ry;
                ri(j,3) = rz;
                ti(j,:) = t;
            end
            P0 = [P0 rx ry rz t(1) t(2) t(3)];
        end
    end
end

% Add position 1 data (exclude entries with NaN's)
is_nan = (sum(~isnan(data_LCSi(:,:,1)),2)<=1);
n_isnnan = sum(~is_nan);


data_P1 = nan(n_unique,3);
data_Six = zeros(n_unique,1);
for s = 1:n_S;
    is_valid = isnan(data_P1(:,1));
    data_P1(is_valid,:) = data_LCSi(is_valid,s,:);
    data_Six(is_valid,:) = s;
end

ix = numel(P0)+1;
P0 = [P0 data_P1(~is_nan,1)' data_P1(~is_nan,2)' data_P1(~is_nan,3)'];

%
% Initial WCS
data_WCSi = nan(n_unique,n_S,3);

for s = 1:n_S;
    temp_isnotnan = ~isnan(data_LCSi(:,s));
    temp_nisnotnan = sum(temp_isnotnan);
    data_WCSi(temp_isnotnan,s,:) = (best_Rieff{s}*squeeze(data_LCSi(temp_isnotnan,s,:))'+...
        repmat(best_tieff{s},[1,temp_nisnotnan]))';
end

%{
% Initial WCS - paths
clear legend_str
figure;
legend_str{1} = 'Truth points';
plot3(Px,Py,Pz,'+k','markersize',10,...
    'markerfacecolor',[0 0 0], 'linewidth',1.2);xlabel('x');ylabel('y');zlabel('z');
view(0,90);
hold on
for s = 1:n_S
    plot3(data_WCSi(:,s,1),data_WCSi(:,s,2),data_WCSi(:,s,3),'ok','markersize',5,...
        'markerfacecolor',color(s,:));
    legend_str{s+1} = sprintf('Camera %g',s);
end
legend(legend_str,'location','best');
axis equal
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
title('Points in WCS before LM');
grid on
%}

% Levenberg-Marquardt (nested)
options = optimset('lsqnonlin');%'Display','iter');
options.Display = 'iter-detailed';
options.MaxIter = 15;
[P,~,~] = LM_toy_nest5(P0,data_LCSi(~is_nan,:,:),data_Six(~is_nan),R_empty,G_path, options );

clear rx ry rz t s
%% Reconstruct

rf = zeros(n_S,3);
tf = zeros(n_S,3);
%rri = zeros(n_S,n_S,3);
%tti = zeros(n_S,n_S,3);
%foo = zeros(n_S,n_S);
ctr = 1;
ix = 1;
match_Rf = cell(n_S,n_S);
match_tf = cell(n_S,n_S);
for i = 1:n_S-1;
    for j = i+1:n_S;
        if ~R_empty(i,j);
            foo(i,j) = ctr;
            ctr = ctr  + 1;
            rx = P(ix);
            ry = P(ix+1);
            rz = P(ix+2);
            t = P(ix+3:ix+5)';
            match_Rf{i,j} = compose_rotation(rx,ry,rz);
            match_tf{i,j} = t;
            ix = ix + 6;
            if i == 1;
                rf(j,1) = rx;
                rf(j,2) = ry;
                rf(j,3) = rz;
                tf(j,:) = t;
            end
        end
    end
end

% Make match_Rf non-directed
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(match_Rf{i,j});
            match_Rf{j,i} = match_Rf{i,j}';
            match_tf{j,i} = -(match_Rf{i,j}')*match_tf{i,j};
        end
    end
end
match_Rf{1,1} = eye(3);
match_tf{1,1} = zeros(3,1);

% R,t from best path
best_Rfeff = cell(n_S,1);
best_tfeff = cell(n_S,1);
for j = 1:n_S;
    path = best_path{j};
    if isempty(path);
        continue
    end
    Reff = eye(3);
    teff = zeros(3,1);
    for k = 1:numel(path)-1;
        Reff = match_Rf{path(k+1),path(k)}*Reff;
        teff = match_Rf{path(k+1),path(k)}*(teff)+ match_tf{path(k+1),path(k)};
    end
    best_Rfeff{j} = Reff;
    best_tfeff{j} = teff;
end
best_Rfeff{1} = eye(3);
best_tfeff{1} = zeros(3,1);

%% Final WCS
% Final WCS by camera
data_WCSf = nan(n_unique,n_S,3);
for s = 1:n_S
    temp_isnotnan = ~isnan(data_LCSi(:,s));
    temp_nisnotnan = sum(temp_isnotnan);
    data_WCSf(temp_isnotnan,s,:) = (best_Rfeff{s}*squeeze(data_LCSi(temp_isnotnan,s,:))'+...
        repmat(best_tfeff{s},[1,temp_nisnotnan]))';
end

% Final and Initial WCS
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
axis(1.5*[P_xmin P_xmax P_ymin P_ymax P_zmin P_zmax]);
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
%{
% Points in LCS after LM optimization
data_LCSf = nan(n_unique,n_S,3);
is_S1 = (data_Six == 1);
temp_S1 = reshape(P(ix:end),[(numel(P)-ix+1)/3,3]);
data_LCSf(is_S1,1,:) = temp_S1(is_S1,:);
for s = 2:n_S;
    R = compose_rotation(rf(s,1),rf(s,2),rf(s,3));
    T = tf(s,:)';
    temp_isnotnan = ~isnan(data_LCSi(:,s));
    temp_nisnotnan = sum(temp_isnotnan);
    data_LCSf(temp_isnotnan,s,:) = (R'*squeeze(data_WCSf(temp_isnotnan,s,:))'-repmat(T,[1,temp_nisnotnan]))';
end

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
sqerrori_p = nansum((data_WCSi - data_WCSt).^2,3);
sqerrorf_p = nansum((data_WCSf - data_WCSt).^2,3);

rmsei_p_point = sqrt(nanmean(sqerrori_p,2));
rmsef_p_point = sqrt(nanmean(sqerrorf_p,2));

rmsei_p_sensor = sqrt(nanmean(sqerrori_p,1));
rmsef_p_sensor = sqrt(nanmean(sqerrorf_p,1));

rmsei_p_all = sqrt(nanmean(sqerrori_p(:)));
rmsef_p_all = sqrt(nanmean(sqerrorf_p(:)));

rmsei_t_sensor = sqrt(nanmean((ti - tt).^2,2));
rmsef_t_sensor = sqrt(nanmean((tf - tt).^2,2));

rmsei_t_all = nanmean(rmsei_t_sensor);
rmsef_t_all = nanmean(rmsef_t_sensor);

rmsei_r_sensor = rad2deg(sqrt(nanmean((ri - rt).^2,2)));
rmsef_r_sensor = rad2deg(sqrt(nanmean((rf - rt).^2,2)));

rmsei_r_all = nanmean(rmsei_r_sensor);
rmsef_r_all = nanmean(rmsef_r_sensor);

fprintf('\nAll errors\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
fprintf('%3.3f \t%3.3f \tRMSE between points [meters]\n', rmsei_p_all, rmsef_p_all);
fprintf('%3.3f \t%3.3f \tRMSE of rotation comp. to truth [degrees]\n', rmsei_r_all, rmsef_r_all);
fprintf('%3.3f \t%3.3f \tRMSE of translation comp. to truth [m]\n', rmsei_t_all, rmsef_t_all);

fprintf('\nErrors by sensor\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for s = 1:n_S;
    fprintf('%3.3f \t%3.3f \tRMSE for sensor %g\n', rmsei_p_sensor(s), rmsef_p_sensor(s), s);
end


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

fprintf('\nErrors by point\n');
fprintf('***************\n');
fprintf('Before \tAfter \tValue\n');
for p = 1:n_isnnan;
    fprintf('%3.3f \t%3.3f \tRMSE between points %g\n', rmsei_p_point(p), rmsef_p_point(p), p);
end
foo =1;

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


