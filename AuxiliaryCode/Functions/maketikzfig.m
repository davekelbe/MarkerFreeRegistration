%% Make tikz figure

% Make data 
n_tree = 50;
seed = 5;
rng(seed);
x_min = -40;
x_max = 40;
tree_x = (x_max - x_min)*rand(n_tree,1) + x_min;
seed = 6;
rng(seed);
y_min = -40;
y_max = 40;
tree_y = (y_max - y_min)*rand(n_tree,1) + y_min;

d_min = .1;
d_max = 1;
seed = 9;
    rng(seed);
    tree_d = 1.5*(d_max - d_min)*rand(n_tree,1) + d_min;
    
    %{
    figure
    hold on
    for t = 1:n_tree;
        str = sprintf('%s %d','\leftarrow',t);
        h = filledCircle([tree_x(t); tree_y(t)]',tree_d(t),500,'r');
        text(tree_x(t), tree_y(t),str);
    end
    axis auto
    axisval = axis;
    xlabel('x Position relative to plot center [m]');
    ylabel('y Position relative to plot center [m]');
    grid on
    title(sprintf('seed = %d',seed));
%}
tree_id = 1:n_tree;
tree_x(8)  = 10;
tree_d(8) = 1;
tree_y(46)  = -20;
tree_d(46) = 1;
tree_d(6) = 1;
tree_y(19)  = 42;
tree_y(17)  = 31;
tree_x(12)  = 25;
tree_y(12)  = 25;
tree_x(13)  = -12;
tree_y(13)  = 4;
tree_x(8)  = 15;
tree_y(8)  = -5;
tree_x(43)  = 35;
tree_y(43)  = 6;
ix_remove = [41 24 1 14 9 11 5 4 40 22 34 28 34 5 29 50 10 42 36 46 15];
is_remove = false(n_tree,1);
is_remove(ix_remove) = true;
tree_x = tree_x(~is_remove);
tree_y = tree_y(~is_remove);
tree_d = tree_d(~is_remove);
tree_id = tree_id(~is_remove);
n_tree = numel(tree_x);

    %{
    figure
    hold on
    for t = 1:n_tree;
        str = sprintf('%s %d','\leftarrow',t);
        h = filledCircle([tree_x(t); tree_y(t)]',tree_d(t),500,'r');
        text(tree_x(t), tree_y(t),str);
    end
    axis auto
    axisval = axis;
    xlabel('x Position relative to plot center [m]');
    ylabel('y Position relative to plot center [m]');
    grid on
   % title(sprintf('seed = %d',seed));
%}

clear ans x_max x_min y_min y_max seed ix_remove is_remove d_min d_max
    % Outputs
    % tree_x                x values of tie points 
    % tree_y                y values of tie points 
    % tree_d                radius values of tie points 
%{
figure
hold on
for t = 1:n_tree;
    str = sprintf('%s %d','\leftarrow',t);
    h = filledCircle([tree_x(t); tree_y(t)]',tree_d(t),500,'r');
    text(tree_x(t), tree_y(t),str);
end
axis auto
axisval = axis;
xlabel('x Position relative to plot center [m]');
ylabel('y Position relative to plot center [m]');
grid on
%}
%% Make two scans 

S_center = [-10 -10 0; 10 10 0];
n_S = size(S_center,1);
radius = 30;
is_scan = false(n_S,n_tree);

P_WCS = cell(n_S,1);
P_rad = cell(n_S,1);
P_id = cell(n_S,1);
for s = 1:n_S;
    tree_dist = sqrt((tree_x - S_center(s,1)).^2 + (tree_y - S_center(s,2)).^2);
    is_scan(s,:) = tree_dist<radius;  
    P_WCS{s} = [tree_x(is_scan(s,:))'; tree_y(is_scan(s,:))'; zeros(1,sum(is_scan(s,:)))];
    P_rad{s} = tree_d(is_scan(s,:));
    P_id{s} = tree_id(is_scan(s,:));
end
n_tree = cellfun(@(x) size(x,2),P_WCS);

    t_rad = 0.2; %0.06^2; % Trees can be matched if their radius is within threshold
    t_coll = 0.1; % Collinnearity threshold: don't consider triangles that are collinear
    t_eig_error = inf;%1e1;
    t_RANSAC_xyz = 0.4^2;
    t_RANSAC_nsearch = 16;
    P_color = [1 0 0; 0 0 1];

clear tree_dist t s h 
    % Outputs
    % S_center              scan position
    % tree_y                y values of tie points 
    % tree_d                radius values of tie points
%% Place in LCS 
rx = [0 0]; 
ry = [0 0];
%rz = [0 rad2deg(-90)];
rz = [0 deg2rad(90)];
tx = [S_center(1,1) S_center(2,1) ]; 
ty = [S_center(1,2) S_center(2,2) ];
tz = [S_center(1,3) S_center(2,3) ];

Rtrue = cell(n_S,1);
ttrue = cell(n_S,1);
for s = 1:n_S;
Rtrue{s} = compose_rotation(rx(s),ry(s),rz(s));
ttrue{s} = [tx(s) ty(s) tz(s)]';
end

P_LCS = cell(n_S,1);
for s = 1:n_S;
    P_LCS{s} = Rtrue{s}'*(P_WCS{s} - repmat(ttrue{s}, [1,n_tree(s)]));
end

%{
s = 1;
figure; 
scatter(P_WCS{s}(1,:), P_WCS{s}(2,:),100, 'r','filled');
hold on;
s=2;
scatter(P_WCS{s}(1,:), P_WCS{s}(2,:),60, 'b','filled');
s = 1;
figure; 
scatter(P_LCS{s}(1,:), P_LCS{s}(2,:),100, 'r','filled');
hold on;
s=2;
scatter(P_LCS{s}(1,:), P_LCS{s}(2,:),60, 'b','filled');
%}

%% Rearrange data to match registration function
aux.P_LCS = P_LCS;
aux.P_color = P_color;
%aux.P_plot = P_plot;
aux.P_rad = P_rad;
%aux.filepath_ply = filepath_ply;
%aux.info_exp = info_exp;
%aux.info_plot = info_plot;
%aux.info_site = info_site;
%aux.info_slash = info_slash;
%aux.info_suffix = info_suffix;
%aux.info_valid_plot = info_valid_plot;
aux.n_S = n_S;
aux.n_tree = n_tree;
%aux.options_axesval = options_axesval;
%aux.options_imagepoints = options_imagepoints;
%aux.options_initialmatch = options_initialmatch;
%aux.options_loadmatch = options_loadmatch;
%aux.options_loadvar = options_loadvar;
%aux.options_unique = options_unique;
%aux.options_verbose = options_verbose;
%aux.p_std = p_std;
%aux.path_mat = path_mat;
%aux.path_ply = path_ply;
%aux.path_save = path_save;
%aux.path_site = path_site;
%aux.path_tikz = path_tikz;
aux.t_RANSAC_nsearch = t_RANSAC_nsearch;
%aux.t_RANSAC_rad = t_RANSAC_rad;
aux.t_RANSAC_xyz = t_RANSAC_xyz;
aux.t_coll = t_coll;
aux.t_eig_error = t_eig_error;
aux.t_rad = t_rad;

aux.S_center = S_center;
aux.p_radius = radius;
aux.P_id = P_id;
aux.Rtrue = Rtrue;
aux.ttrue = ttrue;
aux.P_WCS = P_WCS;
aux.rztrue = rz;
[match_R, match_t] = kelbe_registration_function_tikzfig( aux );





