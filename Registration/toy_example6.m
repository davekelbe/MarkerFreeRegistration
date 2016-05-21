%% Toy Example6 Real data
% Allows matches between 2,3 2,4 etc.
% Corresponds to LM_toy_nest3
%clear all; close all, clc;
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = false;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
% Initialize program-level options and paths
% Outputs
% options               options
% paths                 paths
clear D ctr d filepath_tree plot t tree
%% Load points
fprintf('\nLoad points\n');

info_exp = 'Harvard';
info_suffix = '03-01';
info_slash = '\'; 
info_site = 31;
path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
D = dir(path_site); 
ctr = 1;
info_valid_plot = {'08','11','13','15','18'};
for d = 1:numel(D);
    if strcmp(D(d).name(1), '.')
        continue
    end
    info_plot = D(d).name;
    if ~any(strcmp(info_plot, info_valid_plot))
        continue
    end
    path_mat = sprintf('%s%s%smat%s',path_site,info_plot,info_slash,info_slash);
    filepath_tree = sprintf('%stree.mat',path_mat);
    if ~exist(filepath_tree,'file')
        continue
    end
    load(filepath_tree);
    n_tree = numel(tree);
    plot = str2num(info_plot);
    P_LCS{ctr} = nan(3,n_tree);
    P_rad{ctr} = nan(n_tree,1);
    P_plot(ctr) = plot;
    for t = 1:n_tree;
        P_LCS{ctr}(:,t) = tree(t).loc(:,1);
        P_rad{ctr}(t) = tree(t).r(1);
    end
    path_ply{ctr} = sprintf('%s%s%sply%s',path_site,info_plot,info_slash,info_slash);
    filepath_ply{ctr} = sprintf('%spoints_full_%03.0f-%02.0f.ply', path_ply{ctr}, info_site, plot);
    ctr = ctr + 1;
end
P_n = cellfun(@numel,P_rad);
n_S = numel(P_LCS);

i_xmin = min(cellfun(@(x) min(x(1,:)),P_LCS));
i_xmax = max(cellfun(@(x) max(x(1,:)),P_LCS));
i_ymin = min(cellfun(@(x) min(x(2,:)),P_LCS));
i_ymax = max(cellfun(@(x) max(x(2,:)),P_LCS));
i_zmin = min(cellfun(@(x) min(x(3,:)),P_LCS));
i_zmax = max(cellfun(@(x) max(x(3,:)),P_LCS));

% Colormap for sensors
color = jet(n_S);

% Individual camera views
if false;%options_verbose && options_imagepoints;
    for s = 1:n_S;
        clear legend_str
        figure
        hold on
        plot3(0,0,0,'^k','markersize',10,...
            'markerfacecolor',color(s,:));
        hdummy = plot3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),'ok','markersize',5,...
            'markerfacecolor',color(s,:));
        set(hdummy, 'visible', 'off');
        for t = 1:numel(P_rad{s});
            h = filledCircle([P_LCS{s}(1,t); P_LCS{s}(2,t)]',P_rad{s}(t),1000,color(s,:));
        end
        %scatter3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),30,...
        %    color_P_index(truth_P_index{s},:),'filled');
        %axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
        axis auto
        axisval = axis;
        %set(gca, 'xtick',
        xlabel('x');ylabel('y');zlabel('z');
        view(0,90);
        grid on
        %titlestr = sprintf('Scan %g',P_plot(s));
        %title(titlestr);
        legend_str{1} = sprintf('Scanner %g',P_plot(s));
        legend_str{2} = sprintf('Stem map in LCS %g', P_plot(s));
        legend(legend_str,'location','northeast');
        legend boxoff       % Hides the legend's axes (legend border and background)
        set(gca, 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        filepath_save = sprintf('%sLCS_%02.0f.eps',path_save, P_plot(s));
        saveas(gcf,filepath_save,'psc2')
    end
end
%}
% Outputs
% P_LCS                 points in local coordinate system (cell)
% P_rad                 radius of points
% P_plot                plot number
% P_n                   number of points
% i_xmin                minimum x value of data
% i_zmax                maximum z value of data
% color                 colormap for scans
% info                  program level information
% options               program level options
% path                  paths
% n                     number of given variable
clear D ctr d filepath_tree plot t tree
%% Run registration code
% temporarily saved variables for efficient coding and debugging
filepath_match_R = 'D:\Users\djk2312\Documents\MATLAB\temp\match_R.mat';
filepath_match_t = 'D:\Users\djk2312\Documents\MATLAB\temp\match_t.mat';
filepath_match_i = 'D:\Users\djk2312\Documents\MATLAB\temp\match_i.mat';
filepath_match_j = 'D:\Users\djk2312\Documents\MATLAB\temp\match_j.mat';

if false;% options_loadmatch == false;
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
    save(filepath_match_R, 'match_R');
    save(filepath_match_t, 'match_t');
    save(filepath_match_i, 'match_i');
    save(filepath_match_j, 'match_j');
else
    load(filepath_match_R);
    load(filepath_match_t);
    load(filepath_match_i);
    load(filepath_match_j);
end
for i = 1:n_S;
    match_R{i,i} = eye(3);
    match_t{i,i} = zeros(3,1);
    match_i{i,i} = 1:numel(P_LCS{i});
    match_j{i,i} = 1:numel(P_LCS{i});
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
clear legend_str
i = 4;
j = 5;
clear legend_str
figure;
hold on
% Source points
h = scatter(P_LCS{i}(1,1), P_LCS{i}(2,1),'ok','markerfacecolor','r');
set(h,'visible', 'off');
h = scatter(P_LCS{j}(1,1), P_LCS{j}(2,1),'ok','markerfacecolor','b');
set(h,'visible', 'off');
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
% Outputs
% match_R               (i,j) is pairwise rotation from j into i
% match_t               (i,j) is pairwise translation from j into i
% match_i               (i,j) is pairwise matches from i
% match_j               (i,j) is pairwise matches from j

clear i j h 
%% Remove bad matches
t_n_minmatch = 4;
t_rmin = 20;

match_nmatch = cellfun(@numel,match_i);
match_nmatch = match_nmatch + tril(match_nmatch.',-1);
match_rx = 180*cellfun(@decompose_rotation_rx,match_R)/pi;
match_ry = 180*cellfun(@decompose_rotation_ry,match_R)/pi;
match_rz = 180*cellfun(@decompose_rotation_rz,match_R)/pi;

is_valid = (match_nmatch > t_n_minmatch) & ...
    (match_rx < t_rmin) & (match_ry < t_rmin) & (match_rz < t_rmin);
% is valid should be symmetric
is_valid = (is_valid & is_valid');
match_R(~is_valid) = {[]};
match_t(~is_valid) = {[]};
match_i(~is_valid) = {[]};
match_j(~is_valid) = {[]};

% Plot paths
%{
P_plot_str = cell(n_S,1);
for s = 1:n_S;
    P_plot_str{s} = sprintf('Plot %g',P_plot(s));
end
view(biograph(triu(is_valid,1),P_plot_str,'ShowArrows','off' ));
%}

% Outputs
% is_valid              adjacency matrix for scan connections
% t_n_minmatch          minimum number of matches required
% t_rmin                minimum radius for a match
% match_nmatch          number of matches for a given link 
clear i j match_rx match_ry match_rz
%% Effective (WCS) R and t

% Shortest path to WCS
% Weighted undirected adjacency matrix
G = cellfun(@(x) sqrt(sum(x.^2)), match_t);
G(logical(eye(size(G)))) = 0;
Gsparse = sparse(G);
G_path = cell(1,n_S);
for j = 2:n_S;
    [~,G_path{1,j},~] = graphshortestpath(Gsparse,j,1);
end
G_path{1,1} = 1;

% Effective R and t back to WCS
match_Reff = cell(1,n_S);
match_teff = cell(1,n_S);
i = 1;
for j = i:n_S;
    path = G_path{1,j};
    Rtemp = eye(3);
    ttemp = zeros(3,1);
    for k = 1:numel(path)-1;
        if any(size(match_R{i,k})~=size(Rtemp));
            break
        end
        Rtemp = match_R{path(k+1),path(k)}*Rtemp;
        ttemp = match_R{path(k+1),path(k)}*(ttemp)+ match_t{path(k+1),path(k)};
    end
    match_Reff{i,j} = Rtemp;
    match_teff{i,j} = ttemp;
end

% All points transformed to WCS
match_Pi_all = cell(1,n_S);
match_rad_all = cell(1,n_S);
for j = 1:n_S;
    if ~isempty(match_Reff{1,j});
        match_Pi_all{1,j} = match_Reff{1,j}*P_LCS{j} + repmat(match_teff{1,j},[1,P_n(j)]);
        match_rad_all{1,j} = P_rad{j};
    end
end

% Visualization
if options_verbose && options_initialmatch;
    figure;
    view(0,90);
    hold on
    ix = 1;
    %i = 1; % Match_Pi is the points which match after effective RT to WCS
    for j = 1:n_S;
        if ~isempty(match_Pi_all{1,j});
            legend_str{ix} = sprintf('Points %g in WCS',j);
            hij = plot3(match_Pi_all{1,j}(1,:),match_Pi_all{1,j}(2,:),match_Pi_all{1,j}(3,:),'ok','markersize',5,...
                'markerfacecolor',color(j,:));
            ix = ix + 1;
            %set(hij,'visible','off');
        end
    end
    legend(legend_str,'location','best');
    axis equal
    axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
    title('Points in WCS - Initial Transformation');
    grid on
end
%}

%{
% Initial transformation of PLY
match_teff_arr = zeros(n_S,3);
for s = 1:n_S;
    match_teff_arr(s,:) = match_teff{s}';
end
for j = 1:n_S;
    fprintf('\nWriting ply %g of %g\n',j,n_S);
    [vertex, ~] = read_ply(filepath_ply{j});
    data_x2 = vertex(:,1);
    data_y2 = vertex(:,2);
    data_z2 = vertex(:,3);
    color_ply = vertex(:,4:6);
    filepath_ply_reg{j} = sprintf('%spoints_full_%03.0f-%02.0f-%02.0f.ply', ...
        path_ply{j}, info_site, P_plot(j),P_plot(1));
    xyz2t = (match_Reff{1,j}*[data_x2 data_y2 data_z2]') + ...
        repmat(match_teff{1,j},1,numel(data_x2));
    write2ply(filepath_ply_reg{j},xyz2t', color);
end
clear match_teff_arr vertex data_x2 data_y2 data_z2 color_ply xyz2t
%}

G_path_circ = cell(1,n_S);
for s = 1:n_S;
G_path_circ{s} = [fliplr(G_path{s}) (G_path{s}(2:end))];
end


% Outputs
% match_Pi_all          all points of j transformed into WCS
% match_rad_all         radii of points in WCS
% G_path                (j) best path from j into WCS
% G_path_circ           (j) circular path from j into WCS (for validation)
% match_Reff            effective rotation from j into WCS
% match_teff            effective translation from j into WCS
clear G Gsparse Rtemp i j k path s ttemp
%% Repeat registration with disjoint subsets

filepath_match_R_sub1 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_R_sub1.mat';
filepath_match_t_sub1 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_t_sub1.mat';
filepath_match_i_sub1 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_i_sub1.mat';
filepath_match_j_sub1 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_j_sub1.mat';
filepath_match_i_ix1 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_i_ix1.mat';
filepath_match_R_sub2 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_R_sub2.mat';
filepath_match_t_sub2 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_t_sub2.mat';
filepath_match_i_sub2 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_i_sub2.mat';
filepath_match_j_sub2 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_j_sub2.mat';
filepath_match_i_ix2 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_i_ix2.mat';
filepath_match_all_sub1 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_all_sub1.mat';
filepath_match_all_sub2 = 'D:\Users\djk2312\Documents\MATLAB\temp\match_all_sub2.mat';

if false; %options_loadmatch == false;
    fprintf('\nRun registration code\n');
    match_i_sub1 = cell(n_S); %base
    match_j_sub1 = cell(n_S); %mobile
    match_R_sub1 = cell(n_S);
    match_t_sub1 = cell(n_S);
    match_i_ix1 = cell(n_S);
    match_i_sub2 = cell(n_S); %base
    match_j_sub2 = cell(n_S); %mobile
    match_R_sub2 = cell(n_S);
    match_t_sub2 = cell(n_S);
    match_i_ix2 = cell(n_S);
    match_all_sub1 = cell(n_S);
    match_all_sub2 = cell(n_S);
  
    for i = 1:n_S;
        for j = 1:n_S;
            if i==j;
                continue
            end
            fprintf('\n\tMatching %g to %g\n',j,i);
            ix_randi = randperm(P_n(i));
            mid = floor(P_n(i)/2);
            ix_sub1 = ix_randi(1:mid);
            ix_sub2 = ix_randi(mid + 1:end);
            match_i_ix1{i,j} = ix_sub1;
            match_i_ix2{i,j} = ix_sub2;
            [ a,b, c, d ] = toy_registrationfunction(P_LCS{i}(:,ix_sub1)',P_LCS{j}',P_rad{i}(ix_sub1),P_rad{j});
            [ e,f, g, h ] = toy_registrationfunction(P_LCS{i}(:,ix_sub2)',P_LCS{j}',P_rad{i}(ix_sub2),P_rad{j});
            if ~isempty(a)&&~isempty(b)&&~isempty(c)&&~isempty(d);
                match_R_sub1{i,j} = a;
                match_t_sub1{i,j} = b;
                match_i_sub1{i,j} = ix_sub1(c);
                match_j_sub1{i,j} = d;
                match_all_sub1{i,j} = ix_sub1;
            end
            if ~isempty(e)&&~isempty(f)&&~isempty(g)&&~isempty(h);
                match_R_sub2{i,j} = e;
                match_t_sub2{i,j} = f;
                match_i_sub2{i,j} = ix_sub2(g);
                match_j_sub2{i,j} = h;
                match_all_sub2{i,j} = ix_sub2;
            end
        end
    end
    save(filepath_match_R_sub1, 'match_R_sub1');
    save(filepath_match_t_sub1, 'match_t_sub1');
    save(filepath_match_i_sub1, 'match_i_sub1');
    save(filepath_match_j_sub1, 'match_j_sub1');
    save(filepath_match_i_ix1, 'match_i_ix1');
    save(filepath_match_R_sub2, 'match_R_sub2');
    save(filepath_match_t_sub2, 'match_t_sub2');
    save(filepath_match_i_sub2, 'match_i_sub2');
    save(filepath_match_j_sub2, 'match_j_sub2');
    save(filepath_match_i_ix2, 'match_i_ix2');
    save(filepath_match_all_sub1, 'match_all_sub1');
    save(filepath_match_all_sub2, 'match_all_sub2');    
else
    load(filepath_match_R_sub1);
    load(filepath_match_t_sub1);
    load(filepath_match_i_sub1);
    load(filepath_match_j_sub1);
    load(filepath_match_i_ix1);
    load(filepath_match_R_sub2);
    load(filepath_match_t_sub2);
    load(filepath_match_i_sub2);
    load(filepath_match_j_sub2);
    load(filepath_match_i_ix2);
    load(filepath_match_all_sub1);
    load(filepath_match_all_sub1);
end
for i = 1:n_S;
    match_R_sub1{i,i} = eye(3);
    match_t_sub1{i,i} = zeros(3,1);
    match_i_sub1{i,i} = 1:numel(P_LCS{i});
    match_j_sub1{i,i} = 1:numel(P_LCS{i});
    match_R_sub2{i,i} = eye(3);
    match_t_sub2{i,i} = zeros(3,1);
    match_i_sub2{i,i} = 1:numel(P_LCS{i});
    match_j_sub2{i,i} = 1:numel(P_LCS{i});
end

% Remove bad matches 
t_n_minmatch = 4;
t_rmin = 20;

match_nmatch_sub1 = cellfun(@numel,match_i_sub1);
match_nmatch_sub2 = cellfun(@numel,match_i_sub2);

[match_rx_sub1, match_ry_sub1, match_rz_sub1] = cellfun(@(x) decompose_rotation(x),match_R_sub1);
[match_rx_sub2, match_ry_sub2, match_rz_sub2] = cellfun(@(x) decompose_rotation(x),match_R_sub2);
match_rx_sub1 = rad2deg(match_rx_sub1);
match_ry_sub1 = rad2deg(match_ry_sub1);
match_rz_sub1 = rad2deg(match_rz_sub1);
match_rx_sub2 = rad2deg(match_rx_sub2);
match_ry_sub2 = rad2deg(match_ry_sub2);
match_rz_sub2 = rad2deg(match_rz_sub2);

% Valid i,j pairs are those which both subsets are valid 
is_valid_sub1 = (match_nmatch_sub1 > t_n_minmatch) & ...
    (abs(match_rx_sub1) < t_rmin) & (abs(match_ry_sub1) < t_rmin) & (abs(match_rz_sub1) < t_rmin);
is_valid_sub2 = (match_nmatch_sub2 > t_n_minmatch) & ...
    (abs(match_rx_sub2) < t_rmin) & (abs(match_ry_sub2) < t_rmin) & (abs(match_rz_sub2) < t_rmin);
is_valid_sub = ~((~is_valid_sub1) |(~is_valid_sub2));

match_R_sub1(~is_valid_sub) = {[]};
match_t_sub1(~is_valid_sub) = {[]};
match_i_sub1(~is_valid_sub) = {[]};
match_j_sub1(~is_valid_sub) = {[]};
match_R_sub2(~is_valid_sub) = {[]};
match_t_sub2(~is_valid_sub) = {[]};
match_i_sub2(~is_valid_sub) = {[]};
match_j_sub2(~is_valid_sub) = {[]};

% Outputs
% filepath_match_*_sub  filepath of matlab data related to subsets 
% match_*_sub           pairwise matching rotation, translation, indices
% is_valid_sub          adjacency matrix for subset scan connections (directed)
% match_nmatch_sub      number of matches for a link (subset)
clear G Gsparse Rtemp i j k path s ttemp is_valid_sub1 is_valid_sub2
%% Make match_R 4x4 

match_R4 = cell(n_S);
match_R4_sub1 = cell(n_S);
match_R4_sub2 = cell(n_S);
for i = 1:n_S;
    for j = 1:n_S;
        if is_valid(i,j);
            match_R4{i,j} = [[match_R{i,j}; 0 0 0],[match_t{i,j}; 1]];
        end
        if is_valid_sub(i,j);
            match_R4_sub1{i,j} = [[match_R_sub1{i,j}; 0 0 0],[match_t_sub1{i,j}; 1]];
            match_R4_sub2{i,j} = [[match_R_sub2{i,j}; 0 0 0],[match_t_sub2{i,j}; 1]];
        end
    end
end

% Outputs 
% match_R4              pairwise 4x4 [R|t]
% match_R4_sub1         pairwise 4x4 [R|t] for sub1 
% match_R4_sub2         pairwise 4x4 [R|t] for sub2 
clear s i j 
%% Compare R and t from subsets and full set 

% Compare R and t from subsets to full set
match_rx = nan(n_S);
match_ry = nan(n_S);
match_rz = nan(n_S);
match_tx = nan(n_S);
match_ty = nan(n_S);
match_tz = nan(n_S);
%{
match_rx_sub1 = nan(n_S);
match_ry_sub1 = nan(n_S);
match_rz_sub1 = nan(n_S);
match_tx_sub1 = nan(n_S);
match_ty_sub1 = nan(n_S);
match_tz_sub1 = nan(n_S);
match_rx_sub2 = nan(n_S);
match_ry_sub2 = nan(n_S);
match_rz_sub2 = nan(n_S);
match_tx_sub2 = nan(n_S);
match_ty_sub2 = nan(n_S);
match_tz_sub2 = nan(n_S);
%}

for i = 1:n_S;
    for j = 1:n_S;
        if ~isempty(match_R{i,j});
            [match_rx(i,j), match_ry(i,j), match_rz(i,j)] = decompose_rotation(match_R{i,j});
            match_tx(i,j) = match_t{i,j}(1);
            match_ty(i,j) = match_t{i,j}(2);
            match_tz(i,j) = match_t{i,j}(3);
        end
        %{
        if ~isempty(match_R_sub1{i,j});
            [match_rx_sub1(i,j), match_ry_sub1(i,j), match_rz_sub1(i,j)] = decompose_rotation(match_R_sub1{i,j});
            match_tx_sub1(i,j) = match_t_sub1{i,j}(1);
            match_ty_sub1(i,j) = match_t_sub1{i,j}(2);
            match_tz_sub1(i,j) = match_t_sub1{i,j}(3);
        end
        if ~isempty(match_R_sub2{i,j});
            [match_rx_sub2(i,j), match_ry_sub2(i,j), match_rz_sub2(i,j)] = decompose_rotation(match_R_sub2{i,j});
            match_tx_sub2(i,j) = match_t_sub2{i,j}(1);
            match_ty_sub2(i,j) = match_t_sub2{i,j}(2);
            match_tz_sub2(i,j) = match_t_sub2{i,j}(3);
        end
        %}
    end
end

%{
% Examine a match
for i = 1:n_S;
    for j = 1:n_S;
        if i==j;
            continue;
        end
        clear legend_str
        figure;
        hold on
        % Source points
        h = plot(P_LCS{i}(1,1),P_LCS{i}(2,1),'ok','markersize',5,'markerfacecolor','r');
        set(h,'visible', 'off');
        h = plot(P_LCS{i}(1,1),P_LCS{i}(2,1),'ok','markersize',5,'markerfacecolor','g');
        set(h,'visible', 'off');
        h = plot(P_LCS{i}(1,1),P_LCS{i}(2,1),'ok','markersize',5,'markerfacecolor','b');
        set(h,'visible', 'off');
        legend_str{1} = 'Local points i subset 1';
        legend_str{2} = 'Local points i subset 2';
        legend_str{3} = 'Local points j';
        ix1 = match_i_ix1{i,j};
        ix2 = match_i_ix2{i,j};
        for t = 1:numel(ix1);
            h = filledCircle([P_LCS{i}(1,ix1(t)); P_LCS{i}(2,ix1(t))]',...
                P_rad{i}(ix1(t)),1000,'r');
        end
        for t = 1:numel(ix2);
            h = filledCircle([P_LCS{i}(1,ix2(t)); P_LCS{i}(2,ix2(t))]',...
                P_rad{i}(ix2(t)),1000,'g');
        end
        for t = 1:numel(P_rad{j});
            h = filledCircle([P_LCS{j}(1,t); P_LCS{j}(2,t)]',P_rad{j}(t),1000,'b');
        end
        axis equal
        grid on
        legend(legend_str,'location','best');
        titlestr = sprintf('Camera %g into %g',j,i);
        title(titlestr);
        for m = 1:numel(match_j{i,j}); % Matches
            plot3([P_LCS{i}(1,match_i{i,j}(m)) P_LCS{j}(1,match_j{i,j}(m))],...
                [P_LCS{i}(2,match_i{i,j}(m)) P_LCS{j}(2,match_j{i,j}(m))],...
                [P_LCS{i}(3,match_i{i,j}(m)) P_LCS{j}(3,match_j{i,j}(m))],...
                '-k','linewidth',1.5);
        end
        for m = 1:numel(match_j_sub1{i,j});
            plot3([P_LCS{i}(1,match_i_sub1{i,j}(m)) P_LCS{j}(1,match_j_sub1{i,j}(m))],...
                [P_LCS{i}(2,match_i_sub1{i,j}(m)) P_LCS{j}(2,match_j_sub1{i,j}(m))],...
                [P_LCS{i}(3,match_i_sub1{i,j}(m)) P_LCS{j}(3,match_j_sub1{i,j}(m))],...
                '-.r','linewidth',2);
        end
        for m = 1:numel(match_j_sub2{i,j});
            plot3([P_LCS{i}(1,match_i_sub2{i,j}(m)) P_LCS{j}(1,match_j_sub2{i,j}(m))],...
                [P_LCS{i}(2,match_i_sub2{i,j}(m)) P_LCS{j}(2,match_j_sub2{i,j}(m))],...
                [P_LCS{i}(3,match_i_sub2{i,j}(m)) P_LCS{j}(3,match_j_sub2{i,j}(m))],...
                '-.g','linewidth',2);
        end
        foo = 1;
    end
end 
%}
% Outputs
% no outputs 
clear match_rx match_ry match_rz match_tx match_ty match_tz
clear match_rx_sub1 match_ry_sub1 match_rz_sub1
clear match_tx_sub1 match_ty_sub1 match_tz_sub1
clear match_rx_sub2 match_ry_sub2 match_rz_sub2
clear match_tx_sub2 match_ty_sub2 match_rz_sub2
%% Find all possible permutations

% Selects paths used for effective R,t
t_max_cxn_to1 = n_S;
%t_max_cxndist_to1 = inf;

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
    G_perm{s} = all_perm; %[all_perm all_perm(:,1)]; % append to make circular
    n_perm(s) = size(G_perm{s},1);
end
n_cp = n_comb.*n_perm;

n_all = sum(n_comb.*n_perm);
G_all = nan(n_all,n_S);
ix = cumsum(n_comb.*n_perm);
for s = 2:n_S;
    comb_rep = repmat(G_comb{s},[n_perm(s),1]);
    perm_rep = zeros(n_cp(s),s);
    for p = 1:n_perm(s);
        perm_rep(p*n_comb(s)-n_comb(s) + 1:(p+1)*n_comb(s)-n_comb(s),:) = ...
            repmat(G_perm{s}(p,:), [n_comb(s),1]);
    end
    idx = sub2ind([n_cp(s),s],repmat((1:n_cp(s))',[1,s]),perm_rep);
    temp_all = comb_rep(idx);
    G_all(ix(s-1)+1:ix(s),1:s-1) = temp_all(:,1:s-1);
    G_all(ix(s-1)+1:ix(s),end) = temp_all(:,end);
end

% Output 
% G_all                 cell array of all paths (based on perm and comb)
% n_all                 number of all paths
clear G_comb G_perm n_comb n_perm s n_cp ix idx temp_all comb_rep perm_rep p
%% Make circular paths   

G_val = [];
G_i = [];
G_j = [];
for i = 1:n_S;
    for j = 1:n_S;
        is_i = (G_all(:,1)==i);
        is_j = (G_all(:,end) == j);
        is = is_i & is_j;
        V_temp = G_all(is,:);
        n_vtemp =size(V_temp,1);
        for v = 1:n_vtemp;
        G_val = [G_val; [ V_temp repmat(fliplr(V_temp(v,1:end-1)),n_vtemp,1)]];
        G_j = [G_j; repmat(i,[n_vtemp,1])];
        G_i = [G_i; repmat(j,[n_vtemp,1])];
        end
    end
end
n_v = size(G_val,1);

% Output 
% G_val                 array of circular validation paths ('NaN' in empty)
% G_i                   first node of path 
% G_j                   last node of path (should == G_i)
% n_v                   number of all circular validation paths
clear is_i is_j i j V_temp n_vtemp
%% Create cell array and determine additional parameters 
G_npath = sum(~isnan(G_val),2);
G_path = cell(n_v,1);
G_edge = cell(n_v,1);
G_duplicate = cell(n_v,1);
G_order = cell(n_v,1);

for v = 1:n_v;
    path = G_val(v,~isnan(G_val(v,:)));
    edges = [path(1:end-1); circshift(path(1:end-1),[1,-1])];
    [edges_sort, Isort] = sort(edges,1);
    [~,~,ixdup] = unique(edges_sort','rows','stable');
    ixdup = ixdup';
    duplicate = zeros(1,G_npath(v)-1);
    order = zeros(1,G_npath(v)-1);
    for d = 1:n_S;
        isdup = ixdup==d;
        if sum(isdup)>1;
            duplicate(isdup) = d;
            ix = find(isdup);
            order(ix(1)) = 1;
            order(ix(2)) = 2;
        end
    end
    G_path{v} = path;
    G_edge{v} = edges;
    G_duplicate{v} = duplicate;
    G_order{v} = order;
end 

match_teff_arr = zeros(n_S,3);
for s = 1:n_S;
    match_teff_arr(s,:) = match_teff{s}'; % From full match
end

G_dist = zeros(n_v,1);
for v = 1:n_v;
   loc1 = match_teff_arr(G_edge{v}(1,:),:);
   loc2 = match_teff_arr(G_edge{v}(2,:),:);
   G_dist(v) = sum(sqrt(sum((loc1 - loc2).^2,2))); 
end

% Output 
% G_path                cell array of path nodes  
% G_npath               number of nodes in each path
% G_edge                cell array of 2x#nodes edges (for computation)
% G_duplicate           logical where edge traverse is repeated 
% G_order               order in which repeated edge is discovered
% G_j                   last node of path (should == G_i)
% n_v                   number of all circular validation paths
clear v path edges edges_sort Isort ixdup duplicate order d isdup ix 
clear match_teff_arr s loc1 loc2 

%% Traverse path and (i) compare R,t 
G_R4 = cell(n_v,1);
G_nmatch = cell(n_v,1);
G_valid = false(n_v,1);
for v = 1:n_v;
    R4temp = eye(4);
    G_nmatch{v} = zeros(G_npath(v)-1,1);
    bad_edge = false;
    for e = 1:G_npath(v)-1;
        if bad_edge;
            continue
        end
        j = G_edge{v}(1,e);
        i = G_edge{v}(2,e);
        if G_duplicate{v}(e) && (G_order{v}(e)==1);
            if ~is_valid_sub(i,j);
                bad_edge = true;
                continue
            end
            R4curr = match_R4_sub1{i,j};
            G_nmatch{v}(e) = match_nmatch_sub1(i,j);
        elseif G_duplicate{v}(e) && (G_order{v}(e)==2);
             if ~is_valid_sub(i,j);
                 bad_edge = true;
                continue
            end
            R4curr = match_R4_sub2{i,j};
            G_nmatch{v}(e) = match_nmatch_sub2(i,j);
        else
            if ~is_valid(i,j);
                bad_edge = true; 
                continue
            end
            R4curr = match_R4{i,j};
            G_nmatch{v}(e) = match_nmatch(i,j);
        end
        R4temp = R4curr*R4temp;
    end
    if ~bad_edge;
        G_R4{v} = R4temp;
        G_valid(v) = true;
    end
end

G_R4 = G_R4(G_valid);
G_dist = G_dist(G_valid);
G_path = G_path(G_valid);
G_npath = G_npath(G_valid);
G_nmatch = G_nmatch(G_valid);
G_R3 = cellfun(@(x) x(1:3,1:3), G_R4,'uniformoutput', false);
G_t = cellfun(@(x) x(1:3,4), G_R4,'uniformoutput', false);
G_rx = (180/pi)*cellfun(@decompose_rotation_rx, G_R3);
G_ry = (180/pi)*cellfun(@decompose_rotation_ry, G_R3);
G_rz = (180/pi)*cellfun(@decompose_rotation_rz, G_R3);
G_tx = cellfun(@(x) x(1), G_t);
G_ty = cellfun(@(x) x(2), G_t);
G_tz = cellfun(@(x) x(3), G_t);
G_start = cellfun(@(x) x(1), G_path);
n_g = numel(G_path);


% Identify duplicate permutations
%{
G_is_duplicate = cell(n_S,1);
for s = 2:t_max_cxn_to1;
    G_is_duplicate{s} = mod(1:n_cp(s),2);
end

% Paths and camera locations
G_cp = cell(n_S,1);
G_t = cell(n_S,1);
for s = 2:t_max_cxn_to1;
    G_cp{s} = nan(n_cp(s),s+1);
    G_t{s} = nan(n_cp(s),s+1,3);
    ix = 1;
    for c = 1:n_comb(s);
        for p = 1:n_perm(s);
            path = G_comb{s}(c,G_perm{s}(p,:));
            G_cp{s}(ix,:) = path;
            G_t{s}(ix,:,:) = match_teff_arr(G_comb{s}(c,G_perm{s}(p,:)),:);
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
%}
% Remove paths which are too far (based on distance from all)
%{
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
%}
% Remove paths which are not connected *note combined adjacency matrix
%{
G_valid = cell(n_S,1);
% Make non-directed
match_empty_sub1 = cellfun(@isempty,match_i_sub1);
match_empty_sub2 = cellfun(@isempty,match_i_sub1);
match_empty_all = triu(cellfun(@isempty, match_i),1);
match_empty_all = match_empty_all | match_empty_all';
for g = 2:t_max_cxn_to1;
    is_reject = false(n_path(g),1);
    for p = 1:n_path(g);
        path = G_path{g}(p,:);
        if g ==3;
            if all(path == [3 5 2 3]);
            foo = 1;
            end
        end
        % First step is sub1 
        if match_empty_sub1(path(2),path(1));
                is_reject(p) = true;
                continue
        end
        for k = 2:numel(path)-2; % +1 -1
            if match_empty_all(path(k+1),path(k));
                is_reject(p) = true;
                continue
            end
        end
        if match_empty_sub2(path(end),path(end-1));
            is_reject(p) = true;
            continue
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
G_npath = cellfun(@(x) size(x,1),G_path);
%}

% Output 
% G_R4                  cell array of [R|t] for circular paths
% G_R3                  cell array of [R] for circular paths 
% G_t                   cell array of [t] for circular paths 
% G_rx                  array of x angles [degrees]
% G_ry                  array of y angles [degrees]
% G_rz                  array of z angles [degrees]
% G_tx                  array of x translations [m]
% G_ty                  array of y translations [m]
% G_tz                  array of z translations [m]
% G_*                   updated based on logical valid array
% G_nmatch              number of matches at each connection 
% n_g                   number of paths 
clear G_valid v R4temp bad_edge j i R4curr 
%% Validation: transformation parameters 
n_p = numel(unique(G_npath));
Gmaxpath = max(G_npath);
% Rotation vs. # nodes 
%{
Gminpath = min(G_npath);
Gmaxpath = max(G_npath);
tick = Gminpath:Gmaxpath; 
n_p = numel(tick);
figure; 
hold on
plot(G_npath,G_rx,'ok','markerfacecolor','r');
plot(G_npath,G_ry,'ok','markerfacecolor','g');
plot(G_npath,G_rz,'ok','markerfacecolor','b');
grid on
set(gca, 'xtick',tick);
xlabel('Number of nodes');
ylabel('Error in Degrees');
legend_str{1} = 'rx';
legend_str{2} = 'ry';
legend_str{3} = 'rz';
legend(legend_str);
%title('Error in Rotation');
filename = sprintf('Error-rxyz_vs_nodes_s%02.0f',info_site);
filepath_save = sprintf('%s%s.eps',path_save,filename);
%saveas(gcf,filepath_save,'psc2')
%}

% Make tikz figure 
%{
data_1 = cell(n_p,1);
data_2 = cell(n_p,1);
data_3 = cell(n_p,1);
labels = cell(n_p,1);
ctr = 1;
for p = 1:Gmaxpath;
    is_valid = (G_npath==p);
    if sum(is_valid);
        data_1{ctr} = G_rx(is_valid);
        data_2{ctr} = G_ry(is_valid);
        data_3{ctr} = G_rz(is_valid);
        labels{ctr} = sprintf('%g',p);
       ctr = ctr + 1;
    end
end

tikzxlabel = 'Number of nodes';
tikzylabel = 'Error in deg';
legendlabels = {'\gls{rx}', '\gls{ry}', '\gls{rz}'};
ymin = min([G_rx; G_ry; G_rz]);
ymax = max([G_rx; G_ry; G_rz]);
ylimval = [ymin ymax];

filename = sprintf('Error-rxyz_vs_nodes_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%% Translation vs. # nodes 
%{
figure; 
hold on
for g = 1:n_S
    plot(G_npath,G_tx,'ok','markerfacecolor','r');
    plot(G_npath,G_ty,'ok','markerfacecolor','g');
    plot(G_npath,G_tz,'ok','markerfacecolor','b');
end
grid on
set(gca, 'xtick',2:max(G_npath));
xlabel('Number of nodes');
ylabel('Error in meters');
legend_str{1} = 'tx';
legend_str{2} = 'ty';
legend_str{3} = 'tz';
legend(legend_str);
%title('Error in Translation');
filename = sprintf('Error-txyz_vs_nodes_s%02.0f',info_site);
filepath_save = sprintf('%s%s.eps',path_save,filename);
%saveas(gcf,filepath_save,'psc2')
%} 

% Make tikz figure 
%{
data_1 = cell(n_p,1);
data_2 = cell(n_p,1);
data_3 = cell(n_p,1);
labels = cell(n_p,1);
ctr = 1;
for p = 1:Gmaxpath;
    is_valid = (G_npath==p);
    if sum(is_valid);
        data_1{ctr} = G_tx(is_valid);
        data_2{ctr} = G_ty(is_valid);
        data_3{ctr} = G_tz(is_valid);
        labels{ctr} = sprintf('%g',p);
        ctr = ctr + 1;
    end
end

tikzxlabel = 'Number of nodes';
tikzylabel = 'Error in $m$';
legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};

ymin = min([G_tx; G_ty; G_tz]);
ymax = max([G_tx; G_ty; G_tz]);
ylimval = [ymin ymax];

filename = sprintf('Error-txyz_vs_nodes_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%% Rotation vs path distance 
%{
figure; 
hold on
plot(G_dist,G_rx,'ok','markerfacecolor','r');
plot(G_dist,G_ry,'ok','markerfacecolor','g');
plot(G_dist,G_rz,'ok','markerfacecolor','b');
grid on
xlabel('Total path distance in meters');
ylabel('Error in degrees');
legend_str{1} = 'rx';
legend_str{2} = 'ry';
legend_str{3} = 'rz';
legend(legend_str);
%title('Error in rotation');
filepath_save = sprintf('%sError-rxyz_vs_dist_s%02.0f.eps',path_save, info_site);
%saveas(gcf,filepath_save,'psc2')
%}

% Make tikz figure 
%{
dist_axes = 10:10:90;
n_d = numel(dist_axes) - 1;
data_1 = cell(n_d,1);
data_2 = cell(n_d,1);
data_3 = cell(n_d,1);
labels = cell(n_d,1);
ctr = 1;
for d = 1:n_d;
    is_valid = (G_dist>=dist_axes(d)) & (G_dist<dist_axes(d+1));
    if sum(is_valid);
        data_1{ctr} = G_rx(is_valid);
        data_2{ctr} = G_ry(is_valid);
        data_3{ctr} = G_rz(is_valid);
        labels{ctr} = sprintf('%g-%g',dist_axes(d), dist_axes(d+1));
        ctr = ctr + 1;
    end
end

tikzxlabel = 'Distance in meters';
tikzylabel = 'Error in degrees';
legendlabels = {'\gls{rx}', '\gls{ry}', '\gls{rz}'};

ymin = min([G_rx; G_ry; G_rz]);
ymax = max([G_rx; G_ry; G_rz]);
ylimval = [ymin ymax];

filename = sprintf('Error-rxyz_vs_dist_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%% Translation vs. path distance 
%{
figure; 
hold on
plot(G_dist,G_tx,'ok','markerfacecolor','r');
plot(G_dist,G_ty,'ok','markerfacecolor','g');
plot(G_dist,G_tz,'ok','markerfacecolor','b');
grid on
xlabel('Total path distance in meters');
ylabel('Error in meters');
legend_str{1} = 'tx';
legend_str{2} = 'ty';
legend_str{3} = 'tz';
legend(legend_str);
%title('Error in translation');
filepath_save = sprintf('%sError_txyz_vs_dist_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}
% Make tikz figure 
%{
dist_axes = 10:10:90;
n_d = numel(dist_axes) - 1;
data_1 = cell(n_d,1);
data_2 = cell(n_d,1);
data_3 = cell(n_d,1);
labels = cell(n_d,1);
ctr = 1;
for d = 1:n_d;
    is_valid = (G_dist>=dist_axes(d)) & (G_dist<dist_axes(d+1));
    if sum(is_valid);
        data_1{ctr} = G_tx(is_valid);
        data_2{ctr} = G_ty(is_valid);
        data_3{ctr} = G_tz(is_valid);
        labels{ctr} = sprintf('%g-%g',dist_axes(d), dist_axes(d+1));
        ctr = ctr + 1;
    end
end

tikzxlabel = 'Distance in meters';
tikzylabel = 'Error in $m$';
legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};

ymin = min([G_tx; G_ty; G_tz]);
ymax = max([G_tx; G_ty; G_tz]);
ylimval = [ymin ymax];

filename = sprintf('Error-txyz_vs_dist_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%%
       
G_min_nmatch = cellfun(@min,G_nmatch);
G_meanr = sqrt(G_rx.^2 + G_ry.^2 + G_rz.^2);
G_meant = sqrt(G_tx.^2 + G_ty.^2 + G_tz.^2);
min_m = min(G_min_nmatch);
max_m = max(G_min_nmatch);
m_axes = min_m:1:max_m;
n_m = numel(m_axes);
% Translation vs. min tie points  
%{
figure; 
clear legend_str
hold on
n_maxpath = max(G_npath);
color_npath = jet(n_maxpath);
ctr = 1;
for i = 1:n_maxpath;
    is_valid = (G_npath == i);
    if sum(is_valid);
        plot(G_min_nmatch(is_valid),G_meanr(is_valid),'ok','markerfacecolor',color_npath(i,:));
        legend_str{ctr} = sprintf('%g nodes',i);
        ctr = ctr + 1;
    end
end
grid on
xlabel('Minimum number of tie points');
ylabel('Root Square Error in Degrees');
legend(legend_str);
%title('Error in rotation');
filepath_save = sprintf('%sError-r_vs_nties_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}

% Make tikz figure 
%{
data_1 = cell(n_m,1);
data_2 = cell(n_m,1);
data_3 = cell(n_m,1);
labels = cell(n_m,1);
for m = 1:n_m;
    labels{m} = sprintf('%g',m_axes(m));
    is_valid = (G_min_nmatch==m_axes(m));
    if sum(is_valid);
        data_1{m} = G_rx(is_valid);
        data_2{m} = G_ry(is_valid);
        data_3{m} = G_rz(is_valid);
    end
end

tikzxlabel = 'Minimum number of tie points';
tikzylabel = 'Error in deg';
legendlabels = {'\gls{rx}', '\gls{ry}', '\gls{rz}'};
ymin = min([G_rx; G_ry; G_rz]);
ymax = max([G_rx; G_ry; G_rz]);
ylimval = [ymin ymax];

filename = sprintf('Error-rxyz_vs_nties_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%%

% Translation xyz vs. min tie points 
%{
figure; 
clear legend_str
hold on
n_maxpath = max(G_npath);
color_npath = jet(n_maxpath);
ctr = 1;
for i = 1:n_maxpath;
    is_valid = (G_npath == i);
    if sum(is_valid);
        plot(G_min_nmatch(is_valid),G_meant(is_valid),'ok','markerfacecolor',color_npath(i,:));
        legend_str{ctr} = sprintf('%g nodes',i);
        ctr = ctr + 1;
    end
end
grid on
xlabel('Minimum number of tie points');
ylabel('Root square error in m');
legend(legend_str);
%title('Error in translation');
filepath_save = sprintf('%sError-t_vs_nties_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}

% Make tikz figure 
%{
data_1 = cell(n_m,1);
data_2 = cell(n_m,1);
data_3 = cell(n_m,1);
labels = cell(n_m,1);
for m = 1:n_m;
    labels{m} = sprintf('%g',m_axes(m));
    is_valid = (G_min_nmatch==m_axes(m));
    if sum(is_valid);
        data_1{m} = G_tx(is_valid);
        data_2{m} = G_ty(is_valid);
        data_3{m} = G_tz(is_valid);
    end
end

tikzxlabel = 'Minimum number of tie points';
tikzylabel = 'Error in meters';
legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};
ymin = min([G_tx; G_ty; G_tz]);
ymax = max([G_tx; G_ty; G_tz]);
ylimval = [ymin ymax];

filename = sprintf('Error-txyz_vs_nties_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%%

% Rotation xyz vs. mean connections
%{
figure; 
clear legend_str
hold on
plot(G_min_nmatch(is_valid),G_rx(is_valid),'ok','markerfacecolor','r');
plot(G_min_nmatch(is_valid),G_ry(is_valid),'ok','markerfacecolor','g');
plot(G_min_nmatch(is_valid),G_rz(is_valid),'ok','markerfacecolor','b');
grid on
xlabel('Minimum number of tie points');
ylabel('Error in degrees');
legend_str{1} = 'rx';
legend_str{2} = 'ry';
legend_str{3} = 'rz';
legend(legend_str);
%title('Error in rotation');
filepath_save = sprintf('%sError-rxy_vs_ntie_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}

% Translation xyz vs. mean connections
%{
figure; 
clear legend_str
hold on
plot(G_min_nmatch(is_valid),G_tx(is_valid),'ok','markerfacecolor','r');
plot(G_min_nmatch(is_valid),G_ty(is_valid),'ok','markerfacecolor','g');
plot(G_min_nmatch(is_valid),G_tz(is_valid),'ok','markerfacecolor','b');
grid on
xlabel('Minimum number of tie points');
ylabel('Error in m');
legend_str{1} = 'tx';
legend_str{2} = 'ty';
legend_str{3} = 'tz';
legend(legend_str);
%title('Error in translation');
filepath_save = sprintf('%sError-rxy_vs_ntie_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}
% Outputs 
% G_min_nmatch          minimum number of tie points 
clear filepath_save legend_str G_meanr G_meant 
clear n_maxpath color_npath ctr i is_valid 
%% Compare tree locations 

G_tree_hat = cell(n_g,1);
G_tree_true = cell(n_g,1);
GR4_eye = cellfun(@double,repmat({[eye(3) zeros(3,1); zeros(1,3) 1]},n_g,1),'Un',0);
for s = 1:n_S;
    is_valid = (G_start==s);
    P4_temp = [P_LCS{s}; ones(1,size(P_LCS{s},2))];
    GR4_temp = G_R4(is_valid);
    Geye_temp = GR4_eye(is_valid);
    G_tree_hat(is_valid) = cellfun(@(x) x*P4_temp, GR4_temp, 'uniformoutput',false);
    G_tree_true(is_valid) =  cellfun(@(x) x*P4_temp, Geye_temp, 'uniformoutput',false);
end

n_inG = cellfun(@(x) size(x,2), G_tree_hat);
n_all = sum(n_inG);

all_tree_hat = zeros(n_all,3);
all_tree_true = zeros(n_all,3); 
all_npath = zeros(n_all,1);
all_dist = zeros(n_all,1);
all_minmatch = zeros(n_all,1);
start_ix = [0; cumsum(n_inG)]+1;
end_ix = circshift(start_ix, [-1,0])-1;
start_ix = start_ix(1:end-1);
end_ix = end_ix(1:end-1);

for g = 1:n_g;
    all_tree_hat(start_ix(g):end_ix(g),:) = G_tree_hat{g}(1:3,:)';
    all_tree_true(start_ix(g):end_ix(g),:) = G_tree_true{g}(1:3,:)';
    all_npath(start_ix(g):end_ix(g)) = repmat(G_npath(g),[n_inG(g),1]);
    all_dist(start_ix(g):end_ix(g)) = repmat(G_dist(g),[n_inG(g),1]);
    all_minmatch(start_ix(g):end_ix(g)) = repmat(G_min_nmatch(g),[n_inG(g),1]);
end

all_tree_exy = sqrt(sum((all_tree_hat(:,1:2) - all_tree_true(:,1:2)).^2,2)); 
all_tree_ex = all_tree_hat(:,1) - all_tree_true(:,1); 
all_tree_ey = all_tree_hat(:,2) - all_tree_true(:,2); 
all_tree_ez = all_tree_hat(:,3) - all_tree_true(:,3); 

%% Tree location vs. number of nodes 

%{
clear legend_str
figure;
hold on
plot(all_npath, all_tree_exy,'+','markeredgecolor','k');
axisval = axis;
legend_str{1} = 'x';
legend_str{2} = 'y';
%legend_str{3} = 'z';
legend(legend_str);
xlabel('Numbre of nodes');
ylabel('Error in xy');
%title('Tree position accuracy after registration');
filename = sprintf('Error-tloc_vs_nodes_s%02.0f',info_site);
filepath_save = sprintf('%s%s.eps',path_save,filename);
%saveas(gcf,filepath_save,'psc2')
%}

% Make tikz figure 
%{
n_p = max(all_npath) - min(all_npath);
data1 = cell(n_p,1);
data2 = cell(n_p,1);
data3 = cell(n_p,1);
labels = cell(n_p,1);
data_isvalid = false(n_p,1);
ctr = 1;
for p = 1:max(all_npath);
    is_valid = (all_npath==p);
    if sum(is_valid);
        data_1{p} = all_tree_ex(is_valid);
        data_2{p} = all_tree_ey(is_valid);
        data_3{p} = all_tree_ez(is_valid);
        data_isvalid(p) = true;
        labels{ctr} = sprintf('%g',p);
        ctr = ctr + 1;
    end
end
data_1 = data_1(data_isvalid);
data_2 = data_2(data_isvalid);
data_3 = data_3(data_isvalid);

tikzxlabel = 'Number of nodes';
tikzylabel = 'Error';
legendlabels = {'$\gls{tree:xg}(x) - \hat{\gls{tree:xg}}(x)$',...
    '$\gls{tree:xg}(y) - \hat{\gls{tree:xg}}(y)$',...
    '$\gls{tree:xg}(z) - \hat{\gls{tree:xg}}(z)$'};

ymin = min([all_tree_ex; all_tree_ey; all_tree_ez]);
ymax = max([all_tree_ex; all_tree_ey; all_tree_ez]);
ylimval = [ymin ymax];

filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%% Tree location vs path distance 

%{
clear legend_str
figure;
hold on
plot(all_dist,all_tree_ex,'ok','markerfacecolor','r');
plot(all_dist,all_tree_ey,'ok','markerfacecolor','g');
%plot(all_dist,all_tree_ez,'ok','markerfacecolor','b');
axisval = axis;
legend_str{1} = 'x';
legend_str{2} = 'y';
%legend_str{3} = 'z';
legend(legend_str);
xlabel('Total path distance');
ylabel('Error');
%title('Tree position accuracy after registration');
filename = sprintf('Error-tloc_vs_dist_s%02.0f',info_site);
filepath_save = sprintf('%s%s.eps',path_save,filename);
saveas(gcf,filepath_save,'psc2');
%}

% Make tikz figure 
%{
dist_axes = 10:10:90;
n_d = numel(dist_axes) - 1;
data_1 = cell(n_d,1);
data_2 = cell(n_d,1);
data_3 = cell(n_d,1);
labels = cell(n_d,1);
ctr = 1;
for d = 1:n_d;
    is_valid = (all_dist>=dist_axes(d)) & (all_dist<dist_axes(d+1));
    if sum(is_valid);
        data_1{ctr} = all_tree_ex(is_valid);
        data_2{ctr} = all_tree_ey(is_valid);
        data_3{ctr} = all_tree_ez(is_valid);
        labels{ctr} = sprintf('%g-%g',dist_axes(d), dist_axes(d+1));
        ctr = ctr + 1;
    end
end

tikzxlabel = 'Distance in meters';
tikzylabel = 'Error in $m$';
legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};

ymin = min([all_tree_ex; all_tree_ey; all_tree_ez]);
ymax = max([all_tree_ex; all_tree_ey; all_tree_ez]);
ylimval = [ymin ymax];

filename = sprintf('Error-tloc_vs_dist_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%%
%
clear legend_str
figure;
hold on
for g = 1:n_g;
plot(G_tree_true{g}(1,:),G_tree_hat{g}(1,:),'+','markeredgecolor',[.5 .5 .5]);
plot(G_tree_true{g}(2,:),G_tree_hat{g}(2,:),'x','markeredgecolor',[.8 .8 .8]);
%plot(G_tree_true{g}(3,:),G_tree_hat{g}(3,:),'+','markeredgecolor','b');
end
axisval = axis;
plot(axisval(1:2),axisval(3:4),'-k');
legend_str{1} = 'x';
legend_str{2} = 'y';
%legend_str{3} = 'z';
legend(legend_str);
xlabel('Truth location');
ylabel('Estimated location');
%title('Tree position accuracy after registration');
filepath_save = sprintf('%sError-treeloc_s%02.0f.eps',path_save, info_site);
%saveas(gcf,filepath_save,'psc2')
%}

    for p = 1:G_npath(g);
        path = G_path{g}(p,:);
        i = path(1);
        j = path(2);
        loop_P_true{g}{p} = P_LCS{i}(:,match_i_sub2{i,j});
        loop_P_hat{g}{p} = loop_R{g}{p}*P_LCS{i}(:,match_i_sub2{i,j}) +...
            repmat(loop_t{g}{p},[1,numel(match_i_sub2{i,j})]);
    end


all_true = [];
all_hat = [];
all_g1 = [];
all_npath = [];
for g = 1:n_S;
    for p = 1:G_npath(g);
       all_true = [all_true loop_P_true{g}{p}];
       all_hat = [all_hat loop_P_hat{g}{p}];
    end
end
all_error = sqrt((all_true - all_hat).^2); 

fprintf('\nRMSE x = %3.3f m\n', mean(all_error(1,:)))
fprintf('RMSE y = %3.3f m\n', mean(all_error(2,:)))
fprintf('RMSE z = %3.3f m\n', mean(all_error(3,:)))

%
figure;
hold on
plot(all_true(1,:),all_hat(1,:),'ok','markerfacecolor','r');
plot(all_true(2,:),all_hat(2,:),'ok','markerfacecolor','g');
plot(all_true(3,:),all_hat(3,:),'ok','markerfacecolor','b');
axisval = axis;
plot(axisval(1:2),axisval(3:4),'-k');
legend_str{1} = 'x';
legend_str{2} = 'y';
legend_str{3} = 'z';
legend(legend_str);
xlabel('Truth location');
ylabel('Estimated location');
title('Tree position accuracy after registration');
filepath_save = sprintf('%sError-treeloc_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}
%% Validation: Point cloud error 
% Initial transformation of PLY
for g = 1:n_S;
    for p = 1:G_npath(g);
        fprintf('\nWriting ply error\n');
        path = G_path{g}(p,:);
        i = path(1);
        tmp_plot = P_plot(i);
        [vertex, ~] = read_ply(filepath_ply{i});
        data_xyz = vertex(:,1:3)';
       % color_ply = vertex(:,4:6);
        data_xyzhat = loop_R{g}{p}*data_xyz+...
            repmat(loop_t{g}{p},[1,size(vertex,1)]);
        data_e = abs(mean((data_xyz - data_xyzhat),1));
        
        figure;
        xbin = linspace(0,0.05,100);
        count = hist(data_e, xbin);
        plot(xbin,count,'-x');
        color_e = vec2cmap2(data_e,'jet', 0,.05);
        filepath_ply_e = sprintf('%serror_%03.0f-%02.0f.ply', ...
        path_ply{j}, info_site, P_plot(i));
        write2ply(filepath_ply_e,data_xyz', color_e);
    end
end

for j = 1:n_S;
    [vertex, ~] = read_ply(filepath_ply{j});

    filepath_ply_reg{j} = sprintf('%spoints_full_%03.0f-%02.0f-%02.0f.ply', ...
        path_ply{j}, info_site, P_plot(j),P_plot(1));
    xyz2t = (match_Reff{1,j}*[data_x2 data_y2 data_z2]') + ...
        repmat(match_teff{1,j},1,numel(data_x2));
    write2ply(filepath_ply_reg{j},xyz2t', color);
end
clear match_teff_arr vertex data_x2 data_y2 data_z2 color_ply xyz2t
%}

%% Validation: Error field 

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



%% HERE BEGINS BUNDLE ADJUSTMENT OPTIMIZATION 
%% Find unique points
fprintf('\nFind unique points\n');

% Gather all points in WCS (i = 1) from each camera
all_Pi_x = [];
all_Pi_y = [];
all_Pi_z = [];
all_Pi_r = [];
for j = 1:n_S;
    if ~isempty(match_Pi_all{1,j});
        all_Pi_x = [all_Pi_x match_Pi_all{1,j}(1,:)]; %#ok<*AGROW>
        all_Pi_y = [all_Pi_y match_Pi_all{1,j}(2,:)];
        all_Pi_z = [all_Pi_z match_Pi_all{1,j}(3,:)];
        all_Pi_r = [all_Pi_r match_rad_all{1,j}']; %match_R_all or P_rad is the same
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
t_xyz = .5.^2;
t_r = 0.2.^2;
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
    scatter3(all_Pi_x, all_Pi_y, all_Pi_z,30,all_unique,'filled','markeredgecolor','k');
    scatter3(unique_x,unique_y,unique_z,120,'k');%,'markersize',20,...
    %  'markerfacecolor',[1 1 1],'alpha',1, 'linewidth',.5);
    legend_str{1} = 'All points';
    legend_str{2} = 'Unique points';
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    legend(legend_str,'location','best');
    axis equal
    axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
    grid on
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

data_LCSi = nan(n_unique,n_S,3);
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

match_teff_arr = zeros(n_S,3);
for s = 1:n_S;
    match_teff_arr(s,:) = match_teff{s}';
end

G_cp = cell(n_S,1);
G_t = cell(n_S,1);
for s = 2:t_max_cxn_to1;
    G_cp{s} = nan(n_cp(s),s);
    G_t{s} = nan(n_cp(s),s,3);
    ix = 1;
    for c = 1:n_comb(s);
        for p = 1:n_perm(s);
            path = G_comb{s}(c,G_perm{s}(p,:));
            G_cp{s}(ix,:) = path;
            G_t{s}(ix,:,:) = match_teff_arr(G_comb{s}(c,G_perm{s}(p,:)),:);
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
for s = 2:t_max_cxn_to1;
    clear legend_str
    color_cxn = jet(n_path(s));
    figure;
    xlabel('x');ylabel('y');zlabel('z');
    view(0,90);
    hold on
    for p = 1:n_path(s);
         plot3(G_t{s}(p,:,1),G_t{s}(p,:,2),G_t{s}(p,:,3),'-.','color',color_cxn(p,:), 'linewidth',2);
    end
    %legend_str{c+2} = sprintf('Paths of %g',s);
    title('Points and camera positions');
    %legend(legend_str,'location','best');
    grid on
end
%}


%% Levenberg Marquardt
fprintf('\nLevenberg Marquardt\n');

% Selects paths used in LM
t_max_cxn = 3;
t_max_cxndist = 18;
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
axis(1.5*[i_xmin i_xmax i_ymin i_ymax i_zmin i_zmax]);
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


