%% KELBE REGISTRATION 
% Clear without breakpoints 
tmp = dbstatus;
save('tmp.mat','tmp');
clear all; close all; clc;
load('tmp.mat');
dbstop(tmp);
clear tmp;
delete('tmp.mat');
% Default settings 
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
% Input to function are stem maps derived from lidar.m
fprintf('\nLoad points\n');

% Filename Lookup for Harvard Data 
info_exp = 'Harvard';
info_suffix = '03-01';
info_slash = '\'; 
info_site = 31;
path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
D = dir(path_site); 
ctr = 1;
% Only load valid plots 
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
    plot = str2double(info_plot);
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
%% Effective (WCS) R and t using Dijkstra's 

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
%% Save variables for validation 

save('validation.mat', 'n_S','P_n','P_LCS', 'P_rad','is_valid',...
    'match_R','match_t', 'match_teff', 'match_nmatch', 'info_site');


