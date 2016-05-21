%% Registration Validation 
% This script loads data from kelbe_registration and makes new graphs which
% are independent for validation of registration parameters 

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
options_loadmatch = true;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_temp = 'D:\Users\djk2312\Documents\MATLAB\temp\';
load('validation.mat');
% Initialize program-level options and paths
% Outputs
% info_site
% P_LCS                 points in local coordinate system (cell)
% P_rad                 radius of points
% P_n                   number of points
% is_valid              adjacency matrix for scan connections
% match_nmatch          number of matches for a given link 
% match_R               (i,j) is pairwise rotation from j into i
% match_t               (i,j) is pairwise translation from j into i
% match_Reff            effective rotation from j into WCS
% match_teff            effective translation from j into WCS
% ************************
% options               options
% paths                 paths
%% Repeat registration with disjoint subsets

filepath_match_R_sub1 = sprintf('%s%s',path_temp,'match_R_sub1.mat');
filepath_match_t_sub1 = sprintf('%s%s',path_temp,'match_t_sub1.mat');
filepath_match_i_sub1 = sprintf('%s%s',path_temp, 'match_i_sub1.mat');
filepath_match_j_sub1 = sprintf('%s%s',path_temp, 'match_j_sub1.mat');
filepath_match_i_ix1 = sprintf('%s%s',path_temp, 'match_i_ix1.mat');
filepath_match_R_sub2 = sprintf('%s%s',path_temp, 'match_R_sub2.mat');
filepath_match_t_sub2 = sprintf('%s%s',path_temp, 'match_t_sub2.mat');
filepath_match_i_sub2 = sprintf('%s%s',path_temp, 'match_i_sub2.mat');
filepath_match_j_sub2 = sprintf('%s%s',path_temp, 'match_j_sub2.mat');
filepath_match_i_ix2 = sprintf('%s%s',path_temp, 'match_i_ix2.mat');
filepath_match_all_sub1 = sprintf('%s%s',path_temp, 'match_all_sub1.mat');
filepath_match_all_sub2 = sprintf('%s%s',path_temp, 'match_all_sub2.mat');

if options_loadmatch == false;
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

% Outputs
% filepath_match_*_sub  filepath of matlab data related to subsets 
% match_*_sub           pairwise matching rotation, translation, indices
clear a b c d e f g h i j ix_randi ix_sub1 ix_sub2 mid 
%% Remove bad matches 
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
% is_valid_sub          adjacency matrix for subset scan connections (directed)
% match_nmatch_sub      number of matches for a link (subset)
clear is_valid_sub1 is_valid_sub2
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
% Only for error checking/ visualization 
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
clear all_perm i j 
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
clear is is_i is_j i j v V_temp n_vtemp  
%% Create cell array and determine additional parameters 
G_npath = sum(~isnan(G_val),2);
G_path = cell(n_v,1);
G_edge = cell(n_v,1);
G_duplicate = cell(n_v,1);
G_order = cell(n_v,1);

for v = 1:n_v;
    path = G_val(v,~isnan(G_val(v,:)));
    edges = [path(1:end-1); circshift(path(1:end-1),[1,-1])];
    [edges_sort, ~] = sort(edges,1);
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
% Use duplicate information to traverse different "legs" of same path
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
clear G_valid v R4temp bad_edge j i R4curr e ans 
foo = 1;

%% These were commented... not sure if relevant? 
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

%% Product: point cloud error 
kelbe_reg_product_pcerror( n_S,G_npath,G_path,P_plot,...
    filepath_ply,loop_R,loop_t)
%% Product: Error field 
kelbe_product_errfield(  )

%% Find matching tree locations 

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
    %all_minmatch(start_ix(g):end_ix(g)) = repmat(G_min_nmatch(g),[n_inG(g),1]);
end

all_tree_exy = sqrt(sum((all_tree_hat(:,1:2) - all_tree_true(:,1:2)).^2,2)); 
all_tree_ex = all_tree_hat(:,1) - all_tree_true(:,1); 
all_tree_ey = all_tree_hat(:,2) - all_tree_true(:,2); 
all_tree_ez = all_tree_hat(:,3) - all_tree_true(:,3); 


save('tempvalidation.mat', 'G_npath', 'G_rx', 'G_ry', 'G_rz',...
    'G_tx', 'G_ty', 'G_tz', 'info_site', 'all_dist', 'all_tree_exy', ...
    'G_tree_true', 'G_tree_hat', 'n_g');
