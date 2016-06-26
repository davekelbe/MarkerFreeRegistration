function [ tree_smdi ] = seg2tree_Jan( seg_z, seg_y, seg_r, seg_anorm, seg_iter,...
    is_valid,t_rsearch, path_mat )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

%%
n_seg = size(seg_r,1);
% Find neighbor for each segment
seg_neighbor = zeros(n_seg,1);
for s = 1:n_seg;
    if ~is_valid(s);
        continue;
    end
    a = seg_anorm(s,:);
    r = mean(seg_r(s,:));
    c = (seg_z(s,:) + seg_y(s,:))/2;
    seg_d = seg_z - repmat(c,[n_seg,1]);
    [U,~,~] = svd(a');
    b1 = U(:,1);
    b2 = U(:,2);
    b3 = U(:,3);
    %pproj = @(u) u*u'/(u'*u);
    %qproj = @(u) eye(length(u))-u*u'/(u'*u);
    seg_d1 = dot(seg_d',repmat(b1,[1,n_seg]))';
    seg_d2 = dot(seg_d',repmat(b2,[1,n_seg]))';
    seg_d3 = dot(seg_d',repmat(b3,[1,n_seg]))';
    seg_d23 = sqrt(seg_d2.^2 + seg_d3.^2);
    seg_d23(~is_valid) = inf;
    temp_eligible = true(n_seg,1);
    temp_eligible(s) = false;
    is_eligible = (sign(seg_d1)==1)&(seg_d23<max(2*r,t_rsearch))&temp_eligible;
    temp_iter = seg_iter(is_eligible);
    temp_d1 = seg_d1(is_eligible);
    if numel(temp_d1);
        [~,imin] = min(temp_d1);
        seg_neighbor(s) = temp_iter(imin);
    end
end
% Find tree ownership for all segments

[~,ix] = sort(seg_z(:,3));
%seg_neighbor_old = seg_neighbor;
seg_neighbor_new = seg_neighbor(ix);
%seg_iter_old = seg_iter;
seg_iter_new = seg_iter(ix);
seg_iter = seg_iter_new;
seg_neighbor = seg_neighbor_new;
seg_tree = zeros(n_seg,1);
t = 1;
is_used = false(n_seg,1);
%is_used(~is_valid) = true;
for s = 1:n_seg;
    % s = seg_iter(dummy);
    restart = true;
    if is_used(s);
        continue
    end
    %        s = s;
    seg_tree(s) = t;
    is_used(s) = true;
    %if (seg_neighbor(s) == 0)
        % Jan: Commented to allow for single-segment trees
        % t = t + 1;
        % restart = true;
        % continue
    %else
        restart = false;
        while ~restart
            %s = seg_neighbor(s);
            if seg_neighbor(s) ~= 0;
                s = find(seg_iter == seg_neighbor(s));
            end
            if numel(s)>1;
                s = s(1);
            end
            if ~is_used(s) && (seg_neighbor(s) ~= 0)
                seg_tree(s) = t;
                is_used(s) = true;
            elseif ~is_used(s) && (seg_neighbor(s) == 0)
                seg_tree(s) = t;
                is_used(s) = true;
                restart = true;
                t = t + 1;
                continue
            else
                restart = true;
                t = t + 1;
            end
        end
    %end
end

[~,ix2] = sort(ix);
seg_iter_new_sort = seg_iter_new(ix2);
seg_tree = seg_tree(ix2);
[n, bin] = histc(seg_tree, unique(seg_tree));
multiple = find(n > 1);
index    = find(ismember(bin, multiple));

clear cyl
cyl(1,:) = seg_y(:,1);
cyl(2,:) = seg_y(:,2);
cyl(3,:) = seg_y(:,3);
cyl(4,:) = seg_r(:,1);
cyl(5,:) = seg_z(:,1);
cyl(6,:) = seg_z(:,2);
cyl(7,:) = seg_z(:,3);
cyl(8,:) = seg_r(:,2);
n_dup = numel(unique(seg_tree(index)));
color = zeros(n_seg,3);
cmap = jet(n_dup);
ix = randperm(n_dup)';
cmap = cmap(ix,:);
for d = 1:n_dup;
    is_dup = (seg_tree == multiple(d));
    color(is_dup,1) = cmap(d,1);
    color(is_dup,2) = cmap(d,2);
    color(is_dup,3) = cmap(d,3);
end
color255 = uint8(255*(color'));
%cyl = cyl*1000;

%{
figure;
hold on;
options_cyl_alpha = .7;
for s = 1:n_seg;
[Coneh,End1h,End2h] = Cone(seg_z(s,:),...
            seg_y(s,:),...
            [seg_r(s,1) seg_r(s,2)],...
            30,...
            color(s,:),1,0);
        set(Coneh,'facealpha',options_cyl_alpha)
        set(End1h,'facealpha',options_cyl_alpha)
        set(End2h,'facealpha',options_cyl_alpha)
end
%}
%{
filepath_segments_valid = sprintf('%sseg_valid_%03.0f-%02.0f.ply',path_ply,info_site,info_plot);
if ~exist(filepath_segments_valid,'file');
cyl2ply(filepath_segments_valid,cyl,color255);
end
path_local = 'Z:\Desktop\Local\';
filepath_segments = sprintf('%sseg_color_%03.0f-%02.0f.ply',path_local,info_site,info_plot);
if ~exist(filepath_segments,'file');
cyl2ply(filepath_segments,cyl,color255);
end
%}


% Make tree structure

tree_unique= unique(seg_tree);
n_tree = numel(tree_unique);
%tree = zeros(n_tree,1);
t = 1;
tree.loc = [];
for i = 1:n_tree;
    is_tree = (seg_tree==tree_unique(i));
    if sum(is_tree)<2;
     %   continue
    end
    tree_z = seg_z(is_tree,:);
    tree_y = seg_y(is_tree,:);
    tree_r = seg_r(is_tree,:);
    tree_l = [tree_z; tree_y];
    tree_r = [tree_r(:,1); tree_r(:,2)];
    [~,ixsort] = sort(tree_l(:,3));
    tree(t).loc = tree_l(ixsort,:)';
    tree(t).r = tree_r(ixsort)';
    tree(t).color =  [127 127 127]';
    t = t + 1;
end
%
if t == 1;
    tree_smdi.loc = [];
    tree_smdi.r = [];
    tree_smdi.color = [];
    return
end
n_tree = numel(tree);

%}

tree_min = zeros(n_tree,3);
tree_max = zeros(n_tree,3);

% Combine trees
for t = 1:n_tree;
    tree_min(t,:) = tree(t).loc(:,1);
    tree_max(t,:) = tree(t).loc(:,end);
end
tree_a = tree_max - tree_min;
%[~,ix] = sort(tree_max(:,3));

t_rsearch = .7;
tree_iter = 1:n_tree;
tree_match = zeros(n_tree,1);
for t = 1:n_tree;
    [U,~,~] = svd(tree_a(t,:)');
    b1 = U(:,1);
    b2 = U(:,2);
    b3 = U(:,3);
    tree_d = tree_min - repmat(tree_max(t,:),[n_tree,1]);
    tree_d1 = dot(tree_d',repmat(b1,[1,n_tree]))';
    tree_d2 = dot(tree_d',repmat(b2,[1,n_tree]))';
    tree_d3 = dot(tree_d',repmat(b3,[1,n_tree]))';
    tree_d23 = sqrt(tree_d2.^2 + tree_d3.^2);
    %tree_d23(~is_valid) = inf;
    %temp_eligible = true(n_tree,1);
    %temp_eligible(s) = false;
    is_eligible = (sign(tree_d1)==1)&(tree_d23<max(4*r,t_rsearch));%&temp_eligible;
    temp_iter = tree_iter(is_eligible);
    temp_d1 = tree_d1(is_eligible);
    if numel(temp_d1);
        if numel(temp_d1) > 1;
            foo = 1;
        end
        [~,imin] = min(temp_d1);
        tree_match(t) = temp_iter(imin);
    end
end

% Match trees #2
is_used = false(n_tree);
i = 1;
for t = 1:n_tree;
    if is_used(t);
        continue;
    end
    if tree_match(t)>0;
        t2 = tree_match(t);
        tree2(i).loc = [tree(t).loc tree(t2).loc];
        tree2(i).r = [tree(t).r tree(t2).r];
        %        tree2(i).color = [tree(t).color tree(t2).color];
        is_used(t2) = true;
        if tree_match(t2)>0;
            t3 = tree_match(t2);
            tree2(i).loc = [tree2(i).loc tree(t3).loc];
            tree2(i).r = [tree2(i).r tree(t3).r];
            is_used(t3) = true;
            %            tree2(i).color = [tree2(i).color tree(t3).loc];
        end
    else
        tree2(i).loc = tree(t).loc;
        tree2(i).r = tree(t).r;
        %   tree2(i).color = tree(t).color;
    end
    i = i + 1;
end

n_tree2 = numel(tree2);
color = jet(n_tree2);
for t = 1:n_tree2;
    tree2(t).color = color(t,:);
end
%{
figure;
hold on;
options_cyl_alpha = .7;
for t = 1:n_tree2;
    for c = 1:numel(tree2(t).r)-2;
        
[Coneh,End1h,End2h] = Cone(tree2(t).loc(:,c),...
            tree2(t).loc(:,c+1),...
            [tree2(t).r(c) tree2(t).r(c+1)],...
            30,...
            tree2(t).color,1,0);
        set(Coneh,'facealpha',options_cyl_alpha)
        set(End1h,'facealpha',options_cyl_alpha)
        set(End2h,'facealpha',options_cyl_alpha)
    end
end
%}

% Interpolate and  Resample
tree = tree2;
n_tree = n_tree2;
%t_z_resample = .25;
%t_z_interp = 0.125;
t_z_resample = 1;
t_z_interp = 0.25;
cmap = jet(n_tree);
cmap = double(uint8(255*cmap));
counter = 1;
for t = 1:n_tree;
    %fprintf('\n Tree %2.0f of %2.0f\n', t, n_tree)
    minz = tree(t).loc(3,1);
    maxz = tree(t).loc(3,end);
    if t == 23;
        foo = 1;
    end
    zstep = minz:t_z_resample:maxz;%+t_z_resample;
    if numel(zstep)==1;
       % continue
    end
    linear_interp_z = (minz:t_z_interp:maxz)';
    linear_interp_x = interp1(tree(t).loc(3,:),tree(t).loc(1,:),linear_interp_z);
    linear_interp_y = interp1(tree(t).loc(3,:),tree(t).loc(2,:),linear_interp_z);
    PP = splinefit(linear_interp_z', [linear_interp_x linear_interp_y]',zstep,0.75);
    xystep = ppval(PP,zstep);
    %{
    figure;
    plot3(tree(t).loc(1,:),tree(t).loc(2,:),tree(t).loc(3,:),'.-k')
    hold on
    plot3(xystep(1,:),xystep(2,:),zstep,'.-r')
     axis equal
    %}
    p = polyfit(tree(t).loc(3,:), tree(t).r(:)',1);
    rfit = polyval(p,zstep);
    %{
    figure;
    plot(tree(t).loc(3,:),tree(t).r(:),'.-k')
    hold on;
    plot(zstep,rfit,'-r');
    %}
    tree_sm(counter).loc = [xystep; zstep];
    tree_sm(counter).r = rfit;
    tree_sm(counter).color = cmap(t,:)';
    counter = counter + 1;
end
n_tree = numel(tree_sm);
%{
path_local = 'Z:\Desktop\Local\';
filepath_tree_sm = sprintf('%stree_sm_%03.0f-%02.0f.ply',path_local,info_site,info_plot);
if ~exist(filepath_tree_sm,'file');
tree2ply(filepath_tree_sm,tree_sm,12);
end
%}


% Intersect with DEM
filepath_dem_qx = sprintf('%s%s%s',path_mat,'dem_qx','.mat');
filepath_dem_qy = sprintf('%s%s%s',path_mat,'dem_qy','.mat');
filepath_dem_qz = sprintf('%s%s%s',path_mat,'dem_qz','.mat');

load(filepath_dem_qx);
load(filepath_dem_qy);
load(filepath_dem_qz);
%dem_qx = dem_qx(1:5:end,1:5:end);
%dem_qy = dem_qy(1:5:end,1:5:end);
%dem_qz = dem_qz(1:5:end,1:5:end);

filepath_tree_dem = sprintf('%s%s%s',path_mat,'tree_dem','.mat');
if false;%exist(filepath_tree_dem,'file');
    load(filepath_tree_dem);
else
    dem_qz_nonan = inpaint_nans(dem_qz);
    tree_dem = nan(n_tree,3);
    [nx, ny] = size(dem_qx);
    xv = 1:nx-1;
    yv = 1:ny-1;
    [xxv yyv] = meshgrid(xv, yv);
    xxv = xxv(:);
    yyv = yyv(:);
    [ vertices, indices ] = ply_qxqyqz_tri( dem_qx,dem_qy,dem_qz_nonan );
    for t = 1:n_tree
        %fprintf('\nTree %d of %d\n', t, n_tree)
        line = [tree_sm(t).loc(:,2); 10*(tree_sm(t).loc(:,1) - tree_sm(t).loc(:,2)) ]';
        inters = intersectLineMesh3d(line, vertices, indices);
        if ~isempty(inters)
            tree_dem(t,:) = inters;
        end
    end
    %
    %     for t = 1:n_tree;
    %
    %         line = [tree_sm(t).loc(:,2); 10*(tree_sm(t).loc(:,1) - tree_sm(t).loc(:,2)) ]';
    %         done = false;
    %         i=1;
    %         while ~done;
    %             ix = xxv(i);
    %             iy = yyv(i);
    %             poly = [ dem_qx(ix,iy) dem_qy(ix,iy), dem_qz_nonan(ix,iy);...
    %                 dem_qx(ix+1,iy) dem_qy(ix+1,iy), dem_qz_nonan(ix+1,iy);...
    %                 dem_qx(ix+1,iy+1) dem_qy(ix+1,iy+1), dem_qz_nonan(ix+1,iy+1);...
    %                 dem_qx(ix,iy+1) dem_qy(ix,iy+1), dem_qz_nonan(ix,iy+1)];
    %             [inter ~] = intersectRayPolygon3d(line,poly);
    %             if ~isnan(inter(1))
    %                 tree_dem(t,:) = inter;
    %                 done = true;
    %             end
    %             i = i+1;
    %         end
    %     end
    %     save(filepath_tree_dem,'tree_dem');
    % end
    
    average_dem = min(tree_dem(:,3));
    for t = 1:n_tree
        if isnan(tree_dem(t,1));
            tree_dem(t,1) = tree_sm(t).loc(1);
            tree_dem(t,2) = tree_sm(t).loc(2);
            tree_dem(t,3) = average_dem;
        end
    end
    
    %Re-interpolate
    
    for t = 1:n_tree;
        if ~isnan(tree_dem(t,1));
            tree_smd(t).loc = [tree_dem(t,:)' tree_sm(t).loc];
            tree_smd(t).r = [tree_sm(t).r(1) tree_sm(t).r];
            tree_smd(t).color = tree_sm(t).color;
        end
    end
    
    % Reinterpolate
    for t = 1:n_tree;
        % fprintf('\n Tree %2.0f of %2.0f\n', t, n_tree)
        minz = tree_smd(t).loc(3,1);
        maxz = tree_smd(t).loc(3,end);
        zstep = minz:t_z_resample:maxz;%+t_z_resample;
        linear_interp_z = (minz:t_z_interp:maxz)';
        linear_interp_x = interp1(tree_smd(t).loc(3,:),tree_smd(t).loc(1,:),linear_interp_z);
        linear_interp_y = interp1(tree_smd(t).loc(3,:),tree_smd(t).loc(2,:),linear_interp_z);
        PP = splinefit(linear_interp_z', [linear_interp_x linear_interp_y]',zstep,0.75);
        xystep = ppval(PP,zstep);
        %{
    figure;
    plot3(tree(t).loc(1,:),tree(t).loc(2,:),tree(t).loc(3,:),'.-k')
    hold on
    plot3(xystep(1,:),xystep(2,:),zstep,'.-r')
     axis equal
        %}
        p = polyfit(tree_smd(t).loc(3,:), tree_smd(t).r(:)',1);
        rfit = polyval(p,zstep);
        %{
    figure;
    plot(tree(t).loc(3,:),tree(t).r(:),'.-k')
    hold on;
    plot(zstep,rfit,'-r');
        %}
        tree_smdi(t).loc = [xystep; zstep];
        tree_smdi(t).r = rfit;
        tree_smdi(t).color = cmap(t,:)';
    end
    
    
    
    %{
info_site = 1;
info_plot = 1;
path_local = 'Z:\Desktop\Local\';
filepath_tree_smd = sprintf('%stree_smd_%03.0f-%02.0f.ply',path_local,info_site,info_plot);
if ~exist(filepath_tree_smd,'file');
tree2ply(filepath_tree_smd,tree_smd,12);
end
filepath_tree_smdi = sprintf('%stree_smdi_%03.0f-%02.0f.ply',path_local,info_site,info_plot);
if ~exist(filepath_tree_smdi,'file');
tree2ply(filepath_tree_smdi,tree_smdi,12);
end
    %}
end

