% Unify stem models for Martin
path_merge = 'D:\Users\djk2312\Documents\Harvard\reg\031\13\merge\';
if ~exist(path_merge,'dir');
    mkdir(path_merge);
end
options_showfig = true;
path_exp = 'D:\Users\djk2312\Documents\Harvard\reg\';
filename_csv = 'D:\Users\djk2312\Documents\thisshouldbecyclone\2012-08-Harvard.csv';
fid = fopen(filename_csv);
%C = textscan(fid, '%u8%u8%s%s%s%u16%u16%u16%u16%u16%u16%u16%u16%u16',...
%    'delimiter', ',',...
%    'headerlines',8);
C = textscan(fid, '%u8%u8%s%s%s%s%s%s%s%s%s%s%s%s',...
    'delimiter', ',',...
    'headerlines',8);
% 1 Site
% 2 Plot
% 3 Lidar Filename
% 4 Date
% 5 Lidar Orientation

%n_scans = numel(C{1});
warning('off', 'arguments:exteriordata');

site_unique = unique(C{1});
n_site = numel(site_unique);
info_site = 31;
is_site = site_unique==info_site;
plot_t = C{2}(C{1}==info_site);
n_plot = numel(plot_t);
info_plot = cell(n_plot,1);
for j = 1:n_plot;
    info_plot{j} = sprintf('%02.0f', j);
end
path_reg = sprintf('%s%03.0f%s', path_exp, info_site, '\25\mat\');
filepath_G_R_MST = sprintf('%sG_R_MST.mat',path_reg);
filepath_G_t_MST = sprintf('%sG_t_MST.mat',path_reg);
load(filepath_G_R_MST);
load(filepath_G_t_MST);

clear C fid filename_csv filepath_G_R_MST filepath_G_t_MST
clear j 
%% Load ply - DEM
%{
path_tree_out = sprintf('%s%03.0f%s', path_exp, info_site, '\13\ply\stems\');
path_dem_out = sprintf('%s%03.0f%s', path_exp, info_site, '\13\ply\dem\');
if ~exist(path_tree_out, 'dir');
    mkdir(path_tree_out);
end
if ~exist(path_dem_out, 'dir');
    mkdir(path_dem_out);
end
all_data_xJ = cell(n_plot,1);
all_data_yJ = cell(n_plot,1);
all_data_zJ = cell(n_plot,1);
all_nel = zeros(n_plot,1);
for j = 1:n_plot;
    path_ply = sprintf('%s%03.0f%s%s%s', path_exp, info_site, '\', info_plot{j},'\ply\');
    % Load DEM
    %
        filepath_dem = sprintf('%sply_dem_%03.0f-%s.ply',path_ply, info_site, info_plot{j});
        [vertexJ,faces] = read_ply(filepath_dem);
        data_xJ = vertexJ(:,1);
        data_yJ = vertexJ(:,2);
        data_zJ = vertexJ(:,3);
        all_nel(j) = numel(data_xJ);
        xyz2t = (G_R_MST{13}'*G_R_MST{j}*[data_xJ data_yJ data_zJ]') + ...
            repmat(...
            G_R_MST{13}'*G_t_MST{j} - ...
            G_R_MST{13}'*G_t_MST{13},...
            [1,numel(data_xJ)]);
        all_data_xJ{j} = xyz2t(1,:);
        all_data_yJ{j} = xyz2t(2,:);
        all_data_zJ{j} = xyz2t(3,:);
        % colorJ = vertexJ(:,4:6);
        % filepath_dem_out = sprintf('%sDtx-%s.ply', path_dem_out, info_plot{j});
        %colorout = colorJ.*repmat(all_color(j,:), [numel(data_xJ), 1]);
        %write2plyfaces(filepath_dem_out, xyz2t', colorJ,faces);
end
clear filepath_dem vertexJ faces data_xJ data_yJ data_zJ colorJ xyz2t k 
%}
%% Synthesize DEM
%{
% Figure showing all DEM
%{
figure;
hold on
color = jet(n_plot);
color = color(randperm(n_plot),:);
for j = 1:n_plot;
    scatter3(all_data_xJ{j},all_data_yJ{j},all_data_zJ{j},10,color(j,:), 'filled');
end
%}

n_synth = sum(all_nel);
synth_data_xJ = zeros(n_synth,1);
synth_data_yJ = zeros(n_synth,1);
synth_data_zJ = zeros(n_synth,1);
ixstart = 1;
for j = 1:n_plot;
    ixend = ixstart + all_nel(j) - 1;
    synth_data_xJ(ixstart:ixend) = all_data_xJ{j};
    synth_data_yJ(ixstart:ixend) = all_data_yJ{j};
    synth_data_zJ(ixstart:ixend) = all_data_zJ{j};
    ixstart = ixend + 1;
end
p_dem_step = 1;
p_dem_res = p_dem_step/2;

xmin = floor(min(synth_data_xJ));
xmax = ceil(max(synth_data_xJ));
ymin = floor(min(synth_data_yJ));
ymax = ceil(max(synth_data_yJ));

xstep = xmin:p_dem_step:xmax;
ystep = ymin:p_dem_step:ymax;
n_xstep = numel(xstep);
n_ystep = numel(ystep);

[ymesh, xmesh] = meshgrid(xstep,ystep); %Note change ?? why?
zmesh_min = nan(size(xmesh));
zmesh_med = nan(size(xmesh));
for x=1:n_xstep;
    isvalidx = (abs(synth_data_xJ-xstep(x))<p_dem_res);
    for y = 1:n_ystep;
        isvalidy = (abs(synth_data_yJ-ystep(y))<p_dem_res);
        minval = min(synth_data_zJ(isvalidx&isvalidy));
        medval = nanmedian(synth_data_zJ(isvalidx&isvalidy));
        if ~isempty(minval)
        zmesh_min(x,y) = minval;
        zmesh_med(x,y) = medval;
        end
    end
end

% Surf both meshes
%{
figure;
surf(xmesh, ymesh, zmesh_min, 'facecolor','r');
hold on
surf(xmesh, ymesh, zmesh_med,'facecolor', 'b');
%}

clear n_synth synth_data_xJ synth_data_yJ synth_data_zJ ixstart 
clear j ixend p_dem p_dem_res xmin ymin xmax ymax xstep ystep n_xstep 
clear n_ystep x isvalidx isvalidy minval medval 
%}
%% save ply 
filepath_demx = sprintf('%sdemx.mat',path_merge);
filepath_demy = sprintf('%sdemy.mat',path_merge);
filepath_demz = sprintf('%sdemz.mat',path_merge);
%{
filepath_ply_dem = 'D:\Users\djk2312\Documents\Harvard\reg\031\13\ply\dem\Dtx_merge.ply';
  isvalid = ~isnan(zmesh_med);
  demx = xmesh(isvalid);
  demy = ymesh(isvalid);
  demz = zmesh_med(isvalid);
  plycolor = vec2cmap2(demz, 'jet');
  dt = DelaunayTri(demx,demy);
  
  write2plyfaces_2( filepath_ply_dem,  ...
        [demx demy demz], ...
        plycolor, ...
        dt.Triangulation )
 

save(filepath_demx, 'demx');
save(filepath_demy, 'demy');
save(filepath_demz, 'demz');

clear filepath_ply_dem plycolor isvalid dt 
clear filepath_demx filepath_demy filepath_demz
%}
load(filepath_demx);
load(filepath_demy);
load(filepath_demz);
clear filepath_demx filepath_demy filepath_demz
%% Save transformed stem models PLY 
%{
for j = 1:n_plot;
    path_ply = sprintf('%s%03.0f%s%s%s', path_exp, info_site, '\', info_plot{j},'\ply\');
    % Load Stem Model
    filepath_tree = sprintf('%stree_%03.0f-%s.ply',path_ply, info_site, info_plot{j});
    filepath_tree_out = sprintf('%sStx-%s.ply', path_tree_out, info_plot{j});
    [vertexJ,faces] = read_ply(filepath_tree);
    data_xJ = vertexJ(:,1);
    data_yJ = vertexJ(:,2);
    data_zJ = vertexJ(:,3);
    colorJ = vertexJ(:,4:6);
    xyz2t = (G_R_MST{13}'*G_R_MST{j}*[data_xJ data_yJ data_zJ]') + ...
        repmat(...
        G_R_MST{13}'*G_t_MST{j} - ...
        G_R_MST{13}'*G_t_MST{13},...
        [1,numel(data_xJ)]);
    %colorout = colorJ.*repmat(all_color(j,:), [numel(data_xJ), 1]);
    write2plyfaces(filepath_tree_out, xyz2t', colorJ,faces);
    
end
clear data_xJ data_yJ data_zJ colorJ xyz2t faces vertexJ filepath_tree
clear filepath_tree_out path_ply j 
%}
%% Load tree and transform to WCS 
all_treetx = cell(n_plot,1);
n_all = 0;
allmin = inf;
allmax = -inf;
all_ntree = zeros(n_plot,1);
all_nel = cell(n_plot,1);
for j = 1:n_plot ;
    path_mat = sprintf('%s%03.0f%s%s%s', path_exp, info_site, '\', info_plot{j},'\mat\');
    filepath_tree = sprintf('%stree.mat',path_mat);
    %filepath_treetx = sprintf('%streetx.mat',path_mat);
    load(filepath_tree);
    n_tree = numel(tree);
    all_ntree(j) = n_tree;
    n_all = n_all + n_tree;
    all_nel{j} = zeros(n_tree,1);
    for t = 1:n_tree;
        n_el = size(tree(t).loc,2);
        all_nel{j}(t) = n_el;
        all_treetx{j}(t).loc = (G_R_MST{13}'*G_R_MST{j}*tree(t).loc) + ...
            repmat(...
            G_R_MST{13}'*G_t_MST{j} - ...
            G_R_MST{13}'*G_t_MST{13},...
            [1,n_el]);
        all_treetx{j}(t).r = tree(t).r;
        if all_treetx{j}(t).loc(3,1) < allmin;
            allmin = all_treetx{j}(t).loc(3,1);
        end
        if all_treetx{j}(t).loc(3,end) > allmax;
            allmax = all_treetx{j}(t).loc(3,end);
        end
    end
end

clear path_mat filepath_tree filepath_treetx n_tree t j n_el tree
% Linear interpolation sanity check
%{
zstep = .1;
zq = allmin:zstep:allmax;
z = treetx{1}(1).loc(3,:);
vx = treetx{1}(1).loc(1,:);
vy = treetx{1}(1).loc(2,:);
xq = interp1(z,vx,zq);
yq = interp1(z,vy,zq);
%}
%% Combine all_treetx into array 
zstep = .1;
zq = allmin:zstep:allmax;
n_z = numel(zq);
tx = nan(n_all,n_z);
ty = nan(n_all,n_z);
tr = nan(n_all,n_z);
ix = 0;
for j = 1:n_plot;
    for t = 1:all_ntree(j);
        ix = ix + 1;
        z = all_treetx{j}(t).loc(3,:);
        vx = all_treetx{j}(t).loc(1,:);
        vy = all_treetx{j}(t).loc(2,:);
        vr = all_treetx{j}(t).r;
        tx(ix,:) = interp1(z,vx,zq);
        ty(ix,:) = interp1(z,vy,zq);
        tr(ix,:) = interp1(z,vr,zq);
    end
end
tz = repmat(zq, [n_all,1]);

% Figure - plot stem models 
%{
foo = 1;
tc = jet(n_all);
tc = tc(randperm(n_all),:);
figure;
hold on
for a = 1:n_all;
scatter3(tx(a,:), ty(a,:), tz(a,:), 10, tc(a,:), 'filled');
end
%}
clear vx vy vr z ix j t all_treetx 
%% Remove stem models which do not overlap with another t_con=95% 

% Determine number of segments in each stem
tnnan = ~isnan(tx);
tnseg = sum(tnnan,2);
tnconn = zeros(n_all,1);
%t_overlap = 50;
for a = 1:n_all;
    dist = sqrt((repmat(tx(a,:), [n_all,1]) - tx).^2 + ...
        (repmat(ty(a,:), [n_all,1]) - ty).^2);
    rad = repmat(tr(a,:), [n_all,1]) + tr;
    isvalid = dist<rad;
    isvalid(a,:) = false;
    tnconn(a) = sum(any(isvalid,1));
end

tpconn = tnconn./tnseg;
% Figure - plot stem models with color = pconn
%{
foo = 1;
tc = double(vec2cmap2(tpconn, 'jet'))./255;
figure;
hold on
for a = 1:n_all;
scatter3(tx(a,:), ty(a,:), tz(a,:), 10, tc(a,:), 'filled');
end
figure;
hold on
for a = 1:n_all;
scatter(tpconn(a),0,20,tc(a,:), 'filled');
end
%}
t_conn = 0.95;
isvalid = tpconn > t_conn;
tx = tx(isvalid,:);
ty = ty(isvalid,:);
tz = tz(isvalid,:);
tr = tr(isvalid,:);
tnseg = tnseg(isvalid);
n_all = size(tx,1);
%{
tc = tc(isvalid,:);
figure;
hold on
for a = 1:n_all;
scatter3(tx(a,:), ty(a,:), tz(a,:), 10, tc(a,:), 'filled');
end
%}

clear dist rad isvalid tnconn tnnan
%% Find correspondences 

isnnan = ~isnan(tx);
n_total = sum(isnnan,2);
used = false(n_all,1);
sx = nan(n_all,n_z);
sy = nan(n_all,n_z);
sz = repmat(zq, [n_all,1]);
sr = nan(n_all,n_z);
for a = 1:n_all;
    if used(a);
        continue
    end
    dist = sqrt((repmat(tx(a,:), [n_all,1]) - tx).^2 + ...
        (repmat(ty(a,:), [n_all,1]) - ty).^2);
    rad = repmat(tr(a,:), [n_all,1]) + tr;
    rad = max(rad, 0.25);
    isvalid_dist = dist<rad;
    tree_isdist = any(isvalid_dist,2);
    tree_poverlap = sum(isvalid_dist,2)./n_total;
    tree_ixpoverlap = find(tree_poverlap > 0.8);
    if numel(tree_ixpoverlap) ==1;
        continue
    end
    used(tree_isdist) = true;
    seg_nmatch = sum(isvalid_dist(tree_ixpoverlap,:),1);
    seg_ismatch = seg_nmatch > 1;
    sx(a,:) = nanmedian(tx(tree_ixpoverlap,:),1);
    sy(a,:) = nanmedian(ty(tree_ixpoverlap,:),1);
    sr(a,:) = nanmedian(tr(tree_ixpoverlap,:),1);
    sx(a,~seg_ismatch) = nan;
    sy(a,~seg_ismatch) = nan;
    sr(a,~seg_ismatch) = nan;
    %{
    figure;
    x = tx(ix_used,:);
    y = ty(ix_used,:);
    z = repmat(zq, [numel(ix_used), 1]);
    scatter3(x(:), y(:), z(:), 10,'r','filled')
    hold on
    scatter3(sx(a,:), sy(a,:), sz(a,:), 10,'b','filled')
    axis equal
    foo = 1;
    %}
end

isvalid = ~all(isnan(sx),2);
sx = sx(isvalid,:);
sy = sy(isvalid,:);
sz = sz(isvalid,:);
sr = sr(isvalid,:);
n_s = size(sx,1);
% Figure plot stem models 
%{
c = jet(n_s);
c = c(randperm(n_s),:);
figure 
hold on;
for s = 1:n_s;
    scatter3(sx(s,:), sy(s,:), sz(s,:), 10, c(s,:), 'filled');
end
%}
clear isnnan n_total used a dist rad isvalid v_isvalid n_valid p_valid 
clear ix_used x y z c s 
%% Extend to DEM and height 
%{
firstzero = (cumsum(~isnan(sx),2)==1);
[~, simin]= max( firstzero~=0, [], 2 );
sxmin = nan(n_s,1);
symin = nan(n_s,1);
szmin = nan(n_s,1);
for s = 1:n_s;
sxmin(s) = sx(s,simin(s));
symin(s) = sy(s,simin(s));
szmin(s) = sz(s,simin(s));
end

sidem = nan(n_s,1);
for s = 1:n_s;
    [~,demix] = min(sqrt((demx - sxmin(s)).^2 + (demy - symin(s)).^2));
    if demz(demix) < szmin(s);
        [~,ix] = min(abs(zq - demz(demix)));
        ix_dem = ix-1;
        sidem(s) = ix_dem;
        sz(s,ix_dem:simin(s)+5) = zq(ix_dem:simin(s)+5);
        dx = (sx(s,simin(s)+5) - sx(s,simin(s)))./5;
        dy = (sy(s,simin(s)+5) - sy(s,simin(s)))./5;
        dr = 1+(-((sr(s,simin(s)+5) - sr(s,simin(s)))./sr(s,simin(s)))./5);
        n_step = simin(s) - ix_dem; 
        % original 
        sx(s,ix_dem) = sx(s,simin(s))-n_step*dx;
        sy(s,ix_dem) = sy(s,simin(s))-n_step*dy;
        sr(s,ix_dem) = sr(s,simin(s))*n_step*dr;
        % vector 
        stepval = fliplr(1:numel(ix_dem:simin(s)+5));
        sx(s,ix_dem:simin(s)+5) = sx(s,simin(s))-stepval*dx;
        sy(s,ix_dem:simin(s)+5) = sy(s,simin(s))-stepval*dy;
        sr(s,ix_dem:simin(s)+5) = repmat(sr(s,simin(s)), [numel(ix_dem:simin(s)+5), 1]);    
        %{
        figure; 
        hold on
        scatter3(sx(s,:), sy(s,:), sz(s,:), 10,'r', 'filled');
        scatter3(sx(s,ix_dem), sy(s,ix_dem), sz(s,ix_dem), 40,'b', 'filled');
        view(0,0);
        axis equal;
        foo = 1;
        %}
    end   
end
%}
%% Make Tree structure 
cmap = jet(n_s);
cmap = double(uint8(255*cmap));
cmap = cmap(randperm(n_s),:);
for s = 1:n_s;
    isnnan = ~isnan(sx(s,:));
    tree(s).loc = [sx(s,isnnan); sy(s,isnnan); sz(s,isnnan)];
    tree(s).r = sr(s,isnnan);
    tree(s).color = cmap(s,:)';
end

%% Interpolate 
n_tree = n_s;
t_z_resample = zstep*4;%.25;
t_z_interp = zstep*2;%0.125;
%t_z_resample = 1;
%t_z_interp = 0.25;
counter = 1;
for t = 1:n_tree;   
     minz = tree(t).loc(3,1);
     maxz = tree(t).loc(3,end);
     if t == 23;
         foo = 1;
     end
     zstep = minz:t_z_resample:maxz;%+t_z_resample;
     if numel(zstep)==1;
         continue
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

%% Save PLY 
%path_local = 'Z:\Desktop\Local\';
filepath_tree_ply = 'D:\Users\djk2312\Documents\Harvard\reg\031\13\ply\stems\Smerge_sm.ply';
tree2ply(filepath_tree_ply,tree_sm,20);
foo = 1;
%% Write CSV 
filepath_csv =  'D:\Users\djk2312\Documents\Harvard\reg\031\13\ply\stems\tree.csv';
fid = fopen(filepath_csv, 'w+');
fprintf(fid, 't,x,y,z,d\n');
for s = 1:n_s;
    for e = 1:numel(tree_sm(s).r);
    fprintf(fid, '%u, %f,%f,%f,%f\n', s, tree_sm(s).loc(1,e),tree_sm(s).loc(2,e),...
        tree_sm(s).loc(3,e), 2*tree_sm(s).r(e));
    end
end
fclose all




