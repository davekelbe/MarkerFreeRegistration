% Unify stem models for Martin

options_loaddata = true;
if ~options_loaddata;
    
    %%
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
    site_unique_orig = unique(C{1});
    
    n_site = numel(site_unique);
    
    %%
    % Load field data
    
    info_exp = 'harvard';
    info_suffix = 'reg';
    info_slash = '\';
    % Make directories
    path_common = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, 'Common',info_slash);
    
    NEON_x = cell(n_site,1);
    NEON_y = cell(n_site,1);
    NEON_dbho = cell(n_site,1);
    NEON_dbhu = cell(n_site,1);
    
    for six =1:n_site;
        fprintf('Site %u of %u\n', six, n_site);
        info_site = site_unique(six);
        info_plot = 25;
        path_top = sprintf('%s%s%s%s%s%03.0f%s%02.0f%s','D:\Users\djk2312\Documents\',...
            info_exp, info_slash, info_suffix,info_slash,info_site,info_slash,info_plot,info_slash);
        path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
        
        % Parse Field Data
        filepath_NEON_DBH_10over_csv = sprintf('%sField_Data\\NEON_DBH_10over.csv',path_common);
        filepath_NEON_DBH_10over = sprintf('%sNEON_DBH_10over.mat',path_mat);
        if exist(filepath_NEON_DBH_10over,'file');
            load(filepath_NEON_DBH_10over);
        elseif ~exist(filepath_NEON_DBH_10over);
            NEON_DBH_10over = parse_NEON_DBH_10over(filepath_NEON_DBH_10over_csv,info_site);
            save(filepath_NEON_DBH_10over,'NEON_DBH_10over');
        end
        
        NEON_x{six} = NEON_DBH_10over.x;
        NEON_y{six} = NEON_DBH_10over.y;
        NEON_dbho{six} = NEON_DBH_10over.dbhm;
        
        filepath_NEON_DBH_10under_csv = sprintf('%sField_Data\\NEON_DBH_10under.csv',path_common);
        filepath_NEON_DBH_10under = sprintf('%sNEON_DBH_10under.mat',path_mat);
        if exist(filepath_NEON_DBH_10under,'file');
            load(filepath_NEON_DBH_10under);
        elseif ~exist(filepath_NEON_DBH_10under);
            NEON_DBH_10under = parse_NEON_DBH_10under(filepath_NEON_DBH_10under_csv,info_site);
            save(filepath_NEON_DBH_10under,'NEON_DBH_10under');
        end
        
        NEON_dbhu{six} = NEON_DBH_10under.dbh/100;
    end
    %% Collate plot level parameters
    
    NEON_BAo = zeros(n_site,1);
    NEON_BAu = zeros(n_site,1);
    NEON_SHo = zeros(n_site,1);
    NEON_SHu = zeros(n_site,1);
    for six = 1:n_site;
        isvalid = NEON_x{six} >= -10 & NEON_x{six} <= 10 & NEON_y{six} >= -10 & NEON_y{six} <= 10;
        NEON_BAo(six) = 10000*nansum(pi*(NEON_dbho{six}(isvalid)/2).^2)./(20^2);
        NEON_BAu(six) = 10000*nansum(pi*(NEON_dbhu{six}/2).^2)./(20^2);
        NEON_SHo(six) = sum(isvalid)*10000/(20^2);
        NEON_SHu(six) = numel(NEON_dbhu{six})*10000/(20^2);
    end
    NEON_BA = NEON_BAo + NEON_BAu;
    NEON_SH = NEON_SHo + NEON_SHu;
    
    %%  Synthesize DEM
    %{
for six =1:n_site;
    fprintf('Site %u of %u\n', six, n_site);
    info_site = site_unique(six);
    path_merge = sprintf('%s%03.0f%s','D:\Users\djk2312\Documents\Harvard\reg\',info_site, '\13\merge\');
    if ~exist(path_merge,'dir');
        mkdir(path_merge);
    end
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
    
    clear fid filename_csv filepath_G_R_MST filepath_G_t_MST
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
        if any([isempty(G_R_MST{j}), isempty(G_R_MST{13})]);
            continue
        end
        path_ply = sprintf('%s%03.0f%s%s%s', path_exp, info_site, '\', info_plot{j},'\ply\');
        % Load DEM
        %
        filepath_dem = sprintf('%sply_dem_%03.0f-%s.ply',path_ply, info_site, info_plot{j});
        [vertexJ,~] = read_ply(filepath_dem);
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
    %
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
    
    %% save ply
    filepath_demx = sprintf('%sdemx.mat',path_merge);
    filepath_demy = sprintf('%sdemy.mat',path_merge);
    filepath_demz = sprintf('%sdemz.mat',path_merge);
    %
    filepath_ply_dem = sprintf('%s%03.0f%s','D:\Users\djk2312\Documents\Harvard\reg\', info_site, '\13\ply\dem\Dtx_merge.ply');
    isvalid = ~isnan(zmesh_med);
    demx = xmesh(isvalid);
    demy = ymesh(isvalid);
    demz = zmesh_med(isvalid);
    if isempty(demz);
        continue
    end
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
%    clear filepath_demx filepath_demy filepath_demz
    
    load(filepath_demx);
    load(filepath_demy);
    load(filepath_demz);
    clear filepath_demx filepath_demy filepath_demz
    %}
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
    isvalid = true(n_site,1);
    
    all_treetx = cell(n_site,1);
    allmin = inf(n_site,1);
    allmax = -inf(n_site,1);
    n_all = zeros(n_site,1);
    all_ntree = cell(n_site,1);
    for six =1:n_site;
        %         if site_unique(six)~= 4;
        %             continue
        %         end
        fprintf('Site %u of %u\n', six, n_site);
        info_site = site_unique(six);
        n_plot = 25;
        all_treetx{six} = cell(n_plot,1);
        all_ntree{six} = zeros(n_plot,1);
        all_nel = cell(n_plot,1);
        path_reg = sprintf('%s%03.0f%s', path_exp, info_site, '\25\mat\');
        filepath_G_R_WMF = sprintf('%sG_R_WMF.mat',path_reg);
        filepath_G_t_WMF = sprintf('%sG_t_WMF.mat',path_reg);
        if ~exist(filepath_G_R_WMF, 'file');
            isvalid(six) = false;
            continue
        end
        load(filepath_G_R_WMF);
        load(filepath_G_t_WMF);
        all_treetx{six} = cell(n_plot,1);
        for j = 1:n_plot ;
            path_mat = sprintf('%s%03.0f%s%02.0f%s', path_exp, info_site, '\', j,'\mat\');
            filepath_tree = sprintf('%stree.mat',path_mat);
            %filepath_treetx = sprintf('%streetx.mat',path_mat);
            load(filepath_tree);
            n_tree = numel(tree);
            all_ntree{six}(j) = n_tree;
            n_all(six) = n_all(six) + n_tree;
            all_nel{six}{j} = zeros(n_tree,1);
            if any(any(isnan(G_R_WMF{j})));
                foo = 1;
                continue
            end
            if any(any(isnan(G_R_WMF{13})));
                G_R_WMF{13} = eye(3);
                G_t_WMF{13} = zeros(3,1);
            end
            for t = 1:n_tree;
                n_el = size(tree(t).loc,2);
                all_nel{six}{j}(t) = n_el;
                all_treetx{six}{j}(t).loc = (G_R_WMF{13}'*G_R_WMF{j}*tree(t).loc) + ...
                    repmat(...
                    G_R_WMF{13}'*G_t_WMF{j} - ...
                    G_R_WMF{13}'*G_t_WMF{13},...
                    [1,n_el]);
                all_treetx{six}{j}(t).r = tree(t).r;
                if all_treetx{six}{j}(t).loc(3,1) < allmin(six);
                    allmin(six) = all_treetx{six}{j}(t).loc(3,1);
                end
                if all_treetx{six}{j}(t).loc(3,end) > allmax(six);
                    allmax(six) = all_treetx{six}{j}(t).loc(3,end);
                end
            end
        end
    end
    all_treetx = all_treetx(isvalid);
    allmin = allmin(isvalid);
    allmax = allmax(isvalid);
    site_unique = site_unique(isvalid);
    n_all = n_all(isvalid);
    all_ntree = all_ntree(isvalid);
    isvalid = ~isinf(allmin);
    all_treetx = all_treetx(isvalid);
    allmin = allmin(isvalid);
    allmax = allmax(isvalid);
    site_unique = site_unique(isvalid);
    n_all = n_all(isvalid);
    all_ntree = all_ntree(isvalid);
    isvalid = site_unique<100;
    all_treetx = all_treetx(isvalid);
    allmin = allmin(isvalid);
    allmax = allmax(isvalid);
    site_unique = site_unique(isvalid);
    n_all = n_all(isvalid);
    all_ntree = all_ntree(isvalid);
    n_site = numel(all_treetx);
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
    tx = cell(n_site,1);
    ty = cell(n_site,1);
    tr = cell(n_site,1);
    for six = 1:n_site;
        zstep = .1;
        zq = allmin(six):zstep:allmax(six);
        n_z = numel(zq);
        tx{six} = nan(n_all(six),n_z);
        ty{six} = nan(n_all(six),n_z);
        tr{six} = nan(n_all(six),n_z);
        ix = 0;
        for j = 1:n_plot;
            for t = 1:all_ntree{six}(j);
                if isempty(all_treetx{six}{j});
                    continue
                end
                ix = ix + 1;
                z = all_treetx{six}{j}(t).loc(3,:);
                vx = all_treetx{six}{j}(t).loc(1,:);
                vy = all_treetx{six}{j}(t).loc(2,:);
                vr = all_treetx{six}{j}(t).r;
                tx{six}(ix,:) = interp1(z,vx,zq);
                ty{six}(ix,:) = interp1(z,vy,zq);
                tr{six}(ix,:) = interp1(z,vr,zq);
            end
        end
        tz{six} = repmat(zq, [n_all(six),1]);
    end
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
    % else
    %     filepath_line418 = sprintf('%sline418.mat', path_common);
    %     if ~exist(filepath_line418, 'file');
    %         save(filepath_line418);
    %     else
    %         load(filepath_line418);
    %     end
    % end
    %% Remove stem models which do not overlap with another t_con=95%
    
    % Determine number of segments in each stem
    tnseg = cell(n_site,1);
    for six = 1:n_site;
        tnnan = ~isnan(tx{six});
        tnseg{six} = sum(tnnan,2);
        tnconn = zeros(n_all(six),1);
        %t_overlap = 50;
        for a = 1:n_all(six);
            dist = sqrt((repmat(tx{six}(a,:), [n_all(six),1]) - tx{six}).^2 + ...
                (repmat(ty{six}(a,:), [n_all(six),1]) - ty{six}).^2);
            rad = repmat(tr{six}(a,:), [n_all(six),1]) + tr{six};
            isvalid = dist<rad;
            isvalid(a,:) = false;
            tnconn(a) = sum(any(isvalid,1));
        end
        
        tpconn = tnconn./tnseg{six};
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
        tx{six} = tx{six}(isvalid,:);
        ty{six} = ty{six}(isvalid,:);
        tz{six} = tz{six}(isvalid,:);
        tr{six} = tr{six}(isvalid,:);
        tnseg{six} = tnseg{six}(isvalid);
        n_all(six) = size(tx{six},1);
        %{
tc = tc(isvalid,:);
figure;
hold on
for a = 1:n_all;
scatter3(tx(a,:), ty(a,:), tz(a,:), 10, tc(a,:), 'filled');
end
        %}
    end
    
    clear dist rad isvalid tnconn tnnan
else
    filepath_line481 = sprintf('%sline481.mat', path_common);
    if ~exist(filepath_line481, 'file');
        save(filepath_line481);
    else
        load(filepath_line481);
    end
end

n_z = cellfun(@(x) size(x,2), tx);

isvalid = ismember(site_unique_orig, site_unique);
NEON_BAo = NEON_BAo(isvalid);
NEON_BAu = NEON_BAu(isvalid);
NEON_SHo = NEON_SHo(isvalid);
NEON_SHu = NEON_SHu(isvalid);
NEON_BA = NEON_BA(isvalid);
NEON_SH = NEON_SH(isvalid);



rad_arr = 0.25:0.25:2.5;
n_rad = numel(rad_arr);
pover_arr = 0.5:.1:3;
n_pover = numel(pover_arr);

E_BAo = nan(n_rad,n_pover);
E_SHo = nan(n_rad,n_pover);
E_BA = nan(n_rad,n_pover);
E_SH = nan(n_rad,n_pover);

for rad_ix = 1:n_rad;
    trad = rad_arr(rad_ix);
    for pover_ix = 1:n_pover;
        tpover = pover_arr(pover_ix);
        %% Find correspondences
        sx = cell(n_site,1);
        sy = cell(n_site,1);
        sz = cell(n_site,1);
        sr = cell(n_site,1);
        n_s = zeros(n_site,1);
        for six = 1:n_site;
            fprintf('Site %u of %u\n', six, n_site);
            isnnan = ~isnan(tx{six});
            n_total = sum(isnnan,2);
            used = false(n_all(six),1);
            sx{six} = nan(n_all(six),n_z(six));
            sy{six} = nan(n_all(six),n_z(six));
            sz{six} = repmat(zq, [n_all(six),1]);
            sr{six} = nan(n_all(six),n_z(six));
            for a = 1:n_all(six);
                if used(a);
                    continue
                end
                dist = sqrt((repmat(tx{six}(a,:), [n_all(six),1]) - tx{six}).^2 + ...
                    (repmat(ty{six}(a,:), [n_all(six),1]) - ty{six}).^2);
                rad = repmat(tr{six}(a,:), [n_all(six),1]) + tr{six};
                rad = max(rad, trad);
                isvalid_dist = dist<rad;
                tree_isdist = any(isvalid_dist,2);
                tree_poverlap = sum(isvalid_dist,2)./n_total;
                tree_ixpoverlap = find(tree_poverlap > tpover);
                if numel(tree_ixpoverlap) ==1;
                    continue
                end
                used(tree_isdist) = true;
                seg_nmatch = sum(isvalid_dist(tree_ixpoverlap,:),1);
                seg_ismatch = seg_nmatch > 1;
                sx{six}(a,:) = nanmedian(tx{six}(tree_ixpoverlap,:),1);
                sy{six}(a,:) = nanmedian(ty{six}(tree_ixpoverlap,:),1);
                sr{six}(a,:) = nanmedian(tr{six}(tree_ixpoverlap,:),1);
                sx{six}(a,~seg_ismatch) = nan;
                sy{six}(a,~seg_ismatch) = nan;
                sr{six}(a,~seg_ismatch) = nan;
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
            
            isvalid = ~all(isnan(sx{six}),2);
            sx{six} = sx{six}(isvalid,:);
            sy{six} = sy{six}(isvalid,:);
            sz{six} = sz{six}(isvalid,:);
            sr{six} = sr{six}(isvalid,:);
            n_s(six) = size(sx{six},1);
        end
        
        % Compute DBH
        TLS_x = cellfun(@(x) nanmean(x,2),sx,'uniformoutput',false);
        TLS_y = cellfun(@(x) nanmean(x,2),sx,'uniformoutput',false);
        TLS_DBH = cellfun(@(x) 2*nanmean(x,2),sr,'uniformoutput',false);
        
        % Calculated basal area and stem density
        TLS_BAo = zeros(n_site,1);
        TLS_BAu = zeros(n_site,1);
        TLS_SHo = zeros(n_site,1);
        TLS_SHu = zeros(n_site,1);
        for six = 1:n_site;
            isvalid = TLS_x{six} >= -10 & TLS_x{six} <= 10 & TLS_y{six} >= -10 & TLS_y{six} <= 10;
            isover = TLS_DBH{six} >= .1;
            TLS_BAo(six) = 10000*nansum(pi*(TLS_DBH{six}(isvalid&isover)/2).^2)./(20^2);
            TLS_BAu(six) = 10000*nansum(pi*(TLS_DBH{six}(isvalid&~isover)/2).^2)./(20^2);
            TLS_SHo(six) = sum(isvalid&isover)*10000/(20^2);
            TLS_SHu(six) = sum(isvalid&~isover)*10000/(20^2);
        end
        TLS_BA = TLS_BAo + TLS_BAu;
        TLS_SH = TLS_SHo + TLS_SHu;
        
        E_BAo(rad_ix, pover_ix) = sqrt(mean(NEON_BAo-TLS_BAo).^2);
        E_BA(rad_ix, pover_ix) = sqrt(mean(NEON_BA-TLS_BA).^2);
        E_SHo(rad_ix, pover_ix) = sqrt(mean(NEON_SHo-TLS_SHo).^2);
        E_SH(rad_ix, pover_ix) = sqrt(mean(NEON_SH-TLS_SH).^2);
        
    end
end

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
% Martin stuff
%{
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
    filepath_tree_ply = sprintf('%sSmerge_sm.ply', path_merge);
    tree2ply(filepath_tree_ply,tree_sm,20);
    foo = 1;
    %% Write CSV
    filepath_tree_csv = sprintf('%stree.csv', path_merge);
    %filepath_csv =  'D:\Users\djk2312\Documents\Harvard\reg\031\13\ply\stems\tree.csv';
    fid = fopen(filepath_tree_csv, 'w+');
    fprintf(fid, 't,x,y,z,d\n');
    for s = 1:n_s;
        for e = 1:numel(tree_sm(s).r);
            fprintf(fid, '%u, %f,%f,%f,%f\n', s, tree_sm(s).loc(1,e),tree_sm(s).loc(2,e),...
                tree_sm(s).loc(3,e), 2*tree_sm(s).r(e));
        end
    end
    fclose all
%}


