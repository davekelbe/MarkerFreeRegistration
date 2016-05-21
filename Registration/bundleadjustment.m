
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


