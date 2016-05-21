%batch grah results
options_showfig = true;
aux.info_exp = 'reg';
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

n_scans = numel(C{1});
warning('off', 'arguments:exteriordata');

site_unique = unique(C{1});
n_site = numel(site_unique);

% Load graph ouputs
% MST
all_R_MST = cell(n_site,1);
all_t_MST = cell(n_site,1);
all_RMSE_MST = cell(n_site,1);
% WMF
all_R_WMF = cell(n_site,1);
all_t_WMF = cell(n_site,1);
all_RMSE_WMF = cell(n_site,1);
% SVD
all_R_SVD = cell(n_site,1);
all_t_SVD = cell(n_site,1);
all_RMSE_SVD = cell(n_site,1);
% Pair
all_R_pair = cell(n_site,1);
all_t_pair = cell(n_site,1);
all_RMSE_pair = cell(n_site,1);
% P_LCS
all_P_LCS = cell(n_site,1);
all_P_rad = cell(n_site,1);
all_ntree = cell(n_site,1);
all_w = zeros(n_site,1);
for site = 1:n_site;
    info_site = site_unique(site);
    is_site = (C{1}==info_site);
    plot_t = C{2}(is_site);
    n_plot = numel(plot_t);
    info_valid_plot = cell(n_plot,1);
    for p = 1:n_plot;
        info_valid_plot{p} = sprintf('%02.0f', p);
    end
    fprintf('\nSite %d\n', info_site)
    path_mat = sprintf('%s%03.0f%s', path_exp, info_site, '\25\mat\');
    % MST
    filepath_R_MST = sprintf('%s%s',path_mat,'G_R_MST.mat');
    filepath_t_MST = sprintf('%s%s',path_mat,'G_t_MST.mat');
    filepath_RMSE_MST = sprintf('%s%s',path_mat,'G_RMSE_MST.mat');
    if ~exist(filepath_R_MST,'file');
        continue
    end
    G_R_MST = [];
    G_t_MST = [];
    G_RMSE_MST = [];
    load(filepath_R_MST);
    load(filepath_t_MST);
    load(filepath_RMSE_MST);
    all_R_MST{site} = G_R_MST;
    all_t_MST{site} = G_t_MST;
    all_RMSE_MST{site} = G_RMSE_MST;
    % WMF
    filepath_R_WMF = sprintf('%s%s',path_mat,'G_R_WMF.mat');
    filepath_t_WMF = sprintf('%s%s',path_mat,'G_t_WMF.mat');
    filepath_RMSE_WMF = sprintf('%s%s',path_mat,'G_RMSE_WMF.mat');
    if ~exist(filepath_R_WMF,'file');
        continue
    end
    G_R_WMF = [];
    G_t_WMF = [];
    G_RMSE_WMF = [];
    load(filepath_R_WMF);
    load(filepath_t_WMF);
    load(filepath_RMSE_WMF);
    all_R_WMF{site} = G_R_WMF;
    all_t_WMF{site} = G_t_WMF;
    all_RMSE_WMF{site} = G_RMSE_WMF;
    % SVD
    filepath_R_SVD = sprintf('%s%s',path_mat,'G_R_SVD.mat');
    filepath_t_SVD = sprintf('%s%s',path_mat,'G_t_SVD.mat');
    filepath_RMSE_SVD = sprintf('%s%s',path_mat,'G_RMSE_SVD.mat');
    if ~exist(filepath_R_SVD,'file');
        continue
    end
    G_R_SVD = [];
    G_t_SVD = [];
    G_RMSE_SVD = [];
    load(filepath_R_SVD);
    load(filepath_t_SVD);
    load(filepath_RMSE_SVD);
    all_R_SVD{site} = G_R_SVD;
    all_t_SVD{site} = G_t_SVD;
    all_RMSE_SVD{site} = G_RMSE_SVD;
    % Reference node
    filepath_w = sprintf('%s%s',path_mat, 'w.mat');
    load(filepath_w);
    all_w(site)= w;
    % Pair
    filepath_match1_R = sprintf('%s%s',path_mat, 'match1_R.mat');
    filepath_match_R = sprintf('%s%s',path_mat, 'match_R.mat');
    if exist(filepath_match_R,'file')
        load(filepath_match_R);
    else
        load(filepath_match1_R);
        match_R = match1_R;
    end
    filepath_match1_t = sprintf('%s%s',path_mat, 'match1_t.mat');
    filepath_match_t = sprintf('%s%s',path_mat, 'match_t.mat');
    if exist(filepath_match_t,'file')
        load(filepath_match_t);
    else
        load(filepath_match1_t);
        match_t = match1_t;
    end
    filepath_match12_RMSE = sprintf('%s%s',path_mat, 'match12_RMSE.mat');
    load(filepath_match12_RMSE);
    all_R_pair{site} = match_R(w,:)';
    all_t_pair{site} = match_t(w,:)';
    all_RMSE_pair{site} = match12_RMSE(w,:)';
    % P_LCS
    filepath_P_LCS = sprintf('%s%03.0f%s%02.0f%s%s', path_exp, info_site, ...
        '\', plot_t(p),'\mat\','P_LCS.mat');
    filepath_P_rad = sprintf('%s%03.0f%s%02.0f%s%s', path_exp, info_site, ...
        '\', plot_t(p),'\mat\','P_rad.mat');
    load(filepath_P_LCS);
    load(filepath_P_rad);
    all_P_LCS{site} = P_LCS;
    all_P_rad{site} = P_rad;
end

%% Graph Reported Error vs Pairwise Reported Error
%{
figure;
title('Real data');
hold on
scatter(cell2mat(all_RMSE_pair),cell2mat(all_RMSE_MST),'.r')
scatter(cell2mat(all_RMSE_pair),cell2mat(all_RMSE_WMF),'.g')
scatter(cell2mat(all_RMSE_pair),cell2mat(all_RMSE_SVD),'.b')
plot([0 1],[0 1], '-k');
axis([0 .5 0 .5]);
xlabel('Reported Pairwise Error [m]');
ylabel('Reported Graph Error [m]');
legend_str = {'MST','WMF','SVD','1:1'};
legend(legend_str);
%}
%% Percent Connected Graph vs Pair
foo = 1;
all_num_pair = zeros(n_site,1);
all_num_MST = zeros(n_site,1);
all_num_WMF = zeros(n_site,1);
all_num_SVD = zeros(n_site,1);
all_num_num = zeros(n_site,1);
for site=1:n_site;
    temp_pair = all_RMSE_pair{site};
    temp_MST = all_RMSE_MST{site};
    temp_WMF = all_RMSE_WMF{site};
    temp_SVD = all_RMSE_SVD{site};
    all_num_pair(site) =  sum((~isnan(temp_pair) & temp_pair<1));
    all_num_MST(site) =  sum((~isnan(temp_MST) & temp_MST<1));
    all_num_WMF(site) =  sum((~isnan(temp_WMF) & temp_WMF<1));
    all_num_SVD(site) =  sum((~isnan(temp_SVD) & temp_SVD<1));
    all_num_num(site) = numel(temp_pair);
end
all_perc_pair = 100*all_num_pair./all_num_num;
all_perc_MST = 100*all_num_MST./all_num_num;
all_perc_WMF = 100*all_num_WMF./all_num_num;
all_perc_SVD = 100*all_num_SVD./all_num_num;
%
figure;
hold on
scatter(all_perc_pair(:), all_perc_MST(:), '.r')
scatter(all_perc_pair(:), all_perc_WMF(:), '.g')
scatter(all_perc_pair(:), all_perc_SVD(:), '.b')
plot([0 100],[0 100],'-k');
axis([0 100 0 100]);
xlabel('Percent detected Pairwise');
ylabel('Percent detected Graphical');
%}

isvalid = ~isnan(all_perc_pair);
[site_unique(isvalid) all_perc_pair(isvalid) all_perc_MST(isvalid) ...
    all_perc_WMF(isvalid) all_perc_SVD(isvalid)];

[site_unique(isvalid) all_perc_pair(isvalid)/4 all_perc_MST(isvalid)/4 ...
    all_perc_WMF(isvalid)/4 all_perc_SVD(isvalid)/4]


%table_RMSE_WMF = 100*cellfun(@nanmean, all_RMSE_WMF(isvalid));
%table_RMSE_pair = 100*cellfun(@nanmean, all_RMSE_pair(isvalid));
RMSE_WMF_valid = all_RMSE_WMF(isvalid);
RMSE_pair_valid = all_RMSE_pair(isvalid);
n_valid = numel(RMSE_WMF_valid);
table_RMSE_WMF = zeros(n_valid,1);
table_RMSE_pair = zeros(n_valid,1);
for i = 1:n_valid;
    isvalidtemp = (RMSE_WMF_valid{i}<1);
    table_RMSE_WMF(i) = 100*nanmean(RMSE_WMF_valid{i}(isvalidtemp));
    isvalidtemp = (RMSE_pair_valid{i}<1);
    table_RMSE_pair(i) = 100*nanmean(RMSE_pair_valid{i}(isvalidtemp));    
end
%[site_unique(isvalid)  table_RMSE_WMF]
isvalid(end-1:end) = false;
table_site = site_unique(isvalid);
table_nscan = all_num_num(isvalid);
table_pair = all_num_pair(isvalid);
table_pairp = 100*table_pair./table_nscan;
table_WMF = all_num_WMF(isvalid);
table_WMFp = 100*table_WMF./table_nscan;
n_plot = numel(table_site);
mean_RMSE_pair = nanmean(table_RMSE_pair);
mean_RMSE_WMF = nanmean(table_RMSE_WMF);
mean_conn_pairp = nanmean(table_pairp);
mean_conn_WMFp = nanmean(table_WMFp);
% Write latex table 
filepath_tex = 'Z:\Documents\Research\Documents\KelbeThesis\thesis\table-fix.tex';
fid = fopen(filepath_tex, 'w+');
fprintf(fid,'\\begin{table}[H]\n');
fprintf(fid,'\\caption{Results.}\n');
fprintf(fid,'\\begin{center}{\n')
fprintf(fid,'\\tabcolsep=0.15cm\n');
fprintf(fid,'\\begin{tabular}{r r r r r r r }\n');
fprintf(fid,'\\hline\n');
fprintf(fid,' & \\multicolumn{3}{c}{Pairwise} & \\multicolumn{3}{c}{Graphical Approach} \\\\\n')
fprintf(fid,'Plot & \\gls{rmse} & \\multicolumn{2}{c}{\\%% Connected} & \\gls{rmse} & \\multicolumn{2}{c}{\\%% Connected} \\\\ \\hline\n')
for p = 1:n_plot;
    if  table_WMFp(p)<= table_pairp(p);
   fprintf(fid,' %d & %3.2f & $%d/%d=$ & $%d\\%%$ & %3.2f & $%d/%d=$ & $%d\\%%$ \\\\ \n', ...
       table_site(p), table_RMSE_pair(p), table_pair(p), ...
       table_nscan(p), table_pairp(p), ...
       table_RMSE_WMF(p), table_WMF(p), ...
       table_nscan(p), table_WMFp(p));
    else
      fprintf(fid,' %d & %3.2f & $%d/%d=$ & $%d\\%%$ & %3.2f & $%d/%d=$ & $\\mathbf{%d\\%%}$ \\\\ \n', ...
       table_site(p), table_RMSE_pair(p), table_pair(p), ...
       table_nscan(p), table_pairp(p), ...
       table_RMSE_WMF(p), table_WMF(p), ...
       table_nscan(p), table_WMFp(p));
    end
end
fprintf(fid,'\\hline\n');
fprintf(fid,' Ave. & %3.2f & & $%2.0f\\%%$ & %3.2f & & $%2.0f\\%%$ \\\\ \n', ...
        mean_RMSE_pair, mean_conn_pairp, ...
       mean_RMSE_WMF, mean_conn_WMFp);
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}}\n');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\label{tab:results}\n');
fprintf(fid,'\\end{table}\n'); 
fclose all ;

%% Collate reported RMSE error metrics 


%% Transform points into WCS

isvalid_MST = ~cellfun(@isempty,all_R_MST);
all_P_LCSw_MST = cell(n_site,1);
all_P_LCSW_MST = cell(n_site,1);
for site=1:n_site;
    n_plot = numel(all_P_LCS{site});
    if isempty(all_R_MST{site});
        continue
    end
    [~,P_ntree] = cellfun(@size,all_P_LCS{site});
    all_P_LCSw_MST{site} = cell(n_plot,1);
    all_P_LCSW_MST{site} = cell(n_plot,1);
    for j=1:n_plot;
        if isempty(all_R_MST{site}{j})
            continue
        end
        all_P_LCSw_MST{site}{j} = all_R_MST{site}{j}*all_P_LCS{site}{j} + ...
            repmat(all_t_MST{site}{j}, [1,P_ntree(j)]);
        %
        if isempty(all_R_MST{site}{13})
            continue
        end
        all_P_LCSW_MST{site}{j}= all_R_MST{site}{13}'*all_R_MST{site}{j}*...
            all_P_LCS{site}{j} + ...
            repmat(...
            all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13},...
            [1 P_ntree(j)]);
        %}
    end
end

% Find range min max for site A10  (For determining slope grade) 
%{
minc = inf;
maxc = -inf;
for i = 1:n_S;
    mint = min(all_P_LCSW_MST{8}{i}(1,:));
    maxt = max(all_P_LCSW_MST{8}{i}(1,:));
    if mint < minc;
        minc = mint;
    end
    if maxt > maxc;
        maxc = maxt;
    end
end
  %}

%% Plot results (all stem maps in WCS) 
%
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
for site = 1:n_site;
    if isempty(all_P_LCSw_MST{site})
        continue
    end
    figure;
    hold on;
    P_color = jet(n_plot);
    for j = 1:n_plot;
        n_tree = size(all_P_LCSw_MST{site}{j},2);
        for t = 1:n_tree;
            h = filledCircle(all_P_LCSw_MST{site}{j}(1:2,t),all_P_rad{site}{j}(t),500,P_color(j,:));
            set(h, 'facealpha',0.5);
        end
        if isempty(all_t_MST{site}{j});
            continue
        end
        x_axist = all_R_MST{site}{j}*x_axis + repmat(all_t_MST{site}{j},[1,2]);
        y_axist = all_R_MST{site}{j}*y_axis + repmat(all_t_MST{site}{j},[1,2]);
        z_axist = all_R_MST{site}{j}*z_axis + repmat(all_t_MST{site}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
    end
end
%}

%% Load true Stem Map
info_exp = 'Harvard';
info_suffix = 'reg';
info_slash = '\';
path_common = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, 'Common',info_slash);

all_NEON_DBH_10over = cell(n_site,1);
all_NEON_DBH_10under = cell(n_site,1);
for site = 1:n_site;
    info_site = site_unique(site);
    path_mat = sprintf('%s%03.0f%s', path_exp, info_site, '\25\mat\');
    filepath_NEON_DBH_10over_csv = sprintf('%sField_Data\\NEON_DBH_10over.csv',path_common);
    filepath_NEON_DBH_10over = sprintf('%sNEON_DBH_10over.mat',path_mat);
    if exist(filepath_NEON_DBH_10over,'file');
        load(filepath_NEON_DBH_10over);
        all_NEON_DBH_10over{site} = NEON_DBH_10over;
    elseif ~exist(filepath_NEON_DBH_10over, 'file');
        NEON_DBH_10over = parse_NEON_DBH_10over(filepath_NEON_DBH_10over_csv,info_site);
        save(filepath_NEON_DBH_10over,'NEON_DBH_10over');
        all_NEON_DBH_10over{site} = NEON_DBH_10over;
    end
    filepath_NEON_DBH_10under_csv = sprintf('%sField_Data\\NEON_DBH_10under.csv',path_common);
    filepath_NEON_DBH_10under = sprintf('%sNEON_DBH_10under.mat',path_mat);
    if exist(filepath_NEON_DBH_10under,'file');
        load(filepath_NEON_DBH_10under);
        all_NEON_DBH_10under{site} = NEON_DBH_10over;
    elseif ~exist(filepath_NEON_DBH_10under, 'file');
        NEON_DBH_10under = parse_NEON_DBH_10under(filepath_NEON_DBH_10under_csv,info_site);
        save(filepath_NEON_DBH_10under,'NEON_DBH_10under');
        all_NEON_DBH_10under{site} = NEON_DBH_10over;
    end
end


all_NEON_x = cell(n_site,1);
all_NEON_y = cell(n_site,1);
all_NEON_r = cell(n_site,1);
all_NEON_n = zeros(n_site,1);
for site = 1:n_site;
    if isempty(all_NEON_DBH_10under{site});
        continue
    end
    all_NEON_x{site} = all_NEON_DBH_10over{site}.x;
    all_NEON_y{site} = all_NEON_DBH_10over{site}.y;
    all_NEON_r{site} = all_NEON_DBH_10over{site}.dbh/(2*100);
    all_NEON_n(site) = numel(all_NEON_x{site});
end


% Plot NEON DBH maps
%{
%filepath_dbhmap = sprintf('%s%s%s',path_eps,'dbhmap','.eps');
%filepath_dbhmap= sprintf('%sdbhmap_%03.0f-%02.0f.eps',path_latex,info_site,info_plot);
%if ~exist(filepath_dbhmap,'file');
for site =1:n_site
    figure('position', [587 95 1026 832]);
    hold on
    for s = 1:all_NEON{site}
        h = filledCircle([all_NEON_x{site}(s) all_NEON_y{site}(s)],all_NEON_r{site}(s),1000,'b');
        set(h,'FaceAlpha',.5);
    end
    h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
    axis equal;
    axis([-16 16 -16 16]);
    xlabel('Easting');
    ylabel('Northing');
    grid on
    foo = 1;
    %    print(gcf,'-depsc','-opengl',filepath_dbhmap)
end
%}
%% Correct Lidar scanner coordinates via rotation matrix 
%{
all_n_plot = cellfun(@numel, all_R_MST);
for site = 1:n_site
    if isempty(all_R_MST{site});
        continue
    end
    for j = 1:n_plot
        if isempty(all_R_MST{site}{j});
            continue
        end
        all_R_MST{site}{j} = [1 0 0 ; 0 -1 0; 0 0 1]*all_R_MST{site}{j};
    end
end
%}
%%
% Plot results
%{
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
for site = 1:n_site;
    if isempty(all_P_LCSW_MST{site})
        continue
    end
    if isempty(all_R_MST{site}{13})
        continue
    end
    figure;
    hold on;
    P_color = jet(n_plot);
    P_color = P_color(randperm(n_plot),:);
    for s = 1:all_NEON_n(site)
        h = filledCircle([all_NEON_x{site}(s), -all_NEON_y{site}(s)],all_NEON_r{site}(s),1000,'k');
        set(h,'FaceAlpha',.9);
    end
    h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
    axis equal;
    axis([-16 16 -16 16]);
    xlabel('Easting');
    ylabel('Northing');    
    for j = 1:n_plot;
        n_tree = size(all_P_LCSW_MST{site}{j},2);
        for t = 1:n_tree;
            h = filledCircle(all_P_LCSW_MST{site}{j}(1:2,t),all_P_rad{site}{j}(t),500,P_color(j,:));
            set(h, 'facealpha',0.5);
        end
        if isempty(all_t_MST{site}{j});
            continue
        end
        x_axist = all_R_MST{site}{13}'*all_R_MST{site}{j}*x_axis + ...
            repmat(all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13}, [1 2]);
        y_axist = all_R_MST{site}{13}'*all_R_MST{site}{j}*y_axis + ...
            repmat(all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13}, [1 2]);
        z_axist = all_R_MST{site}{13}'*all_R_MST{site}{j}*z_axis + ...
            repmat(all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13}, [1 2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
    end
end
%}
%%
site_valid = ~cellfun(@isempty,all_t_MST);
for site = 1:n_site;
    if ~site_valid(site);
        continue
    end
    figure;
    title(sprintf('Site %d',site_unique(site))); 
    hold on; 
    axis([-15 15 -15 15]);
    xlabel('x [m]');
    ylabel('y [m]')
    for j = 1:all_n_plot(site);
        if ~isempty(all_t_pair{site}{j})&&~isempty(all_t_pair{site}{13});
        x_axist = all_R_pair{site}{13}'*all_R_pair{site}{j}*x_axis + ...
            repmat(all_R_pair{site}{13}'*all_t_pair{site}{j} - ...
            all_R_pair{site}{13}'*all_t_pair{site}{13}, [1 2]);
        y_axist = all_R_pair{site}{13}'*all_R_pair{site}{j}*y_axis + ...
            repmat(all_R_pair{site}{13}'*all_t_pair{site}{j} - ...
            all_R_pair{site}{13}'*all_t_pair{site}{13}, [1 2]);
        z_axist = all_R_pair{site}{13}'*all_R_pair{site}{j}*z_axis + ...
            repmat(all_R_pair{site}{13}'*all_t_pair{site}{j} - ...
            all_R_pair{site}{13}'*all_t_pair{site}{13}, [1 2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-m', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-m', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-m', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        end
        if ~isempty(all_t_MST{site}{j})&&~isempty(all_t_MST{site}{13});
        x_axist = all_R_MST{site}{13}'*all_R_MST{site}{j}*x_axis + ...
            repmat(all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13}, [1 2]);
        y_axist = all_R_MST{site}{13}'*all_R_MST{site}{j}*y_axis + ...
            repmat(all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13}, [1 2]);
        z_axist = all_R_MST{site}{13}'*all_R_MST{site}{j}*z_axis + ...
            repmat(all_R_MST{site}{13}'*all_t_MST{site}{j} - ...
            all_R_MST{site}{13}'*all_t_MST{site}{13}, [1 2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-r', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-r', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        end
        if ~isempty(all_t_WMF{site}{j})&&~isempty(all_t_WMF{site}{13});
        x_axist = all_R_WMF{site}{13}'*all_R_WMF{site}{j}*x_axis + ...
            repmat(all_R_WMF{site}{13}'*all_t_WMF{site}{j} - ...
            all_R_WMF{site}{13}'*all_t_WMF{site}{13}, [1 2]);
        y_axist = all_R_WMF{site}{13}'*all_R_WMF{site}{j}*y_axis + ...
            repmat(all_R_WMF{site}{13}'*all_t_WMF{site}{j} - ...
            all_R_WMF{site}{13}'*all_t_WMF{site}{13}, [1 2]);
        z_axist = all_R_WMF{site}{13}'*all_R_WMF{site}{j}*z_axis + ...
            repmat(all_R_WMF{site}{13}'*all_t_WMF{site}{j} - ...
            all_R_WMF{site}{13}'*all_t_WMF{site}{13}, [1 2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-g', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-g', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        end
        if ~isempty(all_t_SVD{site}{j})&&~isempty(all_t_SVD{site}{13});
        x_axist = all_R_SVD{site}{13}'*all_R_SVD{site}{j}*x_axis + ...
            repmat(all_R_SVD{site}{13}'*all_t_SVD{site}{j} - ...
            all_R_SVD{site}{13}'*all_t_SVD{site}{13}, [1 2]);
        y_axist = all_R_SVD{site}{13}'*all_R_SVD{site}{j}*y_axis + ...
            repmat(all_R_SVD{site}{13}'*all_t_SVD{site}{j} - ...
            all_R_SVD{site}{13}'*all_t_SVD{site}{13}, [1 2]);
        z_axist = all_R_SVD{site}{13}'*all_R_SVD{site}{j}*z_axis + ...
            repmat(all_R_SVD{site}{13}'*all_t_SVD{site}{j} - ...
            all_R_SVD{site}{13}'*all_t_SVD{site}{13}, [1 2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-b', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-b', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        end
    end
end
%% Similarity Metric


