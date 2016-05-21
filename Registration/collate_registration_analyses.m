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

info_exp = 'Harvard';
info_suffix = 'reg';
info_slash = '\';

site_unique = unique(C{1});
n_site = numel(site_unique);
all_site = [];
all_plotI = [];
all_plotJ = [];
all_noise_RMSE = [];
all_noise_step = [];
all_pose_match_RMSE = [];
all_pose_true_rx = [];
all_pose_true_ry = [];
all_pose_true_rz = [];
all_pose_true_tx = [];
all_pose_true_ty = [];
all_pose_true_tz = [];
all_sub_RMSE = [];
all_sub_step = [];
all_noise_RMSEin = [];
all_noise_RMSEout = [];
all_dbh_step = [];
all_dbh_RMSE = [];


n_plot = zeros(n_site,1);

exclude_site = [];
for s = 1:11;%n_site;
    info_site = site_unique(s);
    if ismember(exclude_site, info_site);
        continue
    end
    path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
    
    for p = 1:25;
    path_mat = sprintf('%s%02.0f%smat%s',path_site,p,info_slash,info_slash);
    filepath_noise_RMSE= sprintf('%s%s',path_mat, 'anal_noise_RMSE.mat');
    if ~exist(filepath_noise_RMSE, 'file');
        continue
    end
    load(filepath_noise_RMSE);
    anal_noise_RMSE = match_RMSE;
    filepath_noise_step = sprintf('%s%s',path_mat, 'anal_noise_step.mat');
    load(filepath_noise_step);
    anal_noise_step = step_std;
    filepath_pose_match_RMSE = sprintf('%s%s',path_mat, 'anal_pose_match_RMSE.mat');
    load(filepath_pose_match_RMSE);
    anal_pose_match_RMSE = match_RMSE;
    filepath_pose_true_rx = sprintf('%s%s',path_mat, 'anal_pose_true_rx.mat');
    load(filepath_pose_true_rx);
    anal_pose_true_rx = rx;
    filepath_pose_true_ry = sprintf('%s%s',path_mat, 'anal_pose_true_ry.mat');
    load(filepath_pose_true_ry);
    anal_pose_true_ry = ry;
    filepath_pose_true_rz = sprintf('%s%s',path_mat, 'anal_pose_true_rz.mat');
    load(filepath_pose_true_rz);
    anal_pose_true_rz = rz;
    filepath_pose_true_tx = sprintf('%s%s',path_mat, 'anal_pose_true_tx.mat');
    load(filepath_pose_true_tx);
    anal_pose_true_tx = tx;
    filepath_pose_true_ty = sprintf('%s%s',path_mat, 'anal_pose_true_ty.mat');
    load(filepath_pose_true_ty);
    anal_pose_true_ty = ty;
    filepath_pose_true_tz = sprintf('%s%s',path_mat, 'anal_pose_true_tz.mat');
    load(filepath_pose_true_tz);
    anal_pose_true_tz = tz;
    filepath_sub_RMSE = sprintf('%s%s',path_mat, 'anal_sub_RMSE.mat');
    load(filepath_sub_RMSE);
    anal_sub_RMSE = match_RMSE;
    filepath_sub_step = sprintf('%s%s',path_mat, 'anal_sub_step.mat');
    load(filepath_sub_step);
    anal_sub_step = step_sub;
    filepath_noise_RMSEin = sprintf('%s%s',path_mat, 'anal_noise_RMSEin.mat');
    load(filepath_noise_RMSEin);
    anal_noise_RMSEin = match_RMSEin;
    filepath_noise_RMSEout = sprintf('%s%s',path_mat, 'anal_noise_RMSEout.mat');
    load(filepath_noise_RMSEout);
    anal_noise_RMSEout = match_RMSEout;
    filepath_dbh_step = sprintf('%s%s',path_mat, 'anal_dbh_step.mat');
    load(filepath_dbh_step);
    anal_dbh_step = step_std;
    filepath_dbh_RMSE = sprintf('%s%s',path_mat, 'anal_dbh_RMSE.mat');
    load(filepath_dbh_RMSE);
    anal_dbh_RMSE = match_RMSE;
    
    %all_site = [all_site; repmat(info_site, [n_plot(s).^2,1])];
   % all_plotI = [all_plotI; match_plotI(:)];
   % all_plotJ = [all_plotJ; match_plotJ(:)];
    all_noise_RMSE = [all_noise_RMSE; anal_noise_RMSE(:)];
    all_noise_step = [all_noise_step; anal_noise_step(:)];
    all_pose_match_RMSE = [all_pose_match_RMSE; anal_pose_match_RMSE(:)];
    all_pose_true_rx = [all_pose_true_rx; anal_pose_true_rx(:)];
    all_pose_true_ry = [all_pose_true_ry; anal_pose_true_ry(:)];
    all_pose_true_rz = [all_pose_true_rz; anal_pose_true_rz(:)];
    all_pose_true_tx = [all_pose_true_tx; anal_pose_true_tx(:)];
    all_pose_true_ty = [all_pose_true_ty; anal_pose_true_ty(:)];
    all_pose_true_tz = [all_pose_true_tz; anal_pose_true_tz(:)];
    all_sub_RMSE = [all_sub_RMSE; anal_sub_RMSE(:)];
    all_sub_step = [all_sub_step; anal_sub_step(:)];
    all_noise_RMSEin = [all_noise_RMSEin; anal_noise_RMSEin(:)];
    all_noise_RMSEout = [all_noise_RMSEout; anal_noise_RMSEout(:)];
    all_dbh_step = [all_dbh_step; anal_dbh_step(:)];
    all_dbh_RMSE = [all_dbh_RMSE; anal_dbh_RMSE(:)];
    end
    
    %{
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    [x_sph, y_sph, z_sph] = sphere;
    figure
    hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    % Plot all pairwise matches from j
    i = 1;
    for j = 1:n_plot(s);
        % if ~isempty0(i,j);
        x_axist = match12_R{i,j}*x_axis + repmat(match12_t{i,j},[1,2]);
        y_axist = match12_R{i,j}*y_axis + repmat(match12_t{i,j},[1,2]);
        z_axist = match12_R{i,j}*z_axis + repmat(match12_t{i,j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),.1+x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),.1+y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),.1+z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
        % h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
        %     z_sph+match_t{i,j}(3));
        % alpha(0.2)
        % set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        %   end
    end
    %}
    
end

%RMSE vs noise
%{
figure;
xrange = sort(unique(all_noise_step));
RMSE_mean = nan(size(xrange));
for i = 1:numel(xrange);
    RMSE_mean(i) = nanmean(all_noise_RMSE(all_noise_step==xrange(i)));
end
boxplot(all_noise_RMSE, all_noise_step,'widths',[.9 .9 .9 .9 .9]);
hold on
ylim([0 1]);
lw = 1.5;
plot(1:numel(xrange), RMSE_mean, '*k');
plot([1 numel(xrange)],[0 .5], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Noise (\sigma) added to Tie Points  ');
ylabel('RMSE [m]')
set(gca,'xtick', 1:2:numel(xrange));
set(gca,'xticklabel', xrange(1:2:numel(xrange)));
% Regression 
is_valid = ~isnan(RMSE_mean);
p = polyfit(xrange(is_valid), RMSE_mean(is_valid),2);
xint = 1:numel(xrange);
yfit = polyval(p,xrange);
plot(xint, yfit, '--', 'color', [0.5 0.5 0.5], 'linewidth', lw)
clear legend_str;
legend_str{1} = 'Mean';
legend_str{2} = '1:1';
legend_str{3} = sprintf('Fit: RMSE = %3.2f*sigma^2 + %3.2fsigma + %3.2f',...
    p(1), p(2), p(3));
legend(legend_str, 'location', 'best');
filepath_rmsevsnoise = 'Z:\Desktop\rmsevsnoise.tex';
matlab2tikz(filepath_rmsevsnoise)
%}

%{
all_pose_true_rxdeg = rad2deg(all_pose_true_rx);
all_pose_true_rydeg = rad2deg(all_pose_true_ry);
all_pose_true_rzdeg = rad2deg(all_pose_true_rz);
all_pose_true_rmeandeg = mean([all_pose_true_rxdeg all_pose_true_rydeg all_pose_true_rzdeg],2);
figure; 
int = 30;
xrange = 0:int:max([all_pose_true_rxdeg;all_pose_true_rydeg;all_pose_true_rzdeg]);
all_ownership = zeros(size(all_pose_true_rxdeg));
label = cell(numel(xrange),1);
RMSE_mean = zeros(1,numel(xrange)-1);
for i = 1:numel(xrange);
    ix_int = (all_pose_true_rmeandeg>= xrange(i)-int/2 & all_pose_true_rmeandeg< xrange(i)+int/2);
    all_ownership(ix_int) = xrange(i);
    RMSE_mean(i) = nanmean(all_pose_match_RMSE(ix_int));
end
boxplot(all_pose_match_RMSE, all_ownership,'widths',[.9 .9 .9 .9 .9]);
set(gca, 'ylim', [-1e-15 1e-13])
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Mean Rotation Angle [deg]');
ylabel('RMSE [m]')
filepath_rmsevspose = 'Z:\Desktop\rmsevspose.tex';
matlab2tikz(filepath_rmsevspose)
%}

%{
figure
lw = 1.5;
xrange = sort(unique(all_sub_step));
RMSE_mean = nan(size(xrange));
for i = 1:numel(xrange);
    RMSE_mean(i) = nanmean(all_sub_RMSE(all_sub_step==xrange(i)));
end
boxplot(all_sub_RMSE, all_sub_step,'widths',[.9 .9 .9 .9 .9]);
hold on
plot(1:numel(xrange), RMSE_mean, '*k');
plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Number of Subset Tie Points');
ylabel('RMSE [m]')
% Regression on mean 
%{
dummy = [1./xrange ones(size(xrange))];
b = regress(RMSE_mean, dummy);
xrange2 = (min(xrange):1:max(xrange))';
dummy2 = [1./xrange2 ones(size(xrange2))];
yfit2 = dummy2*b;
xint = 1:numel(xrange);
xint2 = linspace(min(xint), max(xint), numel(yfit2));
plot(xint2,yfit2,'--', 'color', [.5 .5 .5], 'linewidth', lw)
%}
dummy = [1./all_sub_step ones(size(all_sub_step))];
b = regress(all_sub_RMSE, dummy);
xrange2 = (min(all_sub_step):1:max(all_sub_step))';
dummy2 = [1./xrange2 ones(size(xrange2))];
yfit2 = dummy2*b;
xint = 1:numel(xrange);
xint2 = linspace(min(xint), max(xint), numel(yfit2));
plot(xint2,yfit2,'--', 'color', [.5 .5 .5], 'linewidth', lw)
clear legend_str
legend_str{1} = 'Mean';
legend_str{2} = 'Added RMSE noise floor';
legend_str{3} = sprintf('Fit: RMSE = %3.2f/N + %3.2f', b(1), b(2));
legend(legend_str);
filepath_rmsevssub = 'Z:\Desktop\rmsevssub.tex';
matlab2tikz(filepath_rmsevssub)
foo = 1;
%}

%RMSEout vs RMSEin
%{
figure;
is_valid = ~isnan(all_noise_RMSEin)&~isnan(all_noise_RMSEout);
axis([0 1 0 2.5]);
hold on
lw = 2.5;
plot([0 1.5], [0 1.5], ':', 'color', [.5 .5 .5], 'linewidth', lw);
plot(all_noise_RMSEin, all_noise_RMSEout, '.k', 'markersize', 5)
%p = polyfit(all_noise_RMSEin(is_valid), all_noise_RMSEout(is_valid),2);
%xint = 0:.1:1.5;
%yfit = polyval(p,xint);
%plot(xint, yfit, '--r');
p = polyfit(all_noise_RMSEin(is_valid), all_noise_RMSEout(is_valid),1);
xint = 0:.1:1.5;
yfit = polyval(p,xint);
%plot(xint, yfit, '--','color', [.85 .85 .85], 'linewidth', lw);
xlabel('RMSE in [m]');
ylabel('RMSE out [m]')
clear legend_str;
legend_str{1} = '1:1';
%legend_str{2} = sprintf('Fit: RMSE = %3.2fsigma + %3.2f',...
%    p(1), p(2));
legend(legend_str, 'location', 'best');
filepath_rmsevsrmse = 'Z:\Desktop\rmsevsrmse.tex';
matlab2tikz(filepath_rmsevsrmse)
%}

% DBH vs RMSE
figure;
scatter(2*all_dbh_step, all_dbh_RMSE);
foo = 1;
all_dbh_step2 = all_dbh_step*100;
figure;
boxplot(all_dbh_RMSE,2*all_dbh_step2,'widths',[.9 .9 .9 .9 .9]);
hold on
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([1 numel(xrange)],[0 .5], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Percent error (\sigma) added to Diameters');
ylabel('RMSE [m]')
filepath_rmsevsdbh = 'Z:\Desktop\rmsevsdbh.tex';
matlab2tikz(filepath_rmsevsdbh)




