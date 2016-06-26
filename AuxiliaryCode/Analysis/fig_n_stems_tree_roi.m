function [ ] = fig_n_stems_tree_roi( site_all, plot_all, n_stems_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results
%{
figure;
n_site_all = numel(site_all);
color = jet(n_site_all);
legend_text = cell(n_site_all,1);
xmax = zeros(n_site_all,3);
hold on; 
%}
n_site_all = numel(site_all);
x = [];
y = [];

for i = 1:n_site_all
    x = [x ; n_stems_all{i}.ROI_all];
    y = [y ; n_stems_all{i}.TREE_all];
end
x = x*25; y = y*25;
% OLS
p = polyfit(x,y,1);
yfit  = polyval(p,x);
fiteq =  @(x) p(1)*x + p(2);
b1 = p(1); 
b0 = p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
xmax = max([x;y]);
xq = [0 xmax];
yq = polyval(p,xq);

n = numel(x);
s = std(yresid);
xbar = mean(x);
Sxx = sum((x -xbar).^2);
dof = n - 2;
sb0 = s*sqrt((1/n) + (xbar.^2/Sxx));
sb1 = s/sqrt(Sxx);
rmse = sqrt(mean(yresid.^2));

test_param(b1, sb1, 1, .05, dof);
test_param(b1, sb1, 1, .01, dof);
test_param(b0, sb0, 0, .05, dof);
test_param(b0, sb0, 0, .01, dof);

alpha = 0.05; 
t_a = tinv(1-(alpha/2),dof);
x_step = linspace(0, 25*max(x));
MSresid = SSresid/(n - 2);
sy = sqrt(MSresid*(1 + (1/n) + (x_step - mean(x_step)).^2/Sxx));
PI_delta = t_a*sy;
PI_u = x_step + PI_delta;
PI_l = x_step - PI_delta;

%for i = 1:n_site_all;
% %   xmax(i,1) = max([n_stems_all{i}.ROI_all n_stems_all{i}.SEG_all]);
%    legend_text{i} = sprintf('Site %02.0f',site_all(i));
%    plot(n_stems_all{i}.ROI_all,n_stems_all{i}.TREE_all,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
%end

maxx = max(x);
figure
plot(25*x,25*y,'+k', 'linewidth',2, 'markersize', 10);
grid off;
axis auto
axis square
hold on 
plot(25*[0 maxx], 25*[fiteq(0) fiteq(maxx)],'-', 'linewidth',2, 'color', [.7 .7 .7] );
plot(25*[0 maxx], 25*[0 maxx],'--k', 'linewidth', 2);
%plot(x_step, PI_u,':', 'linewidth', 1.5, 'color', [.7 .7 .7]);
plot(4,2,'-w');
plot(4,2,'-w');
%plot(x_step, PI_l,':', 'linewidth', 1.5, 'color', [.7 .7 .7]);
legend_str{1} = 'Data';
legend_str{2} = sprintf('y = %2.2fx + %2.2f', b1 , b0);
legend_str{3} = '1:1';
%legend_str{4} = '95% P.I.';
legend_str{4} = sprintf('R^2 = %3.2f',rsq);
legend_str{5} = sprintf('RMSE = %3.2f stems ha^{-1}',25*rmse);
%annotation('textbox', [.85, .2,.2,.2], 'string', sprintf('R^2 = %3.3f',rsq));
%text(14,4,['\fontsize{16}R^2'])
hleg = legend(legend_str);
set(hleg, 'position', [.85 .4 .2 .2]);
legend boxoff
axis([0 25*maxx 0 25*maxx]);
ticky = get(gca, 'ytick');
set(gca, 'xtick', ticky);
set(gca,'box','off') 
xlabel('Visible (ROI) stems ha^{-1}');
ylabel('Estimated (lidar) stems ha^{-1}');
path_latex = 'Z:\Desktop\Paper\';
filename = sprintf('%s%s', path_latex,'n_stems_tree_roi.tex');
matlab2tikz(filename);




%{
plot(x,y,'+k','markersize',10,'linewidth',2);
i = 1;
legend_text{i} = 'Data';
%xmax = max(xmax(:));
legend_text{i+1} = sprintf('Fit');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
legend_text{i+2} = sprintf('1:1');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
plot(xq,yq,'-', 'linewidth', 2, 'color', [.7 .7 .7]);
plot(xq,xq,'--k', 'linewidth', 2);

%{
for i = 1:n_site_all;
    xmax(i,2) = max([n_stems_all{i}.ROI n_stems_all{i}.SEG]);
   % xmax(i,3) = max([n_stems_all{i}.ROI n_stems_all{i}.SEG]);
    plot(n_stems_all{i}.ROI_all, n_stems_all{i}.SEG_all,'o','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
    %plot(n_stems_all{i}.ROI_unde,n_stems_all{i}.SEG_under,'+','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end 
%}


hleg = legend(legend_text);
set(hleg, 'edgecolor','w');
set(hleg, 'position', [.85 .4 .2 .2]);
%set(hleg, 'location', 'bestoutside');

xlim([0 xmax ]);
ylim([0 xmax]);
%htxt = annotation('textbox',[.5 .8 .1 .1], 'string', sprintf('R^2 = %4.2f',rsq));
%set(htxt,'edgecolor','w');
legend boxoff
set(gca,'box','off') 
xlabel('Measured (truth) Number of trees [cm]');
ylabel('Estimated (lidar) tree DBH [cm]');

clear legend_str
[ xfit, yfit, m, b, r2_orthog, stde ] = orthog_fit( x,y );
fiteq =  @(x) m*x + b;
maxx = max(x);
figure
plot(25*x,25*y,'+k', 'linewidth',2, 'markersize', 10);
grid off;
axis auto
axis square
hold on 
plot(25*[0 maxx], 25*[fiteq(0) fiteq(maxx)],'-', 'linewidth',2, 'color', [.7 .7 .7] );
plot(25*[0 maxx], 25*[0 maxx],'--k', 'linewidth', 2);
plot(.5,.01,'-w');
plot(.5,.01,'-w');
legend_str{1} = 'Data';
legend_str{2} = 'Fit';
legend_str{3} = '1:1';
legend_str{4} = sprintf('R^2 = %3.2f',r2_orthog);
legend_str{5} = sprintf('\\sigma_e = %3.2f stems ha^{-1}',25*stde);
%annotation('textbox', [.85, .2,.2,.2], 'string', sprintf('R^2 = %3.3f',rsq));
%text(14,4,['\fontsize{16}R^2'])
hleg = legend(legend_str);
set(hleg, 'position', [.85 .4 .2 .2]);
legend boxoff
set(gca,'box','off') 
xlabel('Visible (ROI) stems ha^{-1}');
ylabel('Estimated (lidar) stems ha^{-1}');
path_latex = 'Z:\Desktop\Paper\';
filename = sprintf('%s%s', path_latex,'n_stems_tree_roi.tex');
matlab2tikz(filename);

%}
end

