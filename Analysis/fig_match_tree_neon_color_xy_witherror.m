function [ ] = fig_match_tree_neon_color_xy_witherror( site_all, plot_all, match_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

figure;



n_site_all = numel(site_all);
color = jet(n_site_all);
%legend_text = cell(n_site_all,1);
hold on; 

x = [];
y = [];
xy = [];
for i = 1:n_site_all;
    x = [x ;2*match_all{i}.NEON_neon_r];
    y = [y ;2*match_all{i}.NEON_tree_r];
    xy = [xy; match_all{i}.NEON_tree_xy];
end

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
rmse = sqrt(mean(yresid.^2));

n = numel(x);
s = std(yresid);
xbar = mean(x);
Sxx = sum((x -xbar).^2);
dof = n - 2;
sb0 = s*sqrt((1/n) + (xbar.^2/Sxx));
sb1 = s/sqrt(Sxx);

test_param(b1, sb1, 1, .05, dof);
test_param(b1, sb1, 1, .01, dof);
test_param(b0, sb0, 0, .05, dof);
test_param(b0, sb0, 0, .01, dof);

alpha = 0.05; 
t_a = tinv(1-(alpha/2),dof);
x_step = linspace(0, max(x));
MSresid = SSresid/(n - 2);
sy = sqrt(MSresid*(1 + (1/n) + (x_step - mean(x_step)).^2/Sxx));
PI_delta = t_a*sy;
PI_u = x_step + PI_delta;
PI_l = x_step - PI_delta;
 
% TLS 
%[xfit_o, yfit_o, m, b, rsq, stde] = orthog_fit(x, y);
%fiteq =  @(x) m*x + b;

[sampling, beam_dia,~,~] = SICK_specs(xy,1);
sampling_dist = sampling*xy;
maxx = max(x);
figure
errorbar(100*x,100*y,beam_dia,'.k');
plot(100*x,100*y,'+k', 'linewidth',2, 'markersize', 10);
grid off;
axis auto
axis square
hold on 
plot([0 100*maxx], 100*[fiteq(0) fiteq(maxx)],'-', 'linewidth',2, 'color', [.7 .7 .7] );
plot([0 100*maxx], 100*[0 maxx],'--k', 'linewidth', 2);
%plot(100*x_step, 100*PI_u,':', 'linewidth', 1.5, 'color', [.7 .7 .7]);
plot(4,2,'-w');
plot(4,2,'-w');
%plot(100*x_step, 100*PI_l,':', 'linewidth', 1.5, 'color', [.7 .7 .7]);
legend_str{1} = 'Data';
legend_str{2} = sprintf('y = %2.2fx + %2.2f', b1 , b0);
legend_str{3} = '1:1';
%legend_str{4} = '95% P.I. ';
legend_str{4} = sprintf('R^2 = %2.2f',rsq);
legend_str{5} = sprintf('RMSE = %3.2f cm',100*rmse);
hleg = legend(legend_str);
set(hleg, 'position', [.75 .4 .2 .2]);
legend boxoff
set(gca,'box','off') 
axis(100*[0 maxx 0 maxx]);
xlabel('Measured (truth) tree DBH [cm]');
ylabel('Estimated (lidar) tree DBH [cm]');
path_latex = 'Z:\Desktop\Paper\';
filename = sprintf('%s%s', path_latex,'dbh_vs_xy.tex');
matlab2tikz(filename);



%{
xy_step = 2;
xy_array = 0:xy_step:max(xy);
n_xy = numel(xy_array)-1;
color = jet(n_xy);
counter = 1;
for i = 1:n_xy;
    is_valid = xy>xy_array(i)&xy<=xy_array(i+1);
   if sum(is_valid)>0;
        legend_text{counter} = sprintf('%2.0f - %2.0f m',xy_array(i),xy_array(i+1));
        plot(x(is_valid),y(is_valid),'x',...
            'markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
        counter = counter + 1;
   end
end

xmax = max([x;y]);
p = polyfit(x,y,1);
yfit  = polyval(p,x);
%yfit = p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;

legend_text{counter} = sprintf('Fit');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
legend_text{counter+1} = sprintf('1:1');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);

hold on;
plot([0 xmax],[polyval(p,0) polyval(p,xmax)],'-k')
plot([0 xmax],[0 xmax],':k');
htxt = annotation('textbox',[.2 .8 .1 .1], 'string', sprintf('R^2 = %4.2f',rsq));
set(htxt,'edgecolor','w');


title('DBH Retrieval ');
hleg = legend(legend_text);
set(hleg, 'edgecolor','w');
set(hleg, 'location', 'bestoutside');
hlegtitle = get(hleg,'title');
set(hlegtitle,'string', 'xy-Range');
xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected [m]');
xlabel('NEON truth data [m]');
%}
end

