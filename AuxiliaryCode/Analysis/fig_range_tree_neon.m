function [ ] = fig_range_tree_neon( site_all, plot_all, match_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

n_site_all = numel(site_all);
color = jet(n_site_all);
%legend_text = cell(n_site_all,1);
hold on; 

x = [];
y = [];
for i = 1:n_site_all;
    x = [x ;match_all{i}.NEON_neon_xy];
    y = [y ;match_all{i}.NEON_tree_xy];
end
%{
for i = 1:n_site_all;
counter = 1;
   % xmax(i,1) = max(max([match_all{i}.NEON_neon_r match_all{i}.NEON_seg_r]));
   if ~isempty(match_all{i}.NEON_neon_r);
        legend_text{counter} = sprintf('Site %02.0f',site_all(i));
        plot(match_all{i}.NEON_neon_xy,match_all{i}.NEON_tree_xy,'x',...
            'markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
        counter = counter + 1;
   end 
end
%}

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

maxx = max(x);
figure
plot(x,y,'+k', 'linewidth',2, 'markersize', 10);
grid off;
axis auto
axis square
hold on 
plot([0 maxx], [fiteq(0) fiteq(maxx)],'-', 'linewidth',2, 'color', [.7 .7 .7] );
plot([0 maxx], [0 maxx],'--k', 'linewidth', 2);
%plot(x_step, PI_u,':', 'linewidth', 1.5, 'color', [.7 .7 .7]);
plot(4,2,'-w');
plot(4,2,'-w');
%plot(x_step, PI_l,':', 'linewidth', 1.5, 'color', [.7 .7 .7]);
legend_str{1} = 'Data';
legend_str{2} = sprintf('y = %2.2fx + %2.2f', b1 , b0);
legend_str{3} = '1:1';
%legend_str{4} = '95% P.I. ';
legend_str{4} = sprintf('R^2 = %3.3f',rsq);
legend_str{5} = sprintf('RMSE = %3.2f cm',100*rmse);
%annotation('textbox', [.85, .2,.2,.2], 'string', sprintf('R^2 = %3.3f',rsq));
%text(14,4,['\fontsize{16}R^2'])
hleg = legend(legend_str);
set(hleg, 'position', [.75 .4 .2 .2]);
legend boxoff
set(gca,'box','off') 
axis([0 maxx 0 maxx]);
xlabel('Measured (truth) tree range [m]');
ylabel('Estimated (lidar) tree range [m]');
path_latex = 'Z:\Desktop\Paper\';
filename = sprintf('%s%s', path_latex,'range_tree_neon.tex');
matlab2tikz(filename);

end

