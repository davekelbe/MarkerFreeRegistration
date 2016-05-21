function [ ] = fig_ba_seg_roi( site_all, plot_all, ba_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

figure;
n_site_all = numel(site_all);
color = jet(n_site_all);
legend_text = cell(n_site_all,1);
xmax = zeros(n_site_all,3);
hold on; 

x = [];
y = [];
for i = 1:n_site_all
    x = [x ; ba_all{i}.ROI_all];
    y = [y ; ba_all{i}.SEG_all];
end
p = polyfit(x,y,1);
yfit  = polyval(p,x);
%yfit = p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
xmax = max([x;y]);
xq = [0 xmax];
yq = polyval(p,xq);

for i = 1:n_site_all;
    legend_text{i} = sprintf('Site %02.0f',site_all(i));
    plot(ba_all{i}.ROI_all,ba_all{i}.SEG_all,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end

%{
for i = 1:n_site_all;
    xmax(i,2) = max([ba_all{i}.ROI_over ba_all{i}.SEG_over]);
    xmax(i,3) = max([ba_all{i}.ROI_under ba_all{i}.SEG_under]);
    plot(ba_all{i}.ROI_over,ba_all{i}.SEG_over,'o','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
    plot(ba_all{i}.ROI_under,ba_all{i}.SEG_under,'+','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end 
%}
legend_text{i+1} = sprintf('Fit');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
legend_text{i+2} = sprintf('1:1');

plot(xq,yq,'-k');
plot(xq,xq,':k');

title('Basal Area ');
hleg = legend(legend_text);
set(hleg, 'edgecolor','w');
set(hleg, 'location','bestoutside');
xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected');
xlabel('ROI truth data');
htxt = annotation('textbox',[.5 .8 .1 .1], 'string', sprintf('R^2 = %4.2f',rsq));
set(htxt,'edgecolor','w');

end

