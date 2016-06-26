function [ ] = fig_n_stems_roi( site_all, plot_all, n_stems_all )
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
    x = [x ; n_stems_all{i}.ROI_all];
    y = [y ; n_stems_all{i}.SEG_all];
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
 %   xmax(i,1) = max([n_stems_all{i}.ROI_all n_stems_all{i}.SEG_all]);
    legend_text{i} = sprintf('Site %02.0f',site_all(i));
    plot(n_stems_all{i}.ROI_all,n_stems_all{i}.SEG_all,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end

%xmax = max(xmax(:));
legend_text{i+1} = sprintf('Fit');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
legend_text{i+2} = sprintf('1:1');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
plot(xq,yq,'-k');
plot(xq,xq,':k');

%{
for i = 1:n_site_all;
    xmax(i,2) = max([n_stems_all{i}.ROI n_stems_all{i}.SEG]);
   % xmax(i,3) = max([n_stems_all{i}.ROI n_stems_all{i}.SEG]);
    plot(n_stems_all{i}.ROI_all, n_stems_all{i}.SEG_all,'o','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
    %plot(n_stems_all{i}.ROI_unde,n_stems_all{i}.SEG_under,'+','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end 
%}


title('Number of Stems ');
hleg = legend(legend_text);
set(hleg, 'edgecolor','w');
%set(hleg, 'position', [.85 .4 .2 .2]);
set(hleg, 'location', 'bestoutside');
xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected');
xlabel('ROI truth data');
htxt = annotation('textbox',[.5 .8 .1 .1], 'string', sprintf('R^2 = %4.2f',rsq));
set(htxt,'edgecolor','w');



end

