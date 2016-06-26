function [ ] = fig_match_seg_roi( site_all, plot_all, match_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

figure;
n_site_all = numel(site_all);
color = jet(n_site_all);
hold on;

x = [];
y = [];
counter = 1;
for i = 1:n_site_all;
    x = [x ;2*match_all{i}.ROI_roi_r];
    y = [y ;2*match_all{i}.ROI_seg_r];
    if ~isempty(match_all{i}.ROI_roi_r);    
        legend_text{counter} = sprintf('Site %02.0f',site_all(i));
        plot(2*match_all{i}.ROI_roi_r,2*match_all{i}.ROI_seg_r,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
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
plot([0 xmax],[0 xmax],'-k');
plot([0 xmax],[polyval(p,0) polyval(p,xmax)],':k')


title('DBH Retrieval ');
hleg = legend(legend_text);
set(hleg,'location','bestoutside');
set(hleg,'edgecolor','w');
xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected [m]');
xlabel('ROI truth data [m]');
htxt = annotation('textbox',[.2 .7 .1 .1], 'string', sprintf('R^2 = %4.2f',rsq));
set(htxt,'edgecolor','w');

end

