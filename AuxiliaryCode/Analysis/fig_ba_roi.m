function [ ] = fig_ba_roi( site_all, plot_all, ba_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

figure;
n_site_all = numel(site_all);
color = jet(n_site_all);
legend_text = cell(n_site_all,1);
xmax = zeros(n_site_all,3);
hold on; 
for i = 1:n_site_all;
    xmax(i,1) = max([ba_all{i}.ROI_all ba_all{i}.SEG_all]);
    legend_text{i} = sprintf('Site %02.0f',site_all(i));
    plot(ba_all{i}.ROI_all,ba_all{i}.SEG_all,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end

for i = 1:n_site_all;
    xmax(i,2) = max([ba_all{i}.ROI_over ba_all{i}.SEG_over]);
    xmax(i,3) = max([ba_all{i}.ROI_under ba_all{i}.SEG_under]);
    plot(ba_all{i}.ROI_over,ba_all{i}.SEG_over,'o','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
    plot(ba_all{i}.ROI_under,ba_all{i}.SEG_under,'+','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end 

xmax = max(xmax(:));

title('Number of Stems ');
legend(legend_text);
xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected');
xlabel('ROI truth data');

end

