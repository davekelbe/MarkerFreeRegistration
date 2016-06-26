function [ ] = fig_dbh_seg_neon_r( site_all, plot_all, match_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

figure;
n_site_all = numel(site_all);
color = jet(n_site_all);
xmax = zeros(n_site_all,3);
hold on; 
counter = 1;
for i = 1:n_site_all;
    if ~isempty(match_all{i}.NEON_neon_xy);
    error = match_all{i}.NEON_seg_r - match_all{i}.NEON_neon_r;
    %percent_error = 100*(match_all{i}.NEON_seg_r - match_all{i}.NEON_neon_r)./match_all{i}.NEON_seg_r;
    %xmax(i,1) = max(max([match_all{i}.NEON_neon_r match_all{i}.NEON_seg_r]));
    legend_text{counter} = sprintf('Site %02.0f',site_all(i));
    plot(match_all{i}.NEON_neon_r,error,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
    counter = counter + 1;
    end
end

xmax = 2*max(xmax(:));

title('DBH Error');
hleg = legend(legend_text);
set(hleg,'location','bestoutside');
set(hleg,'edgecolor','w');
%xlim([0 xmax ]);
%ylim([0 xmax]);
ylabel('Error');
xlabel('DBH');

end

