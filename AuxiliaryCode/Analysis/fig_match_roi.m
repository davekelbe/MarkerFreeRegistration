function [ ] = fig_match_neon( site_all, plot_all, match_all )
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
    xmax(i,1) = max(max([match_all{i}.NEON_neon_r match_all{i}.NEON_seg_r]));
    legend_text{i} = sprintf('Site %02.0f',site_all(i));
    plot(2*match_all{i}.NEON_neon_r,2*match_all{i}.NEON_seg_r,'x','markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
end

xmax = 2*max(xmax(:));

title('DBH Retrieval ');
legend(legend_text);
xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected [m]');
xlabel('NEON truth data [m]');

end

