function [ ] = fig_classification_by_site_all( site_all, plot_all, pixel_results_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results
figure;
hold on
n_site_all = numel(site_all);

for i = 1:n_site_all;
    xtick_str{i} = sprintf('%2.0f',site_all(i));
    plot(i,100*pixel_results_all{i}.p.CR,'xk','linewidth',2);
    plot(i,100*pixel_results_all{i}.p.hit_v,'xg','linewidth',2);
    plot(i,100*pixel_results_all{i}.p.miss_v,'xb','linewidth',2);
    plot(i,100*pixel_results_all{i}.p.FA,'xr','linewidth',2);
end
hleg = legend('Correct rejection','Correct detection','Incorrect rejection','Incorrect detection');
set(hleg, 'position',[1 .8 .4 .4])
legend boxoff 
set(gca, 'xtick',1:n_site_all);
set(gca, 'ytick',0:10:100);
set(gca, 'xticklabel',xtick_str);
%title('Classification Accuracy');
ylim([0 100]);
ylabel('Percentage');
xlabel('Plot Number');
grid on

path_latex = 'Z:\Desktop\Paper\';
filename = sprintf('%s%s', path_latex,'classification_by_site_all.tex');
matlab2tikz(filename);

end

