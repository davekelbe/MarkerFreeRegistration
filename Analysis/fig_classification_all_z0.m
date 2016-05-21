function [ ] = fig_classification_all_z0( site_all, plot_all, pixel_results_all, xlabel_str )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

n_site_all = numel(site_all);
xmin_fa = zeros(n_site_all,1);
xmax_fa = zeros(n_site_all,1);
%xmax_miss = zeros(n_site_all,1);
%xmax_CR = zeros(n_site_all,1);
%xmax_hit = zeros(n_site_all,1);
legend_text = cell(n_site_all,1);


for i = 1:n_site_all;
    xmin_fa(i) = pixel_results_all{i}.n_of_z0.x(1);
    xmax_fa(i) = pixel_results_all{i}.n_of_z0.x(end-1); 
end

[maxval,I] = max(xmax_fa);
x_arr = pixel_results_all{I}.n_of_z0.x;
n_x = numel(x_arr);

x = zeros(n_site_all,n_x);
n_of = zeros(n_site_all,n_x);
n_v = zeros(n_site_all,n_x);
n_o = zeros(n_site_all,n_x);
n_vo = zeros(n_site_all,n_x);
n_n = zeros(n_site_all,n_x);

for i = 1:n_site_all;
    is_valid = false(1,n_x);
    is_valid(1:numel(pixel_results_all{i}.n_of_z0.x)) = true;
    x(i,is_valid) = pixel_results_all{i}.n_of_z0.x;
    n_of_FA(i,is_valid) = pixel_results_all{i}.n_of_z0.FA;
    n_of_miss(i,is_valid) = pixel_results_all{i}.n_of_z0.miss_v;
    n_of_hit(i,is_valid) = pixel_results_all{i}.n_of_z0.hit_v;
    n_of_CR(i,is_valid) = pixel_results_all{i}.n_of_z0.CR;

    n_v(i,is_valid) = pixel_results_all{i}.n_of_z0.v;
    n_o(i,is_valid) = pixel_results_all{i}.n_of_z0.o;
    n_vo(i,is_valid) = pixel_results_all{i}.n_of_z0.vo;
    n_n(i,is_valid) = pixel_results_all{i}.n_of_z0.n;
end

n_of_FA = sum(n_of_FA,1);
n_of_miss = sum(n_of_miss,1);
n_of_hit = sum(n_of_hit,1);
n_of_CR = sum(n_of_CR,1);
n_v_all = sum(n_v,1);
n_o_all = sum(n_o,1);
n_vo_all = sum(n_vo,1);
n_n_all = sum(n_n,1);

p_of_FA = n_of_FA./n_n_all;
p_of_miss = n_of_miss./n_v_all;
p_of_hit = n_of_hit./n_v_all;
p_of_CR = n_of_CR./n_n_all;


xmin = min(xmin_fa);
xmax = max(xmax_fa);

fh = figure;
hold on;
plot(x_arr,p_of_FA,'r','linewidth',2);
plot(x_arr,p_of_hit,'g','linewidth',2);
plot(x_arr,p_of_miss,'b','linewidth',2);
plot(x_arr,p_of_CR,'k','linewidth',2);
title('Classification Accuracy');% - Aggregate over all sites');
hleg = legend('FA','Hit','Miss','CR');
set(hleg, 'location','bestoutside');
xlim([xmin xmax ]);
ylim([0 1]);
ylabel('Percentage of laser pulses');
xlabel(xlabel_str);

end

