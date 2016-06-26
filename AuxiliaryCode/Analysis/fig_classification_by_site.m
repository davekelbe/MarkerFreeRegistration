function [ ] = fig_classification_by_site( site_all, plot_all, pixel_results_all, xlabel_str )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results
fh_fa = figure;
fh_miss = figure;
fh_CR = figure;
fh_hit = figure;

n_site_all = numel(site_all);
color = jet(n_site_all);
xmin_fa = zeros(n_site_all,1);
xmax_fa = zeros(n_site_all,1);
%xmax_miss = zeros(n_site_all,1);
%xmax_CR = zeros(n_site_all,1);
%xmax_hit = zeros(n_site_all,1);
legend_text = cell(n_site_all,1);
for i = 1:n_site_all;
    xmin_fa(i) = pixel_results_all{i}.n_of_r.x(1);
    xmax_fa(i) = pixel_results_all{i}.n_of_r.x(end-1);
    legend_text{i} = sprintf('Site %02.0f',site_all(i));
    figure(fh_fa); hold on; plot(pixel_results_all{i}.n_of_r.x(1:end-1),pixel_results_all{i}.p_of_r.FA,'color',color(i,:),'linewidth',2)
    figure(fh_miss); hold on; plot(pixel_results_all{i}.n_of_r.x(1:end-1),pixel_results_all{i}.p_of_r.miss_v,'color',color(i,:),'linewidth',2)
    figure(fh_CR); hold on; plot(pixel_results_all{i}.n_of_r.x(1:end-1),pixel_results_all{i}.p_of_r.CR,'color',color(i,:),'linewidth',2)
    figure(fh_hit); hold on; plot(pixel_results_all{i}.n_of_r.x(1:end-1),pixel_results_all{i}.p_of_r.hit_v,'color',color(i,:),'linewidth',2)
end

xmin = min(xmin_fa);
xmax = max(xmax_fa);

figure(fh_fa);
title('False Alarms');
legend(legend_text);
xlim([xmin xmax ]);
ylim([0 1]);
ylabel('Percentage of laser pulses classified as False Alarm');
xlabel(xlabel_str);

figure(fh_miss);
title('Miss');
legend(legend_text);
xlim([xmin xmax]);
ylim([0 1]);
ylabel('Percentage of laser pulses classified as Miss');
xlabel(xlabel_str);

figure(fh_CR);
title('Correct Rejection');
legend(legend_text);
xlim([xmin xmax]);
ylim([0 1]);
ylabel('Percentage of laser pulses classified as Correct Rejection');
xlabel(xlabel_str);

figure(fh_hit);
title('Hit');
legend(legend_text);
xlim([xmin xmax]);
ylim([0 1]);
ylabel('Percentage of laser pulses classified as Hit');
xlabel(xlabel_str);

end

