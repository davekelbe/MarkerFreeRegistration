function [ p_of_r ] = Pixel_Classification_seg( I_seg_id,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,seg_dia,independent_var )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r_step = 0.05;
r_array = 0:r_step:max(seg_dia)+r_step;
n_r = numel(r_array)-1;
p_of_r.FA = zeros(n_r,1);
%p_of_r.miss_o = zeros(n_r,1);
%p_of_r.miss_v = zeros(n_r,1);
%p_of_r.miss_vo = zeros(n_r,1);
%p_of_r.CR = zeros(n_r,1);
p_of_r.hit_o = zeros(n_r,1);
p_of_r.hit_v = zeros(n_r,1);
p_of_r.hit_vo = zeros(n_r,1);

% Create image of diameters
[n_row,n_col] = size(I_seg_id);
I_dia = zeros(n_row,n_col);
n_seg = numel(seg_dia);
for s = 1:n_seg;
    I_dia(I_seg_id==s) = seg_dia(s);
end

for r = 1:n_r;
    is_in_r = I_dia>r_array(r)&I_dia<=r_array(r+1); 
    n_FA = sum(sum(I_mask_vo2==1&is_in_r));
    %n_miss_v = sum(sum(I_mask_vo2==2&is_in_r));
    %n_miss_o = sum(sum(I_mask_vo2==3&is_in_r));
    %n_CR = sum(sum(I_mask_vo2==4&is_in_r));
    n_hit_v = sum(sum(I_mask_vo2==5&is_in_r));
    n_hit_o = sum(sum(I_mask_vo2==6&is_in_r));
    n_hit_vo = n_hit_v + n_hit_o;
    %n_miss_vo = n_miss_v + n_miss_o;
    
    n_v = sum(sum(I_ROI_is_v&is_in_r));
    n_o = sum(sum(I_ROI_is_o&is_in_r));
    n_vo = sum(sum(I_ROI_is_vo&is_in_r));
    n_n = sum(sum(~I_ROI_is_vo&is_in_r));

    p_of_r.FA(r) = n_FA./n_n; % Probability of FA 
    %p_of_r.miss_o(r) = n_miss_o./n_o; % Probability of miss given occluded 
    %p_of_r.miss_v(r) = n_miss_v./n_v; % Probability of miss given visible
    %p_of_r.miss_vo(r) = n_miss_vo./n_vo; % Probability of miss given visible and occluded
    %p_of_r.CR(r)  = n_CR./n_n; % Probability of Correct Rejection
    p_of_r.hit_o(r) = n_hit_o./n_o; % Probability of hits given occluded 
    p_of_r.hit_v(r) = n_hit_v./n_v; % Probability of hits given visible 
    p_of_r.hit_vo(r) = n_hit_vo./n_vo; % Probability of hits given visible and occluded    
end
    
figure;
hold on
xlim([0 r_array(end)]);
ylim([0 1]);
plot(r_array(1:end-1),p_of_r.FA,'-r','linewidth',2)
%plot(r_array(1:end-1),p_of_r.miss_v,'-b','linewidth',2)
%plot(r_array(1:end-1),p_of_r.miss_o,':c','linewidth',2)
%plot(r_array(1:end-1),p_of_r.CR,'-k','linewidth',2)
plot(r_array(1:end-1),p_of_r.hit_v,'-g','linewidth',2)
plot(r_array(1:end-1),p_of_r.hit_o,':y','linewidth',2)
%plot(r_array(1:end-1),p_of_r.miss_vo,':b','linewidth',2)
plot(r_array(1:end-1),p_of_r.hit_vo,':g','linewidth',2)
legend('False Alarm','Hit given visible','Hit given occluded','Hit given both');
xlabel(independent_var);
ylabel('Percentage of pixels (laser pulses)'); 
titlestr = sprintf('Classification accuracy as a function of %s',independent_var);
title(titlestr); 


