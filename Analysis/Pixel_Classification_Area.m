function [a_of_r, pa_of_r ] = Pixel_Classification_Area( I12r,I12area,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,independent_var )
%%
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r_step = 2;
r_array = min(0,min(I12r(:))):r_step:max(I12r(:))+r_step;
n_r = numel(r_array)-1;
pa_of_r.FA = zeros(n_r,1);
pa_of_r.miss_o = zeros(n_r,1);
pa_of_r.miss_v = zeros(n_r,1);
pa_of_r.miss_vo = zeros(n_r,1);
pa_of_r.CR = zeros(n_r,1);
pa_of_r.hit_o = zeros(n_r,1);
pa_of_r.hit_v = zeros(n_r,1);
pa_of_r.hit_vo = zeros(n_r,1);

for r = 1:n_r;
    is_in_r = I12r>r_array(r)&I12r<=r_array(r+1); 
    a_of_r.FA = sum(sum(I12area(I_mask_vo2==1&is_in_r)));
    a_of_r.miss_v = sum(sum(I12area(I_mask_vo2==2&is_in_r)));
    a_of_r.miss_o = sum(sum(I12area(I_mask_vo2==3&is_in_r)));
    a_of_r.CR = sum(sum(I12area(I_mask_vo2==4&is_in_r)));
    a_of_r.hit_v = sum(sum(I12area(I_mask_vo2==5&is_in_r)));
    a_of_r.hit_o = sum(sum(I12area(I_mask_vo2==6&is_in_r)));
    a_of_r.hit_vo = a_of_r.hit_v + a_of_r.hit_o;
    a_of_r.miss_vo = a_of_r.miss_v + a_of_r.miss_o;
    
    n_v = sum(sum(I12area(I_ROI_is_v&is_in_r)));
    n_o = sum(sum(I12area(I_ROI_is_o&is_in_r)));
    n_vo = sum(sum(I12area(I_ROI_is_vo&is_in_r)));
    n_n = sum(sum(I12area(~I_ROI_is_vo&is_in_r)));
    
    pa_of_r.FA(r) = a_of_r.FA./n_n; % Probability of FA 
    pa_of_r.miss_o(r) = a_of_r.miss_o./n_o; % Probability of miss given occluded 
    pa_of_r.miss_v(r) = a_of_r.miss_v./n_v; % Probability of miss given visible
    pa_of_r.miss_vo(r) = a_of_r.miss_vo./n_vo; % Probability of miss given visible and occluded
    pa_of_r.CR(r)  = a_of_r.CR./n_n; % Probability of Correct Rejection
    pa_of_r.hit_o(r) = a_of_r.hit_o./n_o; % Probability of hits given occluded 
    pa_of_r.hit_v(r) = a_of_r.hit_v./n_v; % Probability of hits given visible 
    pa_of_r.hit_vo(r) = a_of_r.hit_vo./n_vo; % Probability of hits given visible and occluded    

end

    is_nan =  isnan(pa_of_r.hit_o);
    last_index = (is_nan - circshift(is_nan,[1,0])==1);
      if last_index<5; last_index = numel(is_nan);end
%{
figure;
hold on
xlim([r_array(1) r_array(last_index)]);
plot(r_array(1:end-1),pa_of_r.FA,'-r','linewidth',2)
plot(r_array(1:end-1),pa_of_r.miss_v,'-b','linewidth',2)
%plot(r_array(1:end-1),pa_of_r.miss_o,':c','linewidth',2)
plot(r_array(1:end-1),pa_of_r.CR,'-k','linewidth',2)
plot(r_array(1:end-1),pa_of_r.hit_v,'-g','linewidth',2)
%plot(r_array(1:end-1),pa_of_r.hit_o,':y','linewidth',2)
%plot(r_array(1:end-1),pa_of_r.miss_vo,':b','linewidth',2)
%plot(r_array(1:end-1),pa_of_r.hit_vo,':g','linewidth',2)
%legend('False Alarm','Miss given visible','Miss given occluded',...
%    'Correct Rejection','Hit given visible','Hit given occluded',...
%    'Miss given both', 'Hit given both');
legend('False Alarm','Miss',...
    'Correct Rejection','Hit');
xlabel(independent_var);
ylabel('Percentage of area'); 
titlestr = sprintf('Classification accuracy as a function of %s',independent_var);
title(titlestr); 
      %}


