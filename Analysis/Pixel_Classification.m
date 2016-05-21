function [n_of_r, p_of_r ] = Pixel_Classification( I12r,I_mask_vo2,I_ROI_is_v,I_ROI_is_o,I_ROI_is_vo,independent_var )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r_step = 0.25;
r_array = min(0,min(I12r(:))):r_step:max(I12r(:))+r_step;
n_r = numel(r_array)-1;

n_of_r.v = zeros(n_r,1);
n_of_r.o = zeros(n_r,1);
n_of_r.vo = zeros(n_r,1);
n_of_r.n = zeros(n_r,1);


n_of_r.FA = zeros(n_r,1);
n_of_r.miss_o = zeros(n_r,1);
n_of_r.miss_v = zeros(n_r,1);
n_of_r.miss_vo = zeros(n_r,1);
n_of_r.CR = zeros(n_r,1);
n_of_r.hit_o = zeros(n_r,1);
n_of_r.hit_v = zeros(n_r,1);
n_of_r.hit_vo = zeros(n_r,1);

p_of_r.FA = zeros(n_r,1);
p_of_r.miss_o = zeros(n_r,1);
p_of_r.miss_v = zeros(n_r,1);
p_of_r.miss_vo = zeros(n_r,1);
p_of_r.CR = zeros(n_r,1);
p_of_r.hit_o = zeros(n_r,1);
p_of_r.hit_v = zeros(n_r,1);
p_of_r.hit_vo = zeros(n_r,1);


n_of_r.x = (r_array(1:end-1))';
for r = 1:n_r;
    is_in_r = I12r>r_array(r)&I12r<=r_array(r+1); 
    n_of_r.FA(r) = sum(sum(I_mask_vo2==1&is_in_r));
    n_of_r.miss_v(r) = sum(sum(I_mask_vo2==2&is_in_r));
    n_of_r.miss_o(r) = sum(sum(I_mask_vo2==3&is_in_r));
    n_of_r.CR(r) = sum(sum(I_mask_vo2==4&is_in_r));
    n_of_r.hit_v(r) = sum(sum(I_mask_vo2==5&is_in_r));
    n_of_r.hit_o(r) = sum(sum(I_mask_vo2==6&is_in_r));
    n_of_r.hit_vo(r) = n_of_r.hit_v(r) + n_of_r.hit_o(r);
    n_of_r.miss_vo(r) = n_of_r.miss_v(r) + n_of_r.miss_o(r);
    
    n_of_r.v(r) = sum(sum(I_ROI_is_v&is_in_r));
    n_of_r.o(r) = sum(sum(I_ROI_is_o&is_in_r));
    n_of_r.vo(r) = sum(sum(I_ROI_is_vo&is_in_r));
    n_of_r.n(r) = sum(sum(~I_ROI_is_vo&is_in_r));
    
    p_of_r.FA(r) = n_of_r.FA(r)./n_of_r.n(r); % Probability of FA 
    p_of_r.miss_o(r) = n_of_r.miss_o(r)./n_of_r.o(r); % Probability of miss given occluded 
    p_of_r.miss_v(r) = n_of_r.miss_v(r)./n_of_r.v(r); % Probability of miss given visible
    p_of_r.miss_vo(r) = n_of_r.miss_vo(r)./n_of_r.vo(r); % Probability of miss given visible and occluded
    p_of_r.CR(r)  = n_of_r.CR(r)./n_of_r.n(r); % Probability of Correct Rejection
    p_of_r.hit_o(r) = n_of_r.hit_o(r)./n_of_r.o(r); % Probability of hits given occluded 
    p_of_r.hit_v(r) = n_of_r.hit_v(r)./n_of_r.v(r); % Probability of hits given visible 
    p_of_r.hit_vo(r) = n_of_r.hit_vo(r)./n_of_r.vo(r); % Probability of hits given visible and occluded    
end

  is_nan =  isnan(p_of_r.hit_o);
  last_index = (is_nan - circshift(is_nan,[1,0])==1);
  if last_index<5; last_index = numel(is_nan);end
  %{
figure;
hold on
xlim([r_array(1) r_array(last_index)]);
ylim([0 1]);
plot(r_array(1:end-1),p_of_r.FA,'-r','linewidth',2)
plot(r_array(1:end-1),p_of_r.miss_v,'-b','linewidth',2)
%plot(r_array(1:end-1),p_of_r.miss_o,':c','linewidth',2)
plot(r_array(1:end-1),p_of_r.CR,'-k','linewidth',2)
plot(r_array(1:end-1),p_of_r.hit_v,'-g','linewidth',2)
%plot(r_array(1:end-1),p_of_r.hit_o,':y','linewidth',2)
%plot(r_array(1:end-1),p_of_r.miss_vo,':b','linewidth',2)
%plot(r_array(1:end-1),p_of_r.hit_vo,':g','linewidth',2)
%legend('False Alarm','Miss given visible','Miss given occluded',...
%    'Correct Rejection','Hit given visible','Hit given occluded',...
%    'Miss given both', 'Hit given both');
legend('False Alarm','Miss',...
    'Correct Rejection','Hit');
xlabel(independent_var);
ylabel('Percentage of pixels (laser pulses)'); 
titlestr = sprintf('Classification accuracy as a function of %s',independent_var);
title(titlestr); 
  %}


