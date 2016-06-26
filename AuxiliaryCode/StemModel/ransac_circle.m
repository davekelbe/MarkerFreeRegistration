function [ xc_best,yc_best,r_best, is_inlier_best,E_best ] = ransac_circle( data, t_error )
%RANSAC_CYL determines the parameters of a cylinder to fit data
%load('D:\Users\djk2312\Documents\data.mat');

% Initial parameters
n_data = size(data,2);
n_mss = 3;
p = 0.999;
k = inf;
i_max = 1000;
i = 1;
n_inliers_best = 0;

while i < k && i < i_max;
    
    % Choose a random Minimum Sample Set
    index = randperm(n_data);
    index = index(1:n_mss);
    
    % Rename variables for consistency
    Theta = [data(:,index(1)) data(:,index(2)) data(:,index(3))];
    X = data;
    
    % Fit circle
    [xc,yc,r] = CircleFitByKasa(Theta(1,:),Theta(2,:));
    
    % Determine the error
    E = sqrt(abs(sum((X - repmat([xc;yc],[1,size(X,2)])).^2,1)-r^2));
    
    % Find inliers
    is_inlier = (E<t_error);
    n_inliers = numel(is_inlier(is_inlier));
    
    % Update 
    i = i+1;
    if n_inliers > n_inliers_best;
        n_inliers_best = n_inliers;
        is_inlier_best = is_inlier;
        xc_best = xc;
        yc_best = yc;
        r_best = r;
        E_best = E;
       % omega = n_inliers/n_data;
       % k = log( 1 - p) / log( 1 - omega.^n_mss);
    end
    %
    
    % Plot the results
    %{
    figure;
    hold on;
    scatter(data(1,is_inlier),data(2,is_inlier),10,'b','filled');
    scatter(data(1,~is_inlier),data(2,~is_inlier),10,'r','filled');
    rectangle('position',[xc-r,yc-r,r*2,r*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','k','linewidth',2);
    view(0,90);
    %}
end
     % Plot the final results
    %{
    figure;
    hold on;
    scatter(data(1,is_inlier_best),data(2,is_inlier_best),10,'b','filled');
    scatter(data(1,~is_inlier_best),data(2,~is_inlier_best),10,'r','filled');
    rectangle('position',[xc_best-r_best,yc_best-r_best,r_best*2,r_best*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','k','linewidth',2);
    view(0,90);
    axis equal
    %}


