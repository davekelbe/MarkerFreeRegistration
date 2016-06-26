function [ z1,y1, is_inlier_best ] = ransac_3line( data, t_error )
%RANSAC_CYL determines the parameters of a cylinder to fit data
%load('D:\Users\djk2312\Documents\data.mat');

% Initial parameters
n_data = size(data,2);
n_mss = 2;
p = 0.999;
k = inf;
i_max = 2000;
i = 1;
n_inliers_best = 0;

while i < k && i < i_max;
    
    % Choose a random Minimum Sample Set
    index = randperm(n_data);
    index = index(1:n_mss);
    
    % Rename variables for consistency
    Theta = [data(:,index(1)) data(:,index(2))];
    X = data;
    
    % Determine the normal vector
    qproj = @(u) eye(length(u))-u*u'/(u'*u);
    a1 = Theta(:,1) - Theta(:,2);
    a1norm = a1/norm(a1);
    if a1norm(3)<0;
        a1norm = -a1norm;
        Theta = fliplr(Theta);
    end
    
    % Determine the error
    Xproj = qproj(a1norm)*X;
    centerproj = qproj(a1norm)*Theta(:,1);
    E  = sqrt(sum((Xproj - repmat(centerproj,[1,size(Xproj,2)])).^2,1));
    
    % Find inliers
    is_inlier = (E<t_error);
    n_inliers = numel(is_inlier(is_inlier));
    
    % Update 
    i = i+1;
    if n_inliers > n_inliers_best;
        n_inliers_best = n_inliers;
        is_inlier_best = is_inlier;
        z1 = Theta(:,2);
        y1 = Theta(:,1);
        omega = n_inliers/n_data;
        k = log( 1 - p) / log( 1 - omega.^n_mss);
    end
    %
    
    % Plot the results
    %{
    figure;
    hold on;
    plot3(Theta(1,:),Theta(2,:),Theta(3,:),'-k','linewidth',4) 
    scatter3(X(1,is_inlier),X(2,is_inlier),X(3,is_inlier),10,'b','filled')
    scatter3(X(1,~is_inlier),X(2,~is_inlier),X(3,~is_inlier),10,'r','filled')
    view(0,0);
    %}
end
    % Plot the final results
    %{
    figure;
    hold on;
    plot3([y1(1) z1(1)],[y1(2) z1(2)],[y1(3) z1(3)],'-k','linewidth',4) 
    scatter3(X(1,is_inlier_best),X(2,is_inlier_best),X(3,is_inlier_best),10,'b','filled')
    scatter3(X(1,~is_inlier_best),X(2,~is_inlier_best),X(3,~is_inlier_best),10,'r','filled')
    view(0,0);
    %}
    foo = 1;
    


