function [ z1_ext,y1_ext, is_inlier_best,E_best,status ] = ransac_2line( data, t_error )
%RANSAC_2LINE determines the parameters of a line fit to 2D data 
%
%   Inputs:
%       data            - 2 x m array of points 
%       t_error         - Error threshold for inliers
%       
%   Outliers
%       z1_ext          - Lower point (extended to data) of line
%       y1_ext          - Upper point (extended to data) of line
%       is_inlier_best  - Logical array to inliers
%       E_best          - Error of points 
%
%   Example
%       % number of points
%       N = 300;
%       % inilers percentage
%       p = 0.15;
%       % noise
%       sigma = 0.05;
%       % line parameters
%       m = -0.8;
%       q = 0.3;
%
%       % generate a set of points 
%       Ni = round(p*N);
%       No = N-Ni;
% 
%       % inliers
%       X1i = 3 * 2*(rand(1, Ni)-0.5);
%       X2i = m*X1i+q;
% 
%       % and add some noise
%       X1i = X1i + sigma*randn(1, Ni);
%       X2i = X2i + sigma*randn(1, Ni);
% 
%       % outliers
%       X1o = 3 * 2*(rand(1, No)-0.5);
%       X2o = 3 * 2*(rand(1, No)-0.5);
% 
%       % Combine 
%       X1 = [X1i X1o];
%       X2 = [X2i X2o];
%
%       data = [X1; X2];
%       t_error = .1
%       [ z1,y1, is_inlier,E ] = ransac_2line( data, t_error );
%       figure;
%       scatter(data(1,is_inlier),data(2,is_inlier), 10, 'b','filled');
%       hold on
%       scatter(data(1,~is_inlier),data(2,~is_inlier),10, 'r','filled');
%       plot([z1(1) y1(1)]',[z1(2) y1(2)]','-b','linewidth',2)
%
%
%   (C) David Kelbe, Rochester Institute of Technology 

% Initial parameters
n_data = size(data,2);
n_mss = 2;
%p = 0.999;
%k = inf;
i_max = 500;
i = 1;
n_inliers_best = 0;

status_inlier = zeros(i_max,1);
status_iteration = (1:i_max)';
%{
figure; 
scatter(data(1,:),data(2,:),'r','filled')
%}

    qproj = @(u) eye(length(u))-u*u'/(u'*u);
    uproj = @(u) u*u'/(u'*u);
    
%while i < k && i < i_max;
while i < i_max;    
    % Choose a random Minimum Sample Set
    index = randperm(n_data);
    index = index(1:n_mss);
    
    % Rename variables for consistency
    Theta = [data(:,index(1)) data(:,index(2))];
    X = data;
    
    % Determine the normal vector
    a1 = Theta(:,1) - Theta(:,2);
    a1norm = a1/norm(a1);
    
    % Determine the error
    Xproj = qproj(a1norm)*X;
    centerproj = qproj(a1norm)*Theta(:,1);
    E  = sqrt(sum((Xproj - repmat(centerproj,[1,size(Xproj,2)])).^2,1));
    
    %E =  abs(sum((X - repmat([xc;yc],[1,size(X,2)])).^2,1)-r^2);
    %E =  abs(sqrt(sum((X - repmat([xc;yc],[1,size(X,2)])).^2,1))-r);
     
    % Find inliers
    is_inlier = (E<t_error);
    n_inliers = numel(is_inlier(is_inlier));
    
    % Update 
    i = i+1;
    if n_inliers > n_inliers_best;
        
        status_inlier(i-1) = n_inliers;
        n_inliers_best = n_inliers;
        is_inlier_best = is_inlier;
        z1 = Theta(:,2);
        y1 = Theta(:,1);
        %omega = n_inliers/n_data;
        %k = log( 1 - p) / log( 1 - omega.^n_mss);
        E_best = E;
    end
    %
    
    % Plot the results
    %{
    figure;
    hold on;
    plot(Theta(1,:),Theta(2,:),'-k','linewidth',4) 
    scatter(X(1,is_inlier),X(2,is_inlier),10,'b','filled')
    scatter(X(1,~is_inlier),X(2,~is_inlier),10,'r','filled')
    foo = 1
    %}
end

 % Extend to base and top
 a1norm = z1-y1;
 T = uproj(a1norm);
 Xproj =T*X;
 [~,Imax] = max(Xproj(2,:));
 [~,Imin] = min(Xproj(2,:));
 y_min = X(2,Imin);
 y_max = X(2,Imax);
 m = (z1(2) - y1(2))/(z1(1) - y1(1));
 b = y1(2) - m*y1(1);
 x_max = (y_max-b)/m;
 x_min = (y_min-b)/m;
y1_ext =[ x_min y_min]';
z1_ext = [x_max y_max]';

% Status_inlier 
is_valid = (status_inlier)>0;
status.inlier = status_inlier(is_valid);
status.iteration = status_iteration(is_valid);

 
    % Plot the final results
    %{
    figure;
    hold on;
    plot([y1(1) z1(1)],[y1(2) z1(2)],'-m','linewidth',4) 
    plot([y1_ext(1) z1_ext(1)],[y1_ext(2) z1_ext(2)],'-k','linewidth',2) 
    scatter(X(1,is_inlier_best),X(2,is_inlier_best),40,'b','filled')
    scatter(X(1,~is_inlier_best),X(2,~is_inlier_best),40,'r','filled')
    foo = 1
    %}
    


