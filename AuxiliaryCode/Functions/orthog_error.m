function [ Xfit ] = orthog_fit( X )
%ORTHOG_ERROR Compute the error orthogonal to points of number n x 2 
%   Detailed explanation goes here

%% Fitting an Orthogonal Regression Using Principal Components Analysis
% This example shows how to use Principal Components Analysis (PCA) to fit a 
% linear regression. PCA minimizes the perpendicular distances from the data 
% to the fitted model. This is the linear case of what is known as Orthogonal 
% Regression or Total Least Squares, and is appropriate when there is no 
% natural distinction between predictor and response variables, or when all 
% variables are measured with error. This is in contrast to the usual regression 
% assumption that predictor variables are measured exactly, and only the response 
% variable has an error component.
%
% For example, given two data vectors x and y, you can fit a line that
% minimizes the perpendicular distances from each of the points (x(i), y(i))
% to the line.  More generally, with p observed variables, you can fit an
% r-dimensional hyperplane in p-dimensional space (r < p).  The choice of r is
% equivalent to choosing the number of components to retain in PCA. It may be
% based on prediction error, or it may simply be a pragmatic choice to reduce
% data to a manageable number of dimensions.
%
% In this example, we fit a plane and a line through some data on three
% observed variables.  It's easy to do the same thing for any number of
% variables, and for any dimension of model, although visualizing a fit
% in higher dimensions would obviously not be straightforward.

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2012/04/20 19:46:39 $


%% Fitting a Line to 2-D Data
% First, we generate some trivariate normal data for the example.  Two of
% the variables are fairly strongly correlated.
plot(X(:,1),X(:,2),'bo');
grid on;
axis auto
axis square

%%
% Next, we fit a line to the data using PCA.  The coefficients for the first
% principal component defines the line.
% The second PC is orthogonal to the first, and its coefficients define the
% normal vector of the line.
[coeff,score,roots] = pca(X);
basis = coeff(:,1);
%%
normal = coeff(:,2);
%%
% That's all there is to the fit.  But let's look closer at the results, and
% plot the fit along with the data.

%%
% Because the first component explains as much of the variance in the data
% as is possible with one dimension, the line is the best 1-D linear
% approximation to the data.  Equivalently, the second component explains the
% least amount of variation in the data, and it is the error term in the
% regression.  The latent roots (or eigenvalues) from the PCA define the
% amount of explained variance for each component.
pctExplained = roots' ./ sum(roots);

%%
% The first two coordinates of the principal component scores give the
% projection of each point onto the plane, in the coordinate system of the
% plane.  To get the coordinates of the fitted points in terms of the original
% coordinate system, we multiply each PC coefficient vector by the
% corresponding score, and add back in the mean of the data.  The residuals
% are simply the original data minus the fitted points.
[n,p] = size(X);
meanX = mean(X,1);
Xfit = repmat(meanX,n,1) + score(:,1)*coeff(:,1)';
residuals = X - Xfit;

%%
% The equation of the fitted plane, satisfied by each of the fitted points in
% |Xfit|, is |([x1 x2 x3] - meanX)*normal = 0|.  The plane passes through the
% point |meanX|, and its perpendicular distance to the origin is
% |meanX*normal|. The perpendicular distance from each point in |X| to the
% plane, i.e., the norm of the residuals, is the dot product of each centered
% point with the normal to the plane.  The fitted plane minimizes the sum of
% the squared errors.
error = sqrt(residuals(:,1).^2 + residuals(:,2).^2);
sse = sum(error.^2);
mse = mean(error.^2);
rmse = sqrt(mse);
stdev = std(error);

% Fit 
fit_minx = min(X,1);
fit_maxx = max(X,1);
fit_miny = coeff(1,1)*fit_minx 


figure;
scatter(X(:,1),X(:,2),20,error,'filled')

%%

end

