% Test RANSAC inlier parameter


clear legend
noise = linspace(0,.6,25);
n_iter = 10000;
n_noise = numel(noise);

xmin = -25;
ymin = -25;
xmax = 25;
ymax = 25;
zmin = -2;
zmax = 2;


n_trees = 25;

e_Eucl_mean = zeros(n_noise, n_iter);
e_Eucl_max = zeros(n_noise, n_iter);


for n = 1:n_noise;
    for i = 1:n_iter;
        triIx = xmin + (xmax - xmin).*rand(n_trees,1);
        triIy = ymin + (ymax - ymin).*rand(n_trees,1);
        triIz = zmin + (zmax - zmin).*rand(n_trees,1);
        triI = [triIx triIy triIz]; 
        trinoise = noise(n)*randn(n_trees,3);
        triJ = triI + trinoise;
        
        ixrand = randperm(n_trees,3);
        [Rhat,that] = rigid_transform_3D(triJ(ixrand,:),triI(ixrand,:));
        triJt = (((Rhat*(triJ'))+ repmat(that,1,n_trees)))';
        
        e_Eucl = sum((triJt - triJ).^2,2);
        e_Eucl_mean(n,i) = mean(e_Eucl);
        %e_Eucl_max(n,i) = max(e_Eucl);

    end
end

e_Eucl_meanav = mean(e_Eucl_mean,2);
%e_Eucl_maxav = mean(e_Eucl_max,2);

[f, gof] = fit(noise',e_Eucl_meanav,'poly2');
temp = 0:.01:.08;
thresh_opt = f.p1*(temp.^2) + f.p2*temp;

p = polyfitZero(noise, e_Eucl_meanav, 2);
ixsm = (noise < 0.11);
plinear = polyfitZero(noise(ixsm), e_Eucl_meanav(ixsm), 1);

rmse_est = 0.0644;
thresh_opt = p(1)*(rmse_est.^2) + p(2)*rmse_est;
thresh_optlinear = plinear(1)*rmse_est;


figure;
scatter(noise, e_Eucl_meanav, 'filled');
hold on
%plot(f,noise,e_Eucl_meanav)
plot(noise,polyval(p,noise), '-k')
%plot(noise,polyval(plinear,noise), '-m')
%scatter(noise, e_Eucl_maxav, 'filled');
%legend('Mean', 'Max');
xlabel('Noise N(0,sigma) added to tie points [m]')
ylabel('Mean Euclidean distance between corresponding tie points');
hold on 
scatter(rmse_est,thresh_opt, '*g') 
legend('Data', 'Fit', 'Optimal Threshold for input RMSE = 0.0644 m');
plot([rmse_est rmse_est], [0 thresh_opt], '-g');
plot([0 rmse_est], [thresh_opt thresh_opt], '-g');

