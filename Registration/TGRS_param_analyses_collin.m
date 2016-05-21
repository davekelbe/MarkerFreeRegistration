% Generate data with various covariance matrices
% Number of points (triangle)
n = 3;
% Number of dimensions
d = 2;
% Number of iterations 
n_it  = 100;
sxx = 1;
syy = linspace(.0001, 1, 100);
n_syy = numel(syy);

S_n = 0.01;

RMSE = zeros(n_syy, n_it);
T_normeig = zeros(n_syy, n_it);

for isyy = 1:n_syy
    for it = 1:n_it
        d=2;
        Sigma = [ sxx    0  ; ...
            0   syy(isyy)];
       % rng(42)
        %X = randn(n, d) * chol(Sigma);
        
        X = randn(n, d);
        X = bsxfun(@minus, X, mean(X));
        X = X * inv(chol(cov(X)));
        X = X * chol(Sigma);
        
        % Transform to 3D
        X = [X zeros(n,1)];
        d = 3;
        %{
        XT = [X; X(1,:)];
        figure;
        plot(XT(:,1), XT(:,2),'-or');
        axis equal
        %}
        P_noise = S_n.*randn(d,n);
        Y = X + P_noise;
        %{
        XT = [X; X(1,:)];
        figure;
        hold on
        plot(XT(:,1), XT(:,2),'-or');
        axis equal
        YT = [Y; Y(1,:)];
        plot(YT(:,1), YT(:,2),'-ob');
        %}
        
        [R,t] = rigid_transform_3D(X,Y);
        Xt = ((R*X')+ repmat(t,1,n))'; %Note change
        %{
        XT = [X; X(1,:)];
        figure;
        hold on
        plot(XT(:,1), XT(:,2),'-or');
        axis equal
        YT = [Y; Y(1,:)];
        plot(YT(:,1), YT(:,2),'-ob');
        XtT = [Xt; Xt(1,:)];
        plot(XtT(:,1), XtT(:,2),'-og');
        %}
        RMSE(isyy, it) = mean(sqrt(sum((X-Xt).^2,2)));
        temp = eig(cov(X));
        T_eig = temp(2:3);
        coll   = T_eig./(repmat(sum(T_eig),[1,2]))';
        T_normeig(isyy,it) = coll(1);
        foo = 1;
    end
end


T_normeigmean = mean(T_normeig,2);
RMSEmean = mean(RMSE,2);


figure;
plot(T_normeigmean, RMSEmean, '-or');
xlabel('Collinearity');
ylabel('Error');



