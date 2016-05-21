function [  ] = testparam( x,y,btest )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    alpha = 0.05;

    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    %{
    figure; scatter(x,y); hold on
    xstep = min(x):max(x);
    ystep = polyval(p,xstep);
    hold on;
    plot(xstep, ystep);
    %}
    btestslope = btest(1); 
    btestint = btest(2);  
    n = numel(x);
    yfit = polyval(p,x);
    SSres = sum((y-yfit).^2);
    MSres = SSres./(n-2); % ok 
    Sxx = sum((x - mean(x)).^2);
    t0slope = (p(1) - btestslope)./sqrt(MSres./Sxx);
    seint = sqrt(MSres*((1./n) + (mean(x)^2./Sxx) ));
    t0int = (p(2) - btestint)./seint;
    % t_crit
    alphaup = 1-alpha/2;
    dof = numel(x) - 2;
    tcrit = tinv(alphaup,dof);
    
    
    
end

