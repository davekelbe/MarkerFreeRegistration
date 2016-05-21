function [ F ] = myfun( x )
%MATLAB lsqnonlin example 
%   Detailed explanation goes here

fprintf('*************');
fprintf('\n%6.6f\t\t%6.6f\n',x(1),x(2));
k = 1:10;
F = 2 + 2*k-exp(k*x(1))-exp(k*x(2));
fprintf('\nRMSE = %6.6f\n',sqrt(mean(F.^2)));

end

