function [  ] = test_param( b1, sb1, value, alpha, dof )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% H0: b1 = value vs. H1: b1 ~=value;
t1 = (b1 -value)/sb1;
% Reject H0 when |t1| > t_a/2,n-2
t_a = tinv(1-(alpha/2),dof);
if abs(t1) > t_a; 
    fprintf('\nReject H0 %2.2f confidence BAD\n', 1-alpha);
else
    fprintf('\nFail to reject H0 at %2.2f confidence GOOD\n', 1-alpha);
end