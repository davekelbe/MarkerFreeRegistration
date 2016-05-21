function [  ti ] = get_initial_transformation( info_plot1, info_plot2 )
% GET_REGISTRATION_SEQUENCE Returns inital translation
% 
%
%   (c) David Kelbe, Rochester Institute of Technology 

    t1_arr = zeros(3,25);
    t1_arr(:,1) = [-10 10 0]';     
    t1_arr(:,2) = [-5 10 0]';     
    t1_arr(:,3) = [0 10 0]';     
    t1_arr(:,4) = [5 10 0]';     
    t1_arr(:,5) = [10 10 0]';     
    t1_arr(:,6) = [10 5 0]';     
    t1_arr(:,7) = [5 5 0]';
    t1_arr(:,8) = [0 5 0]';
    t1_arr(:,9) = [-5 5 0]';
    t1_arr(:,10) = [-10 5 0]';
    t1_arr(:,11) = [-10 0 0]';
    t1_arr(:,12) = [-5 0 0]'; 
    t1_arr(:,13) = [0 0 0]';
    t1_arr(:,14) = [5 0 0]';
    t1_arr(:,15) = [10 0 0]';
    t1_arr(:,16) = [10 -5 0]';
    t1_arr(:,17) = [5 -5 0]';
    t1_arr(:,18) = [0 -5 0]';
    t1_arr(:,19) = [-5 -5 0]';
    t1_arr(:,20) = [-10 -5 0]';
    t1_arr(:,21) = [-10 -10 0]';
    t1_arr(:,22) = [-5 -10 0]';
    t1_arr(:,23) = [0 -10 0]';
    t1_arr(:,24) = [5 -10 0]';
    t1_arr(:,25) = [10 -10 0]';

foo = 1;

%ti = -ti;
end

