function [  ti ] = get_initial_transformation( info_plot1, info_plot2 )
% GET_REGISTRATION_SEQUENCE Returns inital translation
% 
%
%   (c) David Kelbe, Rochester Institute of Technology 

    t_arr = zeros(3,25);
    t_arr(:,1) = [-10 10 0]';     
    t_arr(:,2) = [-5 10 0]';     
    t_arr(:,3) = [0 10 0]';     
    t_arr(:,4) = [5 10 0]';     
    t_arr(:,5) = [10 10 0]';     
    t_arr(:,6) = [10 5 0]';     
    t_arr(:,7) = [5 5 0]';
    t_arr(:,8) = [0 5 0]';
    t_arr(:,9) = [-5 5 0]';
    t_arr(:,10) = [-10 5 0]';
    t_arr(:,11) = [-10 0 0]';
    t_arr(:,12) = [-5 0 0]'; 
    t_arr(:,13) = [0 0 0]';
    t_arr(:,14) = [5 0 0]';
    t_arr(:,15) = [10 0 0]';
    t_arr(:,16) = [10 -5 0]';
    t_arr(:,17) = [5 -5 0]';
    t_arr(:,18) = [0 -5 0]';
    t_arr(:,19) = [-5 -5 0]';
    t_arr(:,20) = [-10 -5 0]';
    t_arr(:,21) = [-10 -10 0]';
    t_arr(:,22) = [-5 -10 0]';
    t_arr(:,23) = [0 -10 0]';
    t_arr(:,24) = [5 -10 0]';
    t_arr(:,25) = [10 -10 0]';
    
    %{
    figure;
    axis([-15 15 -15 15]);
    hold on
    for i = 1:25;
    text(t_arr(1,i),t_arr(2,i),sprintf('%g',i));
    end
    %}

ti = t_arr(:,info_plot1) - t_arr(:,info_plot2);
ti = -ti;
end

