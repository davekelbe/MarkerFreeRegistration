function [ h ] = plot_triangle( tri )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

T = [tri; tri(1,:)];
h = plot(gca, T(:,1), T(:,2), '-k','linewidth', 2);

end

