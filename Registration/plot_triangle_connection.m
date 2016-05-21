function [  ] = plot_triangle_connection( TI, TJ )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for n = 1:3;
    plot(gca, [TI(n,1),TJ(n,1)],[TI(n,2),TJ(n,2)], 'Color', [.5 .5 .5],'linewidth', 2);
end

end

