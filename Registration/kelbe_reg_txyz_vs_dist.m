function [  ] = kelbe_reg_txyz_vs_dist( G_dist,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Rotation vs. # nodes
%
%Gminpath = min(G_dist);
%Gmaxpath = max(G_dist); % change to automatically floor/ciel
dist_axes = 10:10:90;  
n_p = numel(dist_axes) - 1;
    
if options_matlabfig
    figure;
    hold on
    plot(G_dist,G_tx,'ok','markerfacecolor','r');
    plot(G_dist,G_ty,'ok','markerfacecolor','g');
    plot(G_dist,G_tz,'ok','markerfacecolor','b');
    grid on
    set(gca, 'xtick',tick);
    xlabel('Number of nodes');
    ylabel('Error [m]');
    legend_str{1} = 'tx';
    legend_str{2} = 'ty';
    legend_str{3} = 'tz';
    legend(legend_str);
    %title('Error in Translation');
    %filename = sprintf('Error-txyz_vs_dist_s%02.0f',info_site);
    %filepath_save = sprintf('%s%s.eps',path_save,filename);
    %saveas(gcf,filepath_save,'psc2')
    %}
end


if options_tikzfig;
    % Make tikz figure
    data_1 = cell(n_p,1);
    data_2 = cell(n_p,1);
    data_3 = cell(n_p,1);
    labels = cell(n_p,1);
    ctr = 1;
    for p = 1:Gmaxpath;
    is_valid = (G_dist>=dist_axes(d)) & (G_dist<dist_axes(d+1));
        if sum(is_valid);
            data_1{ctr} = G_tx(is_valid);
            data_2{ctr} = G_ty(is_valid);
            data_3{ctr} = G_tz(is_valid);
            labels{ctr} = sprintf('%g-%g',dist_axes(d), dist_axes(d+1));
            ctr = ctr + 1;
        end
    end
    
    tikzxlabel = 'Distance [m]';
    tikzylabel = 'Error [m]';
    legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};
    ymin = min([G_tx; G_ty; G_tz]);
    ymax = max([G_tx; G_ty; G_tz]);
    ylimval = [ymin ymax];
    
    filename = sprintf('Error-txyz_vs_dist_s%02.0f',info_site);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
    fid = fopen(filepath_tikz,'w+');
    if save_tikz;
        makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
            'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
            'legendlabels', legendlabels,...
            'boxsep',2.25,...
            'ylim', ylimval);
        fclose all;
    end
end

end

