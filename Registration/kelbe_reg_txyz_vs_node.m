function [  ] = kelbe_reg_txyz_vs_node( G_npath,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Rotation vs. # nodes
%
Gminpath = min(G_npath);
Gmaxpath = max(G_npath);
tick = Gminpath:Gmaxpath;
n_p = numel(tick);

if options_matlabfig
    figure;
    hold on
    plot(G_npath,G_tx,'ok','markerfacecolor','r');
    plot(G_npath,G_ty,'ok','markerfacecolor','g');
    plot(G_npath,G_tz,'ok','markerfacecolor','b');
    grid on
    set(gca, 'xtick',tick);
    xlabel('Number of nodes');
    ylabel('Error [m]');
    legend_str{1} = 'tx';
    legend_str{2} = 'ty';
    legend_str{3} = 'tz';
    legend(legend_str);
    %title('Error in Translation');
    %filename = sprintf('Error-txyz_vs_nodes_s%02.0f',info_site);
    %filepath_save = sprintf('%s%s.eps',path_save,filename);
    %saveas(gcf,filepath_save,'psc2')
    %}
end

if options_tikzfig
    % Make tikz figure
    data_1 = cell(n_p,1);
    data_2 = cell(n_p,1);
    data_3 = cell(n_p,1);
    labels = cell(n_p,1);
    ctr = 1;
    for p = 1:Gmaxpath;
        is_valid = (G_npath==p);
        if sum(is_valid);
            data_1{ctr} = G_tx(is_valid);
            data_2{ctr} = G_ty(is_valid);
            data_3{ctr} = G_tz(is_valid);
            labels{ctr} = sprintf('%g',p);
            ctr = ctr + 1;
        end
    end
    
    tikzxlabel = 'Number of nodes';
    tikzylabel = 'Error [m]';
    legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};
    ymin = min([G_tx; G_ty; G_tz]);
    ymax = max([G_tx; G_ty; G_tz]);
    ylimval = [ymin ymax];
    
    filename = sprintf('Error-txyz_vs_nodes_s%02.0f',info_site);
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

