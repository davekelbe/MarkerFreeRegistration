function [  ] = kelbe_reg_txyz_vs_mintie( G_nmatch,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Rotation vs. # min tie points
%
G_min_nmatch = cellfun(@min,G_nmatch);
G_meant = sqrt(G_tx.^2 + G_ty.^2 + G_tz.^2);
min_m = min(G_min_nmatch);
max_m = max(G_min_nmatch);
m_axes = min_m:1:max_m;
n_m = numel(m_axes);

if options_matlabfig
    figure;
    clear legend_str
    hold on
    n_maxpath = max(G_npath);
    color_npath = jet(n_maxpath);
    ctr = 1;
    for i = 1:n_maxpath;
        is_valid = (G_npath == i);
        if sum(is_valid);
            plot(G_min_nmatch(is_valid),G_meant(is_valid),'ok','markerfacecolor',color_npath(i,:));
            legend_str{ctr} = sprintf('%g nodes',i);
            ctr = ctr + 1;
        end
    end
    grid on
    xlabel('Minimum number of tie points');
    ylabel('Root Square Error [m]');
    legend(legend_str);
    %title('Error in rotation');
    filepath_save = sprintf('%sError-t_vs_nties_s%02.0f.eps',path_save, info_site);
    saveas(gcf,filepath_save,'psc2')
end


if options_tikzfig;
    % Make tikz figure
    data_1 = cell(n_m,1);
    data_2 = cell(n_m,1);
    data_3 = cell(n_m,1);
    labels = cell(n_m,1);
    for m = 1:n_m;
        labels{m} = sprintf('%g',m_axes(m));
        is_valid = (G_min_nmatch==m_axes(m));
        if sum(is_valid);
            data_1{m} = G_rx(is_valid);
            data_2{m} = G_ry(is_valid);
            data_3{m} = G_rz(is_valid);
        end
    end
    
    tikzxlabel = 'Minimum number of tie points';
    tikzylabel = 'Error [m]';
    legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};
    ymin = min([G_tx; G_ty; G_tz]);
    ymax = max([G_tx; G_ty; G_tz]);
    ylimval = [ymin ymax];
    
    filename = sprintf('Error-txyz_vs_nties_s%02.0f',info_site);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
    fid = fopen(filepath_tikz,'w+');
    makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
        'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
        'legendlabels', legendlabels,...
        'boxsep',2.25,...
        'ylim', ylimval);
    fclose all;
end

end

