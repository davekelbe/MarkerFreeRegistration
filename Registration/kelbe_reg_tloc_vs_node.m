function [  ] = kelbe_reg_tloc_vs_node( all_npath,all_tree_exy, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Rotation vs. # nodes
%

if options_matlabfig
    clear legend_str
    figure;
    hold on
    plot(all_npath, all_tree_exy,'+','markeredgecolor','k');
    axisval = axis;
    legend_str{1} = 'x';
    legend_str{2} = 'y';
    %legend_str{3} = 'z';
    legend(legend_str);
    xlabel('Numbre of nodes');
    ylabel('Error in xy [m]');
    %title('Tree position accuracy after registration');
    %filename = sprintf('Error-tloc_vs_nodes_s%02.0f',info_site);
    %filepath_save = sprintf('%s%s.eps',path_save,filename);
    %saveas(gcf,filepath_save,'psc2')
    %}
end

if options_tikzfig
    n_p = max(all_npath) - min(all_npath);
    data1 = cell(n_p,1);
    data2 = cell(n_p,1);
    data3 = cell(n_p,1);
    labels = cell(n_p,1);
    data_isvalid = false(n_p,1);
    ctr = 1;
    for p = 1:max(all_npath);
        is_valid = (all_npath==p);
        if sum(is_valid);
            data_1{p} = all_tree_ex(is_valid);
            data_2{p} = all_tree_ey(is_valid);
            data_3{p} = all_tree_ez(is_valid);
            data_isvalid(p) = true;
            labels{ctr} = sprintf('%g',p);
            ctr = ctr + 1;
        end
    end
    data_1 = data_1(data_isvalid);
    data_2 = data_2(data_isvalid);
    data_3 = data_3(data_isvalid);
    
    tikzxlabel = 'Number of nodes';
    tikzylabel = 'Error';
    legendlabels = {'$\gls{tree:xg}(x) - \hat{\gls{tree:xg}}(x)$',...
        '$\gls{tree:xg}(y) - \hat{\gls{tree:xg}}(y)$',...
        '$\gls{tree:xg}(z) - \hat{\gls{tree:xg}}(z)$'};
    
    ymin = min([all_tree_ex; all_tree_ey; all_tree_ez]);
    ymax = max([all_tree_ex; all_tree_ey; all_tree_ez]);
    ylimval = [ymin ymax];
    
    filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
    fid = fopen(filepath_tikz,'w');
    makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
        'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
        'legendlabels', legendlabels,...
        'boxsep',2.25,...
        'ylim', ylimval);
    fclose all;
end

end

