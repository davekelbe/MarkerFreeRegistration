function [  ] = kelbe_reg_product_pcerror( n_S,G_npath,G_path,P_plot,...
    filepath_ply,loop_R,loop_t)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

for g = 1:n_S;
    for p = 1:G_npath(g);
        fprintf('\nWriting ply error\n');
        path = G_path{g}(p,:);
        i = path(1);
        tmp_plot = P_plot(i);
        [vertex, ~] = read_ply(filepath_ply{i});
        data_xyz = vertex(:,1:3)';
       % color_ply = vertex(:,4:6);
        data_xyzhat = loop_R{g}{p}*data_xyz+...
            repmat(loop_t{g}{p},[1,size(vertex,1)]);
        data_e = abs(mean((data_xyz - data_xyzhat),1));
        
        figure;
        xbin = linspace(0,0.05,100);
        count = hist(data_e, xbin);
        plot(xbin,count,'-x');
        color_e = vec2cmap2(data_e,'jet', 0,.05);
        filepath_ply_e = sprintf('%serror_%03.0f-%02.0f.ply', ...
        path_ply{j}, info_site, P_plot(i));
        write2ply(filepath_ply_e,data_xyz', color_e);
    end
end

for j = 1:n_S;
    [vertex, ~] = read_ply(filepath_ply{j});

    filepath_ply_reg{j} = sprintf('%spoints_full_%03.0f-%02.0f-%02.0f.ply', ...
        path_ply{j}, info_site, P_plot(j),P_plot(1));
    xyz2t = (match_Reff{1,j}*[data_x2 data_y2 data_z2]') + ...
        repmat(match_teff{1,j},1,numel(data_x2));
    write2ply(filepath_ply_reg{j},xyz2t', color);
end


end

