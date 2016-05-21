% Transform all point clouds and save PLY for Paper Figure 

site_arr = [31];
for s = 1:numel(site_arr);
    site = site_arr(s);
    path_upper = 'D:\Users\djk2312\Documents\Harvard\reg\';
    path_mat = sprintf('%s%03.0f%s',path_upper,site,'\25\mat\');
    filepath_G_R_MST = sprintf('%sG_R_MST.mat',path_mat);
    filepath_G_t_MST = sprintf('%sG_t_MST.mat',path_mat);
    load(filepath_G_R_MST);
    load(filepath_G_t_MST);
    if isempty(G_R_MST{13})
            continue
        end
    path_out= sprintf('%s%03.0f%s',path_upper,site,'\13\ply\');
    n_S = 25;
    all_color = jet(n_S);
    all_color = all_color(randperm(n_S),:);
    for j = 1:n_S;
        fprintf('\n%d\n',j);
        if isempty(G_R_MST{j})
            continue
        end
        path_ply = sprintf('%s%03.0f%s%02.0f%s',path_upper,site,'\',j,'\ply\');
        filepath_ply = sprintf('%s%s%03.0f-%02.0f.ply',path_ply,'points_full_',site,j);
        path_paper = sprintf('%s%s',path_out,'Paper\');
        if ~exist(path_paper, 'dir');
            mkdir(path_paper);
        end
        filepath_out = sprintf('%s%sPtx-%02.0f.ply',path_out,'Paper\',j);
        if exist(filepath_out, 'file');
           % continue
        end
        [vertexJ, faces] = read_ply(filepath_ply);
        data_xJ = vertexJ(:,1);
        data_yJ = vertexJ(:,2);
        data_zJ = vertexJ(:,3);
        colorJ = vertexJ(:,4:6);
        xyz2t = (G_R_MST{13}'*G_R_MST{j}*[data_xJ data_yJ data_zJ]') + ...
            repmat(...
            G_R_MST{13}'*G_t_MST{j} - ...
            G_R_MST{13}'*G_t_MST{13},...
            [1,numel(data_xJ)]);
        colorout = colorJ;%.*repmat(all_color(j,:), [numel(data_xJ), 1]);
        write2ply(filepath_out, xyz2t', colorout);
        foo = 1;
    end
end