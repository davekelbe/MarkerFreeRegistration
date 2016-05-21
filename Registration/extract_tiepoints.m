function [  ] = extract_tiepoints( path_up, info_experiment, info_suffix, plot_register, info_sites, info_plots  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

info_slash = '/';

n_site = numel(info_sites);
n_plot = numel(info_plots);
for s = 1:n_site;
    for p = 1:n_plot;
        info_site = info_sites(s);
        info_plot = info_plots(p);
        
        path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
            info_experiment, info_slash, info_suffix,info_slash,info_site, info_slash,info_plot,info_slash);
        path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
        
        filepath_seg_z = sprintf('%s%s%s',path_mat,'seg_z','.mat');
        filepath_seg_y = sprintf('%s%s%s',path_mat,'seg_y','.mat');
        filepath_seg_r = sprintf('%s%s%s',path_mat,'seg_r','.mat');
        
        
        filepath_dem_qx = sprintf('%s%s%s',path_mat,'dem_qx','.mat');
        filepath_dem_qy = sprintf('%s%s%s',path_mat,'dem_qy','.mat');
        filepath_dem_qz = sprintf('%s%s%s',path_mat,'dem_qz','.mat');
        
        load(filepath_dem_qx);
        load(filepath_dem_qy);
        load(filepath_dem_qz);
        
        load(filepath_seg_z);
        load(filepath_seg_y);
        load(filepath_seg_r);
        
        n_tree = size(seg_r,1);
        seg_PLCS = nan(n_tree,3);
        seg_Prad = nan(n_tree,1);
        % Calculate intersection of axis with terrain for each segment
        
        dem_qz_nonan = inpaint_nans(dem_qz);
       % tree_dem = nan(n_tree,3);
       % [nx, ny] = size(dem_qx);
        %xv = 1:nx-1;
        %yv = 1:ny-1;
        %[xxv yyv] = meshgrid(xv, yv);
        %xxv = xxv(:);
        %yyv = yyv(:);
        [ vertices, indices ] = ply_qxqyqz_tri( dem_qx,dem_qy,dem_qz_nonan );
        for t = 1:n_tree
            %fprintf('\nTree %d of %d\n', t, n_tree)
            line = [seg_z(t,:); seg_y(t,:) - seg_z(t,:)]';
            %line = [tree_sm(t).loc(:,2); 10*(tree_sm(t).loc(:,1) - tree_sm(t).loc(:,2)) ]';
            inters = intersectLineMesh3d(line, vertices, indices);
            if ~isempty(inters)
                seg_PLCS(t,:) = inters;
            end
        end
        
        % Linear interpolation to determine stem diameter at ground 
        for t = 1:n_tree;
            pfitval = polyfit([seg_y(t,3) seg_z(t,3)], seg_r(t,:),1);
            rfit = polyval(pfitval,seg_PLCS(t,3));
            seg_Prad(t) = rfit;
        end
        
        filepath_seg_PLCS = sprintf('%s%s%s',path_mat,'seg_PLCS','.mat');
        filepath_seg_Prad = sprintf('%s%s%s',path_mat,'seg_Prad','.mat');
        save(filepath_seg_PLCS, 'seg_PLCS');
        save(filepath_seg_Prad, 'seg_Prad');
      
    end
    
end

