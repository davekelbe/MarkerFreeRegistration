%% Label results good or bad
path_figures = 'Z:\Desktop\registration\figures\';
filename_csv = 'D:\Users\djk2312\Documents\thisshouldbecyclone\2012-08-Harvard.csv';
fid = fopen(filename_csv);
%C = textscan(fid, '%u8%u8%s%s%s%u16%u16%u16%u16%u16%u16%u16%u16%u16',...
%    'delimiter', ',',...
%    'headerlines',8);
C = textscan(fid, '%u8%u8%s%s%s%s%s%s%s%s%s%s%s%s',...
    'delimiter', ',',...
    'headerlines',8);
% 1 Site
% 2 Plot
% 3 Lidar Filename
% 4 Date
% 5 Lidar Orientation

n_scans = numel(C{1});
warning('off', 'arguments:exteriordata');

info_exp = 'Harvard';
info_suffix = 'reg';
info_slash = '\';

site_unique = unique(C{1});
n_site = numel(site_unique);

x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
[x_sph, y_sph, z_sph] = sphere;
exclude_site = [];

t_r = 10;
t_t = 1;

for s = 5:5;%1:n_site;
    info_site = site_unique(s);
    if ismember(exclude_site, info_site);
        continue
    end
    path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
        info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
    
    path_mat = sprintf('%s%02.0f%smat%s',path_site,25,info_slash,info_slash);
    filepath_nit = sprintf('%s%s',path_mat, 'match_nit.mat');
    if ~exist(filepath_nit, 'file');
        continue
    end
    % Load matches
    load(filepath_nit);
    filepath_match12_RMSE = sprintf('%s%s',path_mat, 'match12_RMSE.mat');
    load(filepath_match12_RMSE);
    filepath_match12_R = sprintf('%s%s',path_mat, 'match12_R.mat');
    load(filepath_match12_R);
    filepath_match12_t = sprintf('%s%s',path_mat, 'match12_t.mat');
    load(filepath_match12_t);
    filepath_match1_R = sprintf('%s%s',path_mat, 'match1_R.mat');
    load(filepath_match1_R);
    filepath_match1_t = sprintf('%s%s',path_mat, 'match1_t.mat');
    load(filepath_match1_t);
    filepath_match2_R = sprintf('%s%s',path_mat, 'match2_R.mat');
    load(filepath_match2_R);
    filepath_match2_t = sprintf('%s%s',path_mat, 'match2_t.mat');
    load(filepath_match2_t);
    filepath_match_plotI = sprintf('%s%s',path_mat, 'match_plotI.mat');
    load(filepath_match_plotI);
    filepath_match_plotJ = sprintf('%s%s',path_mat, 'match_plotJ.mat');
    load(filepath_match_plotJ);
    
    match1_tx =  cellfun(@decompose_translation_tx,match1_t);
    match1_ty =  cellfun(@decompose_translation_ty,match1_t);
    match1_tz =  cellfun(@decompose_translation_tz,match1_t);
    match2_tx =  cellfun(@decompose_translation_tx,match2_t);
    match2_ty =  cellfun(@decompose_translation_ty,match2_t);
    match2_tz =  cellfun(@decompose_translation_tz,match2_t);
    match_tdiff = sqrt((match1_tx + match2_tx).^2 + ...
        (match1_ty + match2_ty).^2 + (match1_tz + match2_tz).^2);
    match_txydiff = sqrt((match1_tx + match2_tx).^2 + ...
        (match1_ty + match2_ty).^2);
    %{
%     filepath_ran_rx = sprintf('%s%s',path_mat, 'ran_rx.mat');
%     load(filepath_ran_rx);
%     filepath_ran_ry = sprintf('%s%s',path_mat, 'ran_ry.mat');
%     load(filepath_ran_ry);
%     filepath_ran_rz = sprintf('%s%s',path_mat, 'ran_rz.mat');
%     load(filepath_ran_rz);
%     filepath_ran_tx = sprintf('%s%s',path_mat, 'ran_tx.mat');
%     load(filepath_ran_tx);
%     filepath_ran_ty = sprintf('%s%s',path_mat, 'ran_ty.mat');
%     load(filepath_ran_ty);
%     filepath_ran_tz = sprintf('%s%s',path_mat, 'ran_tz.mat');
%     load(filepath_ran_tz);
%     filepath_ran_conf = sprintf('%s%s',path_mat, 'ran_conf.mat');
%     load(filepath_ran_conf);
%     filepath_ran_R = sprintf('%s%s',path_mat, 'ran_R.mat');
%     load(filepath_ran_R);
%     filepath_ran_t = sprintf('%s%s',path_mat, 'ran_t.mat');
%     load(filepath_ran_t);
%     filepath_ran_t = sprintf('%s%s',path_mat, 'ran_t.mat');
%     load(filepath_ran_t);
        %}
        
        % Load graph estimates
        filepath_G_R = sprintf('%s%s',path_mat, 'G_R.mat');
        load(filepath_G_R);
        filepath_G_t = sprintf('%s%s',path_mat, 'G_t.mat');
        load(filepath_G_t);
        filepath_G_tx = sprintf('%s%s',path_mat, 'G_tx.mat');
        load(filepath_G_tx);
        filepath_G_ty = sprintf('%s%s',path_mat, 'G_ty.mat');
        load(filepath_G_ty);
        filepath_G_tz = sprintf('%s%s',path_mat, 'G_tz.mat');
        load(filepath_G_tz);
        filepath_G_rx = sprintf('%s%s',path_mat, 'G_rx.mat');
        load(filepath_G_rx);
        filepath_G_ry = sprintf('%s%s',path_mat, 'G_ry.mat');
        load(filepath_G_ry);
        filepath_G_rz = sprintf('%s%s',path_mat, 'G_rz.mat');
        load(filepath_G_rz);
        
        n_S = size(match_nit,1);
        
        % Set variables into aux. for function
        aux.info_site = info_site;
        aux.i = i;
        aux.G_R = G_R;
        aux.G_t = G_t;
        aux.G_rx = G_rx;
        aux.G_ry = G_ry;
        aux.G_rz = G_rz;
        aux.G_tx = G_tx;
        aux.G_ty = G_ty;
        aux.G_tz = G_tz;
        
        
        % Rename for RT there or back
        for RT = 1:2;
            if RT==1;
                match_R = match1_R;
                match_t = match1_t;
                n_plot = size(match2_R,1);
                matchret_R = cell(size(match2_R));
                matchret_t = cell(size(match2_t));
                for i= 1:size(match2_R,1);
                    for j = 1:size(match2_R,2);
                        if ~isempty(match2_R{i,j});
                            matchret_R{i,j} = match2_R{i,j}';
                            matchret_t{i,j} = -(match2_R{i,j}'*match2_t{i,j});
                        end
                    end
                end
            elseif RT ==2;
                for i= 1:size(match2_R,1);
                    for j = 1:size(match2_R,2);
                        if ~isempty(match2_R{i,j});
                            match_R{i,j} = match2_R{i,j}';
                            match_t{i,j} = -(match2_R{i,j}'*match2_t{i,j});
                        end
                    end
                end
                matchret_R = match1_R;
                matchret_t = match1_t;
            end
            
            match_rx =  180*cellfun(@decompose_rotation_rx,match_R)/pi;
            match_ry =  180*cellfun(@decompose_rotation_ry,match_R)/pi;
            match_rz =  180*cellfun(@decompose_rotation_rz,match_R)/pi;
            match_tx =  cellfun(@decompose_translation_tx,match_t);
            match_ty =  cellfun(@decompose_translation_ty,match_t);
            match_tz =  cellfun(@decompose_translation_tz,match_t);
            matchret_rx =  180*cellfun(@decompose_rotation_rx,matchret_R)/pi;
            matchret_ry =  180*cellfun(@decompose_rotation_ry,matchret_R)/pi;
            matchret_rz =  180*cellfun(@decompose_rotation_rz,matchret_R)/pi;
            matchret_tx =  cellfun(@decompose_translation_tx,matchret_t);
            matchret_ty =  cellfun(@decompose_translation_ty,matchret_t);
            matchret_tz =  cellfun(@decompose_translation_tz,matchret_t);
            diff_rx = abs(match_rx - G_rx);
            diff_ry = abs(match_ry - G_ry);
            diff_rz = abs(match_rz - G_rz);
            diff_tx = abs(match_tx - G_tx);
            diff_ty = abs(match_ty - G_ty);
            diff_tz = abs(match_tz - G_tz);
            is_guess = (diff_rx < t_r) & (diff_ry < t_r) & (diff_rz < t_r) & ...
                (diff_tx < t_t) & (diff_ty < t_t) & (diff_tz < t_t);
            
            aux.match_R = match_R;
            aux.match_t = match_t;
            aux.match_rx = match_rx;
            aux.match_ry = match_ry;
            aux.match_rz = match_rz;
            aux.match_tx = match_tx;
            aux.match_ty = match_ty;
            aux.match_tz = match_tz;
            aux.matchret_R = matchret_R;
            aux.matchret_t = matchret_t;
            aux.matchret_rx = matchret_rx;
            aux.matchret_ry = matchret_ry;
            aux.matchret_rz = matchret_rz;
            aux.matchret_tx = matchret_tx;
            aux.matchret_ty = matchret_ty;
            aux.matchret_tz = matchret_tz;
            aux.is_guess = is_guess;
            
            isempty0 = cellfun(@isempty, match_R);
            
            n_scan = 25;
            vec_plot = 1:25;
            vec_x = [10:-5:-10, -10:5:10, 10:-5:-10, -10:5:10, 10:-5:-10];
            vec_y = [repmat(-10, [1,5]),repmat(-5, [1,5]),zeros(1,5),repmat(5, [1,5]),repmat(10, [1,5])];
            
            prompt = 'Good (g) or bad (b)?';
            
            aux.n_scan = n_scan;
            aux.vec_plot = vec_plot;
            aux.vec_x = vec_x;
            aux.vec_y = vec_y;
            aux.RT = RT;
            aux.isempty0 = isempty0;
            
            
            for i= 1:n_S;
                path_matplot = sprintf('%s%02.0f%smat%s',path_site,i,info_slash,info_slash);
                filepath_match1_tf = sprintf('%s%s%d%s',path_matplot, 'match',RT,'_tf.mat');
                if exist(filepath_match1_tf, 'file');
                    %  continue
                end
                
                aux.i = i;
                % Plot guess in separate figure
                %hfigguess = label_plot_guess(aux);
                
                aux.match_RMSE = match12_RMSE;
                aux.match_tdiff = match_tdiff;
                aux.match_txydiff = match_txydiff;
                hfigrmse = label_plot_rmse(aux);
                
                
                result = input(prompt, 's');
                if strcmp(result, 'g');
                    isgoodguess = 1;
                elseif strcmp(result,'b');
                    isgoodguess = 0;
                else
                    error('invalid key press');
                end
                
                if ~isgoodguess;
                    
                    % For each i in a site
                    % Plot all #1
                    %err_alpha = match12_RMSE(i,:);
                    %err_color = (double(vec2cmap(err_alpha, 'jet', 0,2 )))./255;
                    %figure('position', [251         392        1156         559]);
                    figure('position', [68         413        1593         538]);
                    subplot(1,3,1)
                    if RT==1;
                        title(sprintf('Pairwise Forward Pose Estimates for Site %3.0f-%g', info_site,i));
                    elseif RT ==2;
                        title(sprintf('Pairwise Reverse Pose Estimates for Site %3.0f-%g', info_site,i));
                    end
                    hold on;
                    xlabel('x');
                    ylabel('y');
                    zlabel('z');
                    axis equal
                    % Plot all pairwise matches from j
                    for j = 1:n_S;
                        if ~isempty0(i,j);
                            x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
                            y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
                            z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
                            plot3(x_axist(1,:),x_axist(2,:),.1+x_axist(3,:),'-r', 'linewidth',2)
                            plot3(y_axist(1,:),y_axist(2,:),.1+y_axist(3,:),'-g', 'linewidth',2)
                            plot3(z_axist(1,:),z_axist(2,:),.1+z_axist(3,:),'-b', 'linewidth',2)
                            textloc = (x_axist(:,2) + y_axist(:,2))/2;
                            textstr = sprintf('%g', j);
                            text(textloc(1), textloc(2), textloc(3), textstr);
                            % h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
                            %     z_sph+match_t{i,j}(3));
                            % alpha(0.2)
                            % set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
                        end
                    end
                    curr1_xlim = xlim;
                    curr1_xlim(1) = curr1_xlim(1) -2;
                    curr1_xlim(2) = curr1_xlim(2) +2;
                    curr1_ylim = ylim;
                    curr1_ylim(1) = curr1_ylim(1) -2;
                    curr1_ylim(2) = curr1_ylim(2) +2;
                    xrange = curr1_xlim(2) - curr1_xlim(1);
                    yrange = curr1_ylim(2) - curr1_ylim(1);
                    if xrange>yrange;
                        diff = xrange-yrange;
                        curr1_ylim(1) = curr1_ylim(1) - diff/2;
                        curr1_ylim(2) = curr1_ylim(2) + diff/2;
                    elseif yrange>xrange;
                        diff = yrange-xrange;
                        curr1_xlim(1) = curr1_xlim(1) - diff/2;
                        curr1_xlim(2) = curr1_xlim(2) + diff/2;
                    end
                    
                    
                    
                    xlimMode = 'manual' ;
                    ylimMode = 'manual';
                    zlimMode = 'manual';
                    
                    % Plot reference map
                    subplot(1,3,2);
                    title('Reference locations');
                    hold on
                    for j = 1:n_S;
                        text(vec_x(j)-vec_x(i), vec_y(j)-vec_y(i),sprintf('%g', vec_plot(j)),...
                            'HorizontalAlignment', 'center', 'fontsize', 14);
                    end
                    curr2_xlim(1) = min(vec_x-vec_x(i)) -2;
                    curr2_xlim(2) = max(vec_x-vec_x(i)) +2;
                    curr2_ylim(1) = min(vec_y-vec_y(i)) -2;
                    curr2_ylim(2) = max(vec_y-vec_y(i)) +2;
                    axis([curr2_xlim curr2_ylim]);
                    
                    
                    % Plot graph best guess
                    subplot(1,3,3);
                    title('Graph guess');
                    hold on;
                    xlabel('x');
                    ylabel('y');
                    zlabel('z');
                    axis equal
                    % Plot all pairwise matches from j
                    for j = 1:n_S;
                        if ~isempty0(i,j);
                            x_axist = G_R{i,j}*x_axis + repmat(G_t{i,j},[1,2]);
                            y_axist = G_R{i,j}*y_axis + repmat(G_t{i,j},[1,2]);
                            z_axist = G_R{i,j}*z_axis + repmat(G_t{i,j},[1,2]);
                            plot3(x_axist(1,:),x_axist(2,:),.1+x_axist(3,:),'-r', 'linewidth',2)
                            plot3(y_axist(1,:),y_axist(2,:),.1+y_axist(3,:),'-g', 'linewidth',2)
                            plot3(z_axist(1,:),z_axist(2,:),.1+z_axist(3,:),'-b', 'linewidth',2)
                            textloc = (x_axist(:,2) + y_axist(:,2))/2;
                            textstr = sprintf('%g', j);
                            text(textloc(1), textloc(2), textloc(3), textstr);
                            % h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
                            %     z_sph+match_t{i,j}(3));
                            % alpha(0.2)
                            % set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
                        end
                    end
                    curr3_xlim = xlim;
                    curr3_xlim(1) = curr3_xlim(1) -2;
                    curr3_xlim(2) = curr3_xlim(2) +2;
                    curr3_ylim = ylim;
                    curr3_ylim(1) = curr3_ylim(1) -2;
                    curr3_ylim(2) = curr3_ylim(2) +2;
                    xrange = curr3_xlim(2) - curr3_xlim(1);
                    yrange = curr3_ylim(2) - curr3_ylim(1);
                    if xrange>yrange;
                        diff = xrange-yrange;
                        curr3_ylim(1) = curr3_ylim(1) - diff/2;
                        curr3_ylim(2) = curr3_ylim(2) + diff/2;
                    elseif yrange>xrange;
                        diff = yrange-xrange;
                        curr3_xlim(1) = curr3_xlim(1) - diff/2;
                        curr3_xlim(2) = curr3_xlim(2) + diff/2;
                    end
                    
                    curr_xlim = [min([curr1_xlim(1),curr2_xlim(1),curr3_xlim(1)]) max([curr1_xlim(2),curr2_xlim(2),curr3_xlim(2)])];
                    curr_ylim = [min([curr1_ylim(1),curr2_ylim(1),curr3_ylim(1)]) max([curr1_ylim(2),curr2_ylim(2),curr3_ylim(2)])];
                    subplot(1,3,1);
                    axis([curr_xlim curr_ylim]);
                    subplot(1,3,2);
                    axis([curr_xlim curr_ylim]);
                    subplot(1,3,3);
                    axis([curr_xlim curr_ylim]);
                    
                    % Request user input and label Good/Bad
                    match1_tf = nan(n_S,1);
                    for j = 1:n_S;
                        if isempty(match_t{i,j});
                            match1_tf(j) = nan;
                            continue
                        end
                        subplot(1,3,1);
                        axis([curr_xlim curr_ylim]);
                        h0 = circle(match_t{i,j}(1:2), 1, 50, '-k');
                        h1 = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                        h2 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,'g');
                        set(h1, 'FaceAlpha', 0.4);
                        set(h2, 'FaceAlpha', 0.4);
                        
                        subplot(1,3,2);
                        axis([curr_xlim curr_ylim]);
                        h7 = circle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)], 1, 50, '-k');
                        h3 = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                        h4 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,'g');
                        set(h3, 'FaceAlpha', 0.4);
                        set(h4, 'FaceAlpha', 0.4);
                        
                        subplot(1,3,3);
                        axis([curr_xlim curr_ylim]);
                        h3a = circle(G_t{i,j}(1:2), 1, 50, '-k');
                        h3b = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                        h3c = filledCircle(G_t{i,j}(1:2),1,1000,'g');
                        set(h3b, 'FaceAlpha', 0.4);
                        set(h3c, 'FaceAlpha', 0.4);
                        
                        
                        result = input(prompt, 's');
                        if strcmp(result, 'g');
                            tf = 1;
                            clr = 'b';
                        elseif strcmp(result,'b');
                            tf = 0;
                            clr = 'r';
                        elseif strcmp(result, 'z');
                            jz = j - 1;
                            % Redo last plot
                            delete(h1);
                            delete(h2);
                            delete(h3);
                            delete(h4);
                            delete(h5);
                            delete(h6);
                            delete(h0);
                            delete(h7);
                            delete(h3a);
                            delete(h3b);
                            delete(h3c);
                            delete(h3d);
                            
                            subplot(1,3,1);
                            axis([curr_xlim curr_ylim]);
                            circle(match_t{i,jz}(1:2), 1, 50, '-k');
                            h1 = filledCircle(match_t{i,jz}(1:2),1,1000,'y');
                            h2 = filledCircle([vec_x(jz)-vec_x(i) vec_y(jz)-vec_y(i)],1,1000,'g');
                            set(h1, 'FaceAlpha', 0.4);
                            set(h2, 'FaceAlpha', 0.4);
                            
                            subplot(1,3,2);
                            axis([curr_xlim curr_ylim]);
                            circle([vec_x(jz)-vec_x(i) vec_y(jz)-vec_y(i)], 1, 50, '-k');
                            h3 = filledCircle(match_t{i,jz}(1:2),1,1000,'y');
                            h4 = filledCircle([vec_x(jz)-vec_x(i) vec_y(jz)-vec_y(i)],1,1000,'g');
                            set(h3, 'FaceAlpha', 0.4);
                            set(h4, 'FaceAlpha', 0.4);
                            
                            subplot(1,3,3);
                            axis([curr_xlim curr_ylim]);
                            circle(G_t{i,j}(1:2), 1, 50, '-k');
                            h3b = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                            h3c = filledCircle(G_t{i,j}(1:2),1,1000,'g');
                            set(h3b, 'FaceAlpha', 0.4);
                            set(h3c, 'FaceAlpha', 0.4);
                            
                            
                            
                            result = input(prompt, 's');
                            if strcmp(result, 'g');
                                tf = 1;
                                clr = 'b';
                            elseif strcmp(result,'b');
                                tf = 0;
                                clr = 'r';
                            end
                            delete(h1);
                            delete(h2);
                            delete(h3);
                            delete(h4);
                            delete(h3b);
                            delete(h3c);
                            
                            match1_tf(jz) = tf;
                            subplot(1,3,1);
                            h = filledCircle(match_t{i,jz}(1:2),1,1000,clr);
                            set(h, 'FaceAlpha', 0.2)
                            subplot(1,3,2);
                            h = filledCircle([vec_x(jz)-vec_x(i) vec_y(jz)-vec_y(i)],1,1000,clr);
                            set(h, 'FaceAlpha', 0.2)
                            subplot(1,3,3);
                            h = filledCircle(G_t{i,j}(1:2),1,1000,clr);
                            set(h, 'FaceAlpha', 0.4);
                            
                            
                            % Repeat for current
                            subplot(1,3,1);
                            axis([curr_xlim curr_ylim]);
                            circle(match_t{i,j}(1:2), 1, 50, '-k');
                            h1 = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                            h2 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,'g');
                            set(h1, 'FaceAlpha', 0.4);
                            set(h2, 'FaceAlpha', 0.4);
                            
                            subplot(1,3,2);
                            axis([curr_xlim curr_ylim]);
                            circle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)], 1, 50, '-k');
                            h3 = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                            h4 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,'g');
                            set(h3, 'FaceAlpha', 0.4);
                            set(h4, 'FaceAlpha', 0.4);
                            
                            subplot(1,3,3);
                            axis([curr_xlim curr_ylim]);
                            circle(G_t{i,j}(1:2), 1, 50, '-k');
                            h3b = filledCircle(match_t{i,j}(1:2),1,1000,'y');
                            h3c = filledCircle(G_t{i,j}(1:2),1,1000,'g');
                            set(h3b, 'FaceAlpha', 0.4);
                            set(h3c, 'FaceAlpha', 0.4);
                            
                            
                            result = input(prompt, 's');
                            if strcmp(result, 'g');
                                tf = 1;
                                clr = 'b';
                            elseif strcmp(result,'b');
                                tf = 0;
                                clr = 'r';
                            end
                        else
                            error('invalid key press');
                        end
                        delete(h1);
                        delete(h2);
                        delete(h3);
                        delete(h4);
                        delete(h3b);
                        delete(h3c);
                        
                        match1_tf(j) = tf;
                        subplot(1,3,1);
                        h5 = filledCircle(match_t{i,j}(1:2),1,1000,clr);
                        set(h5, 'FaceAlpha', 0.2)
                        subplot(1,3,2);
                        h6 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,clr);
                        set(h6, 'FaceAlpha', 0.2)
                        subplot(1,3,3);
                        h3d = filledCircle(G_t{i,j}(1:2),1,1000,clr);
                        set(h3d, 'FaceAlpha', 0.2);
                        
                        
                    end
                    delete(hfigguess);
                    
                else
                    match1_tf = is_guess(i,:);
                end
                
                % Save results
                filepath_match1_fig = sprintf('%s%s%d_%03.0f_%02.0f.png', path_figures,'truth',RT,...
                    info_site, i);
                export_fig(filepath_match1_fig, '-m1');
                path_matplot = sprintf('%s%02.0f%smat%s',path_site,i,info_slash,info_slash);
                filepath_match1_tf = sprintf('%s%s%d%s',path_matplot, 'match',RT,'_tf.mat');
                if RT==1;
                    save(filepath_match1_tf, 'match1_tf');
                elseif RT==2;
                    match2_tf = match1_tf;
                    save(filepath_match1_tf, 'match2_tf');
                end
                close all
                %}
            end
            %{
        % Plot all for #2
        figure('position', [251         392        1156         559]);
        subplot(1,2,1)
        title(sprintf('Pairwise Reverse Pose Estimates for Site %3.0f-%g', info_site,i));
        hold on;
        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis equal
        for j = 1:n_S;
            if ~isempty2(i,j);
                Rflip = match2_R{i,j}';
                tflip = -(match2_R{i,j}')*match2_t{i,j};
                x_axis1t = Rflip*x_axis + repmat(tflip,[1,2]);
                y_axis1t = Rflip*y_axis + repmat(tflip,[1,2]);
                z_axis1t = Rflip*z_axis + repmat(tflip,[1,2]);
                plot3(x_axis1t(1,:),x_axis1t(2,:),x_axis1t(3,:),'-r', 'linewidth',2)
                plot3(y_axis1t(1,:),y_axis1t(2,:),y_axis1t(3,:),'-g', 'linewidth',2)
                plot3(z_axis1t(1,:),z_axis1t(2,:),z_axis1t(3,:),'-b', 'linewidth',2)
                textloc = (x_axis1t(:,2) + y_axis1t(:,2))/2;
                textstr = sprintf('%g-o', j);
                text(textloc(1), textloc(2), textloc(3), textstr);
            end
            % end
        end
        curr1_xlim = xlim;
        curr1_xlim(1) = curr1_xlim(1) -2;
        curr1_xlim(2) = curr1_xlim(2) +2;
        curr1_ylim = ylim;
        curr1_ylim(1) = curr1_ylim(1) -2;
        curr1_ylim(2) = curr1_ylim(2) +2;
        
        % Plot reference
        subplot(1,2,2);
        title('Reference locations');
        hold on
        for j = 1:n_S;
            text(vec_x(j)-vec_x(i), vec_y(j)-vec_y(i),sprintf('%g', vec_plot(j)),...
                'HorizontalAlignment', 'center','fontsize', 14);
        end
        
        curr2_xlim(1) = min(vec_x-vec_x(i)) -2;
        curr2_xlim(2) = max(vec_x-vec_x(i)) +2;
        curr2_ylim(1) = min(vec_y-vec_y(i)) -2;
        curr2_ylim(2) = max(vec_y-vec_y(i)) +2;
        axis([curr2_xlim curr2_ylim]);
        
        curr_xlim = [min(curr1_xlim(1),curr2_xlim(1)) max(curr1_xlim(2),curr2_xlim(2))];
        curr_ylim = [min(curr1_ylim(1),curr2_ylim(1)) max(curr1_ylim(2),curr2_ylim(2))];
        subplot(1,2,1);
        axis([curr_xlim curr_ylim]);
        subplot(1,2,2);
        axis([curr_xlim curr_ylim]);
        
        
        match2_tf = false(n_S,1);
        % Go through each and request user input
        for j = 1:n_S;
            subplot(1,2,1);
            axis([curr_xlim curr_ylim]);
            circle(match_t{i,j}(1:2), 1, 50, '-k');
            h1 = filledCircle(match_t{i,j}(1:2),1,1000,'y');
            h2 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,'g');
            set(h1, 'FaceAlpha', 0.4);
            set(h2, 'FaceAlpha', 0.4);
            subplot(1,2,2);
            axis([curr_xlim curr_ylim]);
            circle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)], 1, 50, '-k');
            h3 = filledCircle(match_t{i,j}(1:2),1,1000,'y');
            h4 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,'g');
            set(h3, 'FaceAlpha', 0.4);
            set(h4, 'FaceAlpha', 0.4);
            result = input(prompt, 's');
            if strcmp(result, 'g');
                tf = true;
                clr = 'b';
            elseif strcmp(result,'b');
                tf = false;
                clr = 'r';
            else
                error('invalid key press');
            end
            delete(h1);
            delete(h2);
            delete(h3);
            delete(h4);
            match2_tf(j) = tf;
            subplot(1,2,1);
            h = filledCircle(tflip(1:2),1,1000,clr);
            set(h, 'FaceAlpha', 0.2)
            subplot(1,2,2);
            h = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,clr);
            set(h, 'FaceAlpha', 0.2)
        end
        set(gcf, 'color', 'white')
        
        % save #2
        filepath_match2_fig = sprintf('%s%s%03.0f_%02.0f.png', path_figures,'truth2_',...
            info_site, i);
        export_fig(filepath_match2_fig, '-m1');
        path_matplot = sprintf('%s%02.0f%smat%s',path_site,i,info_slash,info_slash);
        filepath_match2_tf = sprintf('%s%s',path_matplot, 'match2_tf.mat');
        save(filepath_match2_tf, 'match1_tf');
        close all
            %}
        end
end





foo  =1;
