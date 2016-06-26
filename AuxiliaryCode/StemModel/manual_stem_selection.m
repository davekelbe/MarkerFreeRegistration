function [  ] = manual_stem_selection(path_up, info_site, info_plot, info_exp, info_suffix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

info_slash = '/';
path_top = sprintf('%s%s%s%s%s%03.0f%s%03.0f%s',path_up,...
    info_exp, info_slash, info_suffix,info_slash,info_site, info_slash,info_plot,info_slash);
path_mat = sprintf('%s%s%s',path_top,'mat',info_slash);
path_png = sprintf('%s%s%s',path_top,'png',info_slash);

filepath_seg_col = sprintf('%s%s%s',path_mat,'seg_col','.mat');
filepath_seg_row = sprintf('%s%s%s',path_mat,'seg_row','.mat');
filepath_seg_iter = sprintf('%s%s%s',path_mat,'seg_iter','.mat');
load(filepath_seg_col);
load(filepath_seg_row);
load(filepath_seg_iter);

filepath_I12ieq = sprintf('%s%s_%03.0f-%03.0f%s',path_png,...
    'I12ieq',info_site,info_plot,'.png');
I12ieq = imread(filepath_I12ieq);
filepath_I12ieq_outline = sprintf('%sI12ieq_outline_%03.0f-%03.0f.png',path_png,info_site,info_plot);


filepath_I12x = sprintf('%s%s%s',path_mat,'I12x','.mat');
filepath_I12y = sprintf('%s%s%s',path_mat,'I12y','.mat');
filepath_I12z = sprintf('%s%s%s',path_mat,'I12z','.mat');
filepath_I12r = sprintf('%s%s%s',path_mat,'I12r','.mat');

load(filepath_I12x);
load(filepath_I12y);
load(filepath_I12z);
load(filepath_I12r);


[n_row, n_col] = size(I12ieq);
fh1 = figure;
imagesc(I12ieq);
hold on
for s = 1:numel(seg_col);
    plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-r')
end
set(gca,'Units','normalized','Position',[0 0 1 1]);
[n_row, n_col,~] = size(I12ieq);
set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
is_remove = false(numel(seg_col),1);
save('~/temp_is_remove', 'is_remove');

% Load seg data
filepath_seg_row = sprintf('%s%s%s',path_mat,'seg_row','.mat');
filepath_seg_col = sprintf('%s%s%s',path_mat,'seg_col','.mat');
filepath_seg_iter = sprintf('%s%s%s',path_mat,'seg_iter','.mat');
filepath_seg_z = sprintf('%s%s%s',path_mat,'seg_z','.mat');
filepath_seg_y = sprintf('%s%s%s',path_mat,'seg_y','.mat');
filepath_seg_r = sprintf('%s%s%s',path_mat,'seg_r','.mat');
filepath_seg_index = sprintf('%s%s%s',path_mat,'seg_index','.mat');
filepath_seg_fill = sprintf('%s%s%s',path_mat,'seg_fill','.mat');
filepath_seg_taper = sprintf('%s%s%s',path_mat,'seg_taper','.mat');
filepath_seg_lean = sprintf('%s%s%s',path_mat,'seg_lean','.mat');
filepath_seg_anorm = sprintf('%s%s%s',path_mat,'seg_anorm','.mat');

load(filepath_seg_row);
load(filepath_seg_col);
load(filepath_seg_iter);
load(filepath_seg_z);
load(filepath_seg_y);
load(filepath_seg_r);
load(filepath_seg_index);
load(filepath_seg_fill);
load(filepath_seg_taper);
load(filepath_seg_lean);
load(filepath_seg_anorm);

% Set thresholds
t_r_adj = 2;
t_theta = 25;
t_r1_max = 1;
t_r1_min = 0;
t_taper_min = -5;
t_taper_max = 5;
[info_sampling,~,~,~] = SICK_specs();
t_min_density = 0.1;
%t_min_obj_size = 0; % forced
%t_min_pts = 0; % forced
t_coverage = 0; % forced
sub_axes = 0; % forced
options_verbose_model_cyl = false;
a1norm =[0 0 1]'; % forced
filepath_axis_a = sprintf('%s%s%s',path_mat,'axis_a','.mat');
filepath_axis_e = sprintf('%s%s%s',path_mat,'axis_e','.mat');
load(filepath_axis_a);
load(filepath_axis_e);

%I12ieq

continueflag = true;
while continueflag;
    
    % choice = questdlg('Please zoom if you like', ...
    % 'Hi Jan', ...
    % 'OK', 'OK');
    
    % Construct a questdlg with three options
    choice = questdlg('Would you like to remove, add more segments, or finish?', ...
        '', ...
        'Remove','Digitize New','Finished', 'Remove');
    seg_iterT = seg_iter;
    
    switch choice
        case 'Remove';
            figure(fh1);
                 
            h.mybrush = brush;
            set(h.mybrush, 'Enable', 'on', 'ActionPostCallback', @displayBrushData); 

            k = 0;
            while ~k;
                 k = waitforbuttonpress;
            end

            load('~/temp_is_remove');
            in = is_remove;
            %{
            % User selects within points
            %[colq, rowq] = getpts;

            % Check which segments were selected
           
            n_seg = numel(seg_col);
            in = false(n_seg,1);
            for s = 1:n_seg;
                ptsin = inpolygon(rowq, colq, seg_row{s}, seg_col{s});
                in(s) = any(ptsin);
            end
            %}
            seg_row = seg_row(~in);
            seg_col = seg_col(~in);
            seg_iter = seg_iter(~in);
           % seg_a = seg_a(~in,:);
            seg_anorm = seg_anorm(~in,:);
            seg_fill = seg_fill(~in);
            seg_index = seg_index(~in);
            seg_lean = seg_lean(~in);
            seg_r = seg_r(~in,:);
            seg_taper = seg_taper(~in);
            seg_y = seg_y(~in,:);
            seg_z = seg_z(~in,:);

            n_seg = numel(seg_row);
            %{
figure;
imagesc(I12ieq);
hold on
for s = 1:numel(seg_colT);
    plot([seg_colT{s} seg_colT{s}(1)],[seg_rowT{s} seg_rowT{s}(1)],'-r')
end
            %}
            
            %{
            % Mark in image
            ind = 1:n_col*n_row;
            [pixrow, pixcol] = ind2sub([n_row, n_col], ind);
            n_pts = n_row*n_col;
            ptsin = false(1, n_pts);
            for s = 1:n_segT;
                ptscurr = inpolygon(pixrow, pixcol, seg_rowT{s}, seg_colT{s});
                ptsin = any([ptsin; ptscurr]);
            end
            pixrowT = pixrow(ptsin);
            pixcolT = pixcol(ptsin);
            
            hold on
            scatter(pixcolT, pixrowT, '.r')
            %}
            
            delete(fh1);
            fh1 = figure;
            imagesc(I12ieq);
            hold on
            for s = 1:numel(seg_col);
                plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-r')
            end
            set(gca,'Units','normalized','Position',[0 0 1 1]);
            [n_row, n_col,~] = size(I12ieq);
            set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
            foo = 1;
            
            save(filepath_seg_row, 'seg_row');
            save(filepath_seg_col, 'seg_col');
            save(filepath_seg_iter, 'seg_iter');
            save(filepath_seg_z, 'seg_z');
            save(filepath_seg_y, 'seg_y');
            save(filepath_seg_r, 'seg_r');
            save(filepath_seg_index, 'seg_index');
            save(filepath_seg_fill, 'seg_fill');
            save(filepath_seg_taper, 'seg_taper');
            save(filepath_seg_lean, 'seg_lean');
            save(filepath_seg_anorm, 'seg_anorm');
           % save(filepath_seg_a, 'seg_a');

        case 'Digitize New';
            %BW = roipoly;
            
            % Make mask
            h = impoly;
            mask = createMask(h);
            [coltrunk,rowtrunk] = ginput(1);
            coltrunk = round(coltrunk);
            rowtrunk = round(rowtrunk);
            xtrunk = I12x(rowtrunk, coltrunk);
            ytrunk = I12y(rowtrunk, coltrunk);
            ztrunk = I12z(rowtrunk, coltrunk);
            rtrunk = sqrt(xtrunk.^2 + ytrunk.^2 + ztrunk.^2);
            
         %   temp = I12ieq(rowtrunk-70:rowtrunk+70,coltrunk-70:coltrunk+70,:);
         %   figure; imagesc(temp);

%             pos = getPosition(h);
%             seg_row_temp = pos(:,2);
%             seg_col_temp = pos(:,1);
%             [~,row_sort_ix] = sort(seg_row_temp,'ascend');
%             [~,col_sort_ix] = sort(seg_col_temp,'ascend');
%             
%             ism = ismember(row_sort_ix(1:2), col_sort_ix(3:4));
%             temp = row_sort_ix(1:2);
%             UR = temp(ism);
%             ism = ismember(row_sort_ix(3:4), col_sort_ix(3:4));
%             temp = row_sort_ix(3:4);
%             LR = temp(ism);
%             ism = ismember(row_sort_ix(3:4), col_sort_ix(1:2));
%             temp = row_sort_ix(3:4);
%             LL = temp(ism);
%             ism = ismember(row_sort_ix(1:2), col_sort_ix(1:2));
%             temp = row_sort_ix(1:2);
%             UL = temp(ism);
%             
%             n_seg = numel(seg_row);
%             seg_row{n_seg+1} = [seg_row_temp(UR) seg_row_temp(LR) ...
%                 seg_row_temp(LL) seg_row_temp(UL)];
%             seg_col{n_seg+1} = [seg_col_temp(UR) seg_col_temp(LR) ...
%                 seg_col_temp(LL) seg_col_temp(UL)];
            %n_seg = n_seg + 1;
            %plot(seg_col{n_seg}(4),seg_row{n_seg}(4),'*c', 'markersize', 60)
            
            
            ptsx = I12x(mask);
            ptsy = I12y(mask);
            ptsz = I12z(mask);
            ptsr = sqrt(ptsx.^2 + ptsy.^2 + ptsz.^2);
            isclose = (ptsr<=rtrunk + 0.5);
            ptsx = ptsx(isclose);
            ptsy = ptsy(isclose);
            ptsz = ptsz(isclose);
            ptsr = ptsr(isclose);
                        
            %seg_a = seg_z - seg_y;
            %seg_anorm = seg_a./repmat(sqrt(sum(seg_a.*seg_a,2)),[1,3]);
            %seg_lean = acosd(seg_anorm(:,3));
            
            %ptsr = sqrt(ptsx.^2 + ptsy.^2 + ptsz.^2);
            mean_range = median(ptsr);
            [ ~, ~, t_min_obj_size, t_min_pts ] = SICK_specs( mean_range, t_r_adj);
            
            % Fit cylinder to selected points
            [cyl_ctr,rcirc,is_inlier,E,~,~,~,cyl_left,cyl_right,return_str,coverage,t_inlier, status] = model_cyl_allcomp_manual(ptsx,ptsy,ptsz,...
                t_r_adj,t_theta,t_r1_max,t_r1_min, t_taper_min,...
                t_taper_max,info_sampling,t_min_density,t_min_obj_size,...
                t_min_pts, t_coverage, sub_axes, options_verbose_model_cyl,a1norm,...
                axis_a, axis_e, I12ieq); % was 0.5*t_min_density
            r1 = rcirc; % [bottom radius, top radius]
            z1 = cyl_ctr(:,end); % top
            y1 = cyl_ctr(:,1); % bottom
            a1 = z1 - y1;
            a1norm = a1/norm(a1);
            
            seg_r(end+1,:) = r1';
            seg_y(end+1,:) = y1';
            seg_z(end+1,:) = z1';
         %   seg_a(end+1,:) = a1';
            seg_anorm(end+1,:) = a1norm';
            seg_fill(end+1) = 1;
            seg_index{end+1} = []; % force
            seg_iter(end+1) = 0; % force
            seg_iterT(end+1) = 0; % force
            seg_taper(end+1) =  atand((r1(1)-r1(2))/norm(z1-y1)); % Expect small positive values
            seg_lean(end+1) = 0; % force
            
            clear seg_col_temp
            clear seg_row_temp
            
            % update seg_col and seg_row
                    z_left = cyl_left(:,end);
                    y_left = cyl_left(:,1);
                    z_right = cyl_right(:,end);
                    y_right = cyl_right(:,1);
                    
            [col, row, ~] = andrieu_fix_wraparound(...
                I12r,axis_a,axis_e,seg_iter,...
                z_left, y_left, z_right, y_right);
            seg_row(end+1) = row;
            seg_col(end+1) = col;
%             is_finite = false(numel(seg_col),1);
%             for s = 1:numel(seg_col);
%                 is_finite(s) = ~(any(~isfinite(seg_col{s})) || any(~isfinite(seg_row{s})));
%             end
%             seg_iter = seg_iter(is_finite);
%             seg_row = seg_row(is_finite);
%             seg_col = seg_col(is_finite);
%             save(filepath_seg_col,'seg_col');
%             save(filepath_seg_row,'seg_row');
%             save(filepath_seg_iter,'seg_iter');
            
            
            %foo = 1
            delete(h);
            n_seg = numel(seg_fill);
            plot([seg_col{n_seg} seg_col{n_seg}(1)],[seg_row{n_seg} seg_row{n_seg}(1)],'-b', 'linewidth',2)
            
            save(filepath_seg_row, 'seg_row');
            save(filepath_seg_col, 'seg_col');
            save(filepath_seg_iter, 'seg_iter');
            save(filepath_seg_z, 'seg_z');
            save(filepath_seg_y, 'seg_y');
            save(filepath_seg_r, 'seg_r');
            save(filepath_seg_index, 'seg_index');
            save(filepath_seg_fill, 'seg_fill');
            save(filepath_seg_taper, 'seg_taper');
            save(filepath_seg_lean, 'seg_lean');
            save(filepath_seg_anorm, 'seg_anorm');
        %    save(filepath_seg_a, 'seg_a');

            %             figure;
            %             scatter3(ptsx,ptsy, ptsz);
            %             [ seg_z, seg_y, seg_r, seg_lean, seg_anorm, seg_taper ] = model_manual_subset( ptsx, ptsy, ptsz );
            
            % Subset points
            
            
            %pos = getPosition(h);
            %colmin = min(pos(:,2));
            %colmax = max(pos(:,2));
            
            % Make segment from ROI
            foo = 1;
            
            
        case 'Finished'
            continueflag = false;
            delete(fh1);

            %{
%             filepath_seg_row = sprintf('%s%s%s',path_mat,'seg_row','.mat');
%             filepath_seg_col = sprintf('%s%s%s',path_mat,'seg_col','.mat');
%             filepath_seg_iter = sprintf('%s%s%s',path_mat,'seg_iter','.mat');
%             filepath_seg_z = sprintf('%s%s%s',path_mat,'seg_z','.mat');
%             filepath_seg_y = sprintf('%s%s%s',path_mat,'seg_y','.mat');
%             filepath_seg_r = sprintf('%s%s%s',path_mat,'seg_r','.mat');
%             filepath_seg_index = sprintf('%s%s%s',path_mat,'seg_index','.mat');
%             filepath_seg_fill = sprintf('%s%s%s',path_mat,'seg_fill','.mat');
%             filepath_seg_taper = sprintf('%s%s%s',path_mat,'seg_taper','.mat');
%             filepath_seg_lean = sprintf('%s%s%s',path_mat,'seg_lean','.mat');
%             filepath_seg_anorm = sprintf('%s%s%s',path_mat,'seg_anorm','.mat');
%             
%             load(filepath_seg_row);
%             load(filepath_seg_col);
%             load(filepath_seg_iter);
%             load(filepath_seg_z);
%             load(filepath_seg_y);
%             load(filepath_seg_r);
%             load(filepath_seg_index);
%             load(filepath_seg_fill);
%             load(filepath_seg_taper);
%             load(filepath_seg_lean);
%             load(filepath_seg_anorm);
%             
%             % Subset by manually identified segments
%             seg_isvalid = ismember(seg_iter, seg_iterT, 'rows');
%             
%             seg_row = seg_row(seg_isvalid);
%             seg_col = seg_col(seg_isvalid);
%             seg_iter = seg_iter(seg_isvalid);
%             seg_z = seg_z(seg_isvalid,:);
%             seg_y = seg_y(seg_isvalid,:);
%             seg_r = seg_r(seg_isvalid,:);
%             seg_index = seg_index(seg_isvalid);
%             seg_fill = seg_fill(seg_isvalid);
%             seg_taper = seg_taper(seg_isvalid);
%             seg_lean = seg_lean(seg_isvalid);
%             seg_anorm = seg_anorm(seg_isvalid,:);
%             
%             save(filepath_seg_row, 'seg_row');
%             save(filepath_seg_col, 'seg_col');
%             save(filepath_seg_iter, 'seg_iter');
%             save(filepath_seg_z, 'seg_z');
%             save(filepath_seg_y, 'seg_y');
%             save(filepath_seg_r, 'seg_r');
%             save(filepath_seg_index, 'seg_index');
%             save(filepath_seg_fill, 'seg_fill');
%             save(filepath_seg_taper, 'seg_taper');
%             save(filepath_seg_lean, 'seg_lean');
%             save(filepath_seg_anorm, 'seg_anorm');
%}
    end
end



%    f = getframe(gcf);
%    imwrite(f.cdata,filepath_I12ieq_outline,'png');




function displayBrushData(~, eventdata)
nlines = length(eventdata.Axes.Children)-1;
brushdata = cell(nlines, 1);
%global is_remove
load('~/temp_is_remove');
is_remove_now = false(nlines,1);
for ii = 1:nlines;
    brushdata{ii} = eventdata.Axes.Children(ii).BrushHandles.Children(1).VertexData;
    if size(brushdata{ii},2)~=0;
        is_remove_now(ii) = true;
    end
    fprintf('Line %i\n', ii)
    fprintf('X: %f Y: %f Z: %f\n', brushdata{ii})
end
is_remove_now = flipud(is_remove_now);
is_remove = is_remove | is_remove_now;
save('~/temp_is_remove', 'is_remove');

