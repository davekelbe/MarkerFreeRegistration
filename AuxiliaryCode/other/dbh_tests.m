% DBH Tests
% 
% Code written for IGARSS 2014 Paper
% Also may be used for validation of Paper01
% 
% (c) David Kelbe, Rochester Institute of Technology 
% 

path_source = 'D:\Users\djk2312\Documents\DBH tests\';

set(0,'DefaultFigureRenderer','OpenGL')
set(0,'defaultfigureposition',[966  322  557 472]);
slash = '\';
warning('off','MATLAB:nearlySingularMatrix');


% Filename LUT
filename_lut = 'D:\Users\djk2312\Documents\DBH tests\Lidar_Filenames_DBHTests.txt';
fid = fopen(filename_lut);
C = textscan(fid,'%u %s %s %s %s %s %s %s %s','headerlines',1,'delimiter','\t','collectOutput',true);
fclose all;
range = C{1};
filename_matrix = C{2};
dbh = [4 8.05 12 16 20 24 28.05 32];
n_iter = 200;
n_range = numel(range);
n_dbh = numel(dbh);
alpha = zeros(n_iter,n_range, n_dbh);
kasa = zeros(n_iter,n_range, n_dbh);
pratt = zeros(n_iter,n_range, n_dbh);

% Cropping info
y_ctr = -[ 888 2805 4700 6600 8510 10393 12250 14154 16050 17950 19800 21636 23570]/1000;
x_ctr = [  351 1010 1700 2300 2940 3551  4200  4859  5500  6120  6736  7372  7973]/1000;
z_lo = [   -63 -69  -78  -81  -85  -90   -100  -100  -110  -110  -120  -130  -140]/100;
z_hi = [   50  50   50   50   50   50    50    50    50    50    50    50    50]/100;
for d = 8:n_dbh;
    
        fprintf('\n**************\n');
        fprintf('dbh   = %3.0f \n', dbh(d));
        
    for r = 4:n_range;
        
        %fprintf('\n**************\n');
        %fprintf('range = %3.0f \n', range(r));
        %fprintf('dbh   = %3.0f \n', dbh(d));
        
        % Set up initial info
        info_site = range(r);
        info_plot = dbh(d);
        info_exp = 'DBHTest';
        info_suffix = '01-06';
        info_name = sprintf('%03.0f_%02.0f_%s',info_site,info_plot,info_suffix);
        
        % Make directories
        path_common = sprintf('%s%s%s%s%s','D:\Users\djk2312\Documents\',...
            info_exp, slash, 'Common',slash);
        path_top = sprintf('%s%s%s%s%s%03.0f%s%02.0f%s','D:\Users\djk2312\Documents\',...
            info_exp, slash, info_suffix,slash,info_site,slash,info_plot,slash);
        path_mat = sprintf('%s%s%s',path_top,'mat',slash);
        path_fig = sprintf('%s%s%s',path_top,'fig',slash);
        path_eps = sprintf('%s%s%s',path_top,'eps',slash);
        path_ply = sprintf('%s%s%s',path_top,'ply',slash);
        path_png = sprintf('%s%s%s',path_top,'png',slash);
        % Create directories
        if ~exist(path_common,'dir');
            mkdir(path_common);
        end
        if ~exist(path_top,'dir');
            mkdir(path_top);
        end
        if ~exist(path_mat,'dir');
            mkdir(path_mat);
        end
        if ~exist(path_fig,'dir');
            mkdir(path_fig);
        end
        if ~exist(path_eps,'dir');
            mkdir(path_eps);
        end
        if ~exist(path_ply,'dir');
            mkdir(path_ply);
        end
        if ~exist(path_png,'dir');
            mkdir(path_png);
        end
        
        % Load data
        filepath_data_i  = sprintf('%s%s',path_mat,'data_i.mat');
        filepath_data_ieq  = sprintf('%s%s',path_mat,'data_ieq.mat');
        filepath_data_index  = sprintf('%s%s',path_mat,'data_index.mat');
        filepath_data_x  = sprintf('%s%s',path_mat,'data_x.mat');
        filepath_data_y  = sprintf('%s%s',path_mat,'data_y.mat');
        filepath_data_z  = sprintf('%s%s',path_mat,'data_z.mat');
        filepath_data_n  = sprintf('%s%s',path_mat,'data_n.mat');
        filepath_data_xy = sprintf('%s%s',path_mat,'data_xy.mat');
        filepath_data_r  = sprintf('%s%s',path_mat,'data_r.mat');
        filepath_data_a  = sprintf('%s%s',path_mat,'data_a.mat');
        filepath_data_e  = sprintf('%s%s',path_mat,'data_e.mat');
        filepath_logical_sub  = sprintf('%s%s',path_mat,'logical_sub.mat');
        % Load data if exist
        if exist(filepath_data_x,'file');
            load(filepath_data_i);
            load(filepath_data_x);
            load(filepath_data_y);
            load(filepath_data_z);
            load(filepath_data_n);
            load(filepath_data_xy);
            load(filepath_data_r);
            load(filepath_data_a);
            load(filepath_data_e);
            % Else compute
        elseif ~exist(filepath_data_x, 'file');
            filename_source = char(filename_matrix(r,d));
            filepath_source = sprintf('%sSICK_NONE_2014-01-06_%s.txt',path_source,filename_source);
            [data_i, data_z, data_n, data_r, data_a, data_e]  = ...
                read_lidar2(filepath_source);
            % Convert scan angle & rotation stage
            [data_a, data_e] = az_to_360_2(data_a,data_e, data_z);
            data_a = -(data_a-pi/2);
            data_a(data_a<0) = data_a(data_a<0) + 360;
            [data_x,data_y, data_z] = sph2cart(deg2rad(data_a),deg2rad(data_e),data_r); % Note coordinate flip
            data_xy = sqrt((data_x).^2+(data_y).^2);
            % Crop
            dist_from_ctr = sqrt((data_x-x_ctr(r)).^2 + (data_y-y_ctr(r)).^2 );
            index = (dist_from_ctr<.3)&(data_z>z_lo(r))&(data_z<z_hi(r));
            data_i = data_i(index);
            data_x = data_x(index);
            data_y = data_y(index);
            data_z = data_z(index);
            
            %{
            figure; 
            title(sprintf('Range = %2.0f', range(r)));
            scatter3(data_x,data_y,data_z,20,'b','filled');
            daspect([1 1 1]); view([0 0]);
            title(sprintf('Range = %2.0f', range(r)));
            %}
            
            data_n = data_n(index);
            data_xy = data_xy(index);
            data_r = data_r(index);
            data_a = data_a(index);
            data_e = data_e(index);
            % Save
            save(filepath_data_i, 'data_i');
            save(filepath_data_x, 'data_x');
            save(filepath_data_y, 'data_y');
            save(filepath_data_z, 'data_z');
            save(filepath_data_n, 'data_n');
            save(filepath_data_xy,'data_xy');
            save(filepath_data_r, 'data_r');
            save(filepath_data_a, 'data_a');
            save(filepath_data_e, 'data_e');
            clear index_far index_near index
        end
        
        %{
        h = figure;
        scatter3(data_x,data_y,data_z, 10,data_i,'filled');
        hold on;
        scatter3(x_ctr(r),y_ctr(r),z_ctr,5000,'b')
        colormap('jet')
        title(sprintf('Range %2.0fm; DBH %2.0fcm',range(r),dbh(d)));
        view(0,90)
        xlim([x_ctr(r)-.5 x_ctr(r) + .5]);
        ylim([y_ctr(r)-.5 y_ctr(r) + .5]);
        %campos([0 0 0]);
        daspect([1 1 1]);
        foo = 1;
       % delete(h);
        %}
        
        for iter = 1:n_iter;
            
           % fprintf('\n**************\n');
           % fprintf('iter  = %3.0f \n', iter);
            
            block_axes = [x_ctr(r)-.3 x_ctr(r) + .3 y_ctr(r)-.3 y_ctr(r) + .3 -.65 1];
            tick = [];
            t_error_line = 0.01;
            t_error_circ = 0.08;
            t_far = .5;
            t_theta = 65;
            t_sampling = .0044; %tand(0.25) = 0.0044;
            info_z_step = 0.1;
            t_kasa = 0.03;
            [r_alpha, r_kasa, r_pratt ] = model_cyl10_dbhtests(data_x,data_y,data_z,block_axes,tick,...
                t_error_line,t_error_circ,t_far,t_theta,t_sampling,info_z_step,t_kasa);
            if numel(r_alpha);
                alpha(iter,r,d) = 2*r_alpha;
                kasa(iter,r,d) = 2*r_kasa;
                pratt(iter,r,d) = 2*r_pratt;
            end
            
        end
    end
end


foo = 1
