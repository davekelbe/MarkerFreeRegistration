%% Validation of registration - wiggle
% This script assesses the error of registration independent of error in
% the tie point sets.
% A tie point set is rotated and translated by random angles/distances,
% and wiggled.
% The transformed point set is matched to the original set, giving Rhat
% The extent to which R~=Rhat describes the error in registration
% This should be aggregated for multiple tie points (say, all plots of 31)
% We evaluate:
%   error vs. rotation angle
%   error vs. distance
%
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = false;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
%% Load points
fprintf('\nLoad points\n');

info_exp = 'Harvard';
info_suffix = '03-01';
info_slash = '\';
info_site = 31;
path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
D = dir(path_site);
ctr = 1;
info_valid_plot = {'13'};
for d = 1:numel(D);
    if strcmp(D(d).name(1), '.')
        continue
    end
    info_plot = D(d).name;
    if ~any(strcmp(info_plot, info_valid_plot))
        continue
    end
    path_mat = sprintf('%s%s%smat%s',path_site,info_plot,info_slash,info_slash);
    filepath_tree = sprintf('%stree.mat',path_mat);
    if ~exist(filepath_tree,'file')
        continue
    end
    load(filepath_tree);
    n_tree = numel(tree);
    plot = str2num(info_plot);
    P_LCS{ctr} = nan(3,n_tree);
    P_rad{ctr} = nan(n_tree,1);
    P_plot(ctr) = plot;
    for t = 1:n_tree;
        P_LCS{ctr}(:,t) = tree(t).loc(:,1);
        P_rad{ctr}(t) = tree(t).r(1);
    end
    path_ply{ctr} = sprintf('%s%s%sply%s',path_site,info_plot,info_slash,info_slash);
    filepath_ply{ctr} = sprintf('%spoints_full_%03.0f-%02.0f.ply', path_ply{ctr}, info_site, plot);
    ctr = ctr + 1;
end
P_n = cellfun(@numel,P_rad);
n_S = numel(P_LCS);

i_xmin = min(cellfun(@(x) min(x(1,:)),P_LCS));
i_xmax = max(cellfun(@(x) max(x(1,:)),P_LCS));
i_ymin = min(cellfun(@(x) min(x(2,:)),P_LCS));
i_ymax = max(cellfun(@(x) max(x(2,:)),P_LCS));
i_zmin = min(cellfun(@(x) min(x(3,:)),P_LCS));
i_zmax = max(cellfun(@(x) max(x(3,:)),P_LCS));

% Colormap for sensors
color = jet(n_S);

% Individual camera views
if false;%options_verbose && options_imagepoints;
    for s = 1:n_S;
        clear legend_str
        figure
        hold on
        plot3(0,0,0,'^k','markersize',10,...
            'markerfacecolor',color(s,:));
        hdummy = plot3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),'ok','markersize',5,...
            'markerfacecolor',color(s,:));
        set(hdummy, 'visible', 'off');
        for t = 1:numel(P_rad{s});
            h = filledCircle([P_LCS{s}(1,t); P_LCS{s}(2,t)]',P_rad{s}(t),1000,color(s,:));
        end
        %scatter3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),30,...
        %    color_P_index(truth_P_index{s},:),'filled');
        %axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
        axis auto
        axisval = axis;
        %set(gca, 'xtick',
        xlabel('x');ylabel('y');zlabel('z');
        view(0,90);
        grid on
        %titlestr = sprintf('Scan %g',P_plot(s));
        %title(titlestr);
        legend_str{1} = sprintf('Scanner %g',P_plot(s));
        legend_str{2} = sprintf('Stem map in LCS %g', P_plot(s));
        legend(legend_str,'location','northeast');
        legend boxoff       % Hides the legend's axes (legend border and background)
        set(gca, 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        filepath_save = sprintf('%sLCS_%02.0f.eps',path_save, P_plot(s));
        saveas(gcf,filepath_save,'psc2')
    end
end
%}
% Outputs
% P_LCS                 points in local coordinate system (cell)
% P_rad                 radius of points
% P_plot                plot number
% P_n                   number of points
% i_xmin                minimum x value of data
% i_zmax                maximum z value of data
% color                 colormap for scans
% info                  program level information
% options               program level options
% path                  paths
% n                     number of given variable
clear D ctr d filepath_tree plot t tree
%% Run registration code

% Define array of rotation and translation steps 
% wiggle in rotation angles [m]
rmin = 0;
rstep = 2;
rmax = 10;
% wiggle in transformation parameters [m]
tmin = 0;
tstep = 2;
tmax = 10;
% wiggle in tree locations (range, [m])
wmin = 0.1;
wstep = 0.1; 
wmax = .5;
% Wiggle in percentage of points with wiggle
pmin = 100;
pstep = 10;
pmax = 100;

wminxyz = sqrt(wmin^2/3);
wstepxyz = sqrt(wstep^2/3);
wmaxxyz = sqrt(wmax^2/3);

d_rx = rmin:rstep:rmax;
d_ry = rmin:rstep:rmax;
d_rz = rmin:rstep:rmax;
d_tx = tmin:tstep:tmax;
d_ty = tmin:tstep:tmax;
d_tz = tmin:tstep:tmax;
d_w = wmin:wstep:wmax;
d_wx = wminxyz:wstepxyz:wmaxxyz;
d_wy = wminxyz:wstepxyz:wmaxxyz;
d_wz = wminxyz:wstepxyz:wmaxxyz;
d_p = pmin:pstep:pmax;

n_rx = numel(d_rx);
n_ry = numel(d_ry);
n_rz = numel(d_rz);
n_tx = numel(d_tx);
n_ty = numel(d_ty);
n_tz = numel(d_tz);
n_wx = numel(d_wx);
n_wy = numel(d_wy);
n_wz = numel(d_wz);
n_p = numel(d_p);

P_LCSi = P_LCS{1}; 
P_radi = P_rad{1};
P_radj = P_rad{1};


%% Register wiggled points to unwiggled points 
% No rotation/translation
% 100% of tree locations wiggled 
% No error in diameter estimation 

n_trial = 200;
rhat = zeros(n_wx, n_trial, 3);
that = zeros(n_wx, n_trial, 3);

%for p = 1:n_p; % Percentage of points wiggled 
%        fprintf('\nPercentage %2.0f of %2.0f', p, n_p); 
tic 
    for w = 1:1;%n_wx % Deviation of wiggle 
        fprintf('\n\tWiggle %2.0f of %2.0f', w, n_wx); 
        for t = 1:n_trial;
            fprintf('\n\t\tTrial %2.0f of %2.0f', t, n_trial); 
            %percentage = 1;%d_p(p)/100;
            %n_wiggle = floor(n_tree*percentage);
            stdval = d_wx(w);
            noise = stdval*randn(3,n_wiggle);
            %ix = randperm(n_tree);
            %ix = ix(1:n_wiggle);
            %P_noise = zeros(3,n_tree);
            %P_noise(:,ix) = noise;
            P_LCSj = P_LCSi + P_noise;
            [Rhat,that, ~, ~ ] = toy_registrationfunction(P_LCSi',P_LCSj',P_radi,P_radj);
%             rhat(w,t,1) = rad2deg(decompose_rotation_rx(Rhat));
%             rhat(w,t,2) = rad2deg(decompose_rotation_ry(Rhat));
%             rhat(w,t,3) = rad2deg(decompose_rotation_rz(Rhat));
%             that(w,t,1) = that(1);
%             that(w,t,2) = that(2);
%             that(w,t,3) = that(3);           
        end
    end
fprintf('\nTime = %2.0f\n', toc);

tic 
    for w = 1:1;%n_wx % Deviation of wiggle 
        fprintf('\n\tWiggle %2.0f of %2.0f', w, n_wx); 
        parfor t = 1:n_trial;
            fprintf('\n\t\tTrial %2.0f of %2.0f', t, n_trial); 
            %percentage = 1;%d_p(p)/100;
            %n_wiggle = floor(n_tree*percentage);
            stdval = d_wx(w);
            noise = stdval*randn(3,n_tree);
            %ix = randperm(n_tree);
            %ix = ix(1:n_wiggle);
            %P_noise = zeros(3,n_tree);
            %P_noise(:,ix) = noise;
            P_LCSj = P_LCSi + noise;
            [Rhat,that, ~, ~ ] = toy_registrationfunction(P_LCSi',P_LCSj',P_radi,P_radj);
%             rhat(w,t,1) = rad2deg(decompose_rotation_rx(Rhat));
%             rhat(w,t,2) = rad2deg(decompose_rotation_ry(Rhat));
%             rhat(w,t,3) = rad2deg(decompose_rotation_rz(Rhat));
%             that(w,t,1) = that(1);
%             that(w,t,2) = that(2);
%             that(w,t,3) = that(3);           
        end
    end
fprintf('\n Parallel Time = %2.0f\n', toc);
    
%end

sensor_offset = sqrt(that(:,:,1).^2 + that(:,:,2).^2 + that(:,:,3).^2);
sensor_angle = sqrt(rhat(:,:,1).^2 + rhat(:,:,2).^2 + rhat(:,:,3).^2);
sensor_offset_mean = mean(sensor_offset,2);
sensor_angle_mean = mean(sensor_angle,2);
sensor_offset_std = std(sensor_offset,0,2);
sensor_angle_std = std(sensor_angle,0,2);

figure; 
scatter(d_w, sensor_angle_mean, 'filled');
hold on
errorbar(d_w, sensor_angle_mean, sensor_angle_std);
xlabel('Standard deviation of Tree locations [m]');
ylabel('Error in sensor angle [deg]')
axis([wmin wmax 0 150]);
grid on
%%
