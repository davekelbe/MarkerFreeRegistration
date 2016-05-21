%% Validation of registration
% This script assesses the error of registration independent of error in
% the tie point sets.
% A tie point set is rotated and translated by random angles/distances, R
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
rmin = 0;
rstep = 2;
rmax = 10;
tmin = 0;
tstep = 2;
tmax = 10;
d_rx = rmin:rstep:rmax;
d_ry = rmin:rstep:rmax;
d_rz = rmin:rstep:rmax;
d_tx = tmin:tstep:tmax;
d_ty = tmin:tstep:tmax;
d_tz = tmin:tstep:tmax;
n_rx = numel(d_rx);
n_ry = numel(d_ry);
n_rz = numel(d_rz);
n_tx = numel(d_tx);
n_ty = numel(d_ty);
n_tz = numel(d_tz);

P_LCSi = P_LCS{1}; 
P_radi = P_rad{1};
P_radj = P_rad{1};

% rx
fprintf('\n Working on rx\n');
action.rx = true; action.ry = false; action.rz = false;
action.tx = false; action.ty = false; action.tz = false;
[ Ptrue.rx, Phat.rx ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% ry
fprintf('\n Working on ry\n');
action.rx = false; action.ry = true; action.rz = false;
action.tx = false; action.ty = false; action.tz = false;
[ Ptrue.ry, Phat.ry ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% rz
fprintf('\n Working on rz\n');
action.rx = false; action.ry = false; action.rz = true;
action.tx = false; action.ty = false; action.tz = false;
[ Ptrue.rz, Phat.rz ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% tx
fprintf('\n Working on tx\n');
action.rx = false; action.ry = false; action.rz = false;
action.tx = true; action.ty = false; action.tz = false;
[ Ptrue.tx, Phat.tx ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% ty
fprintf('\n Working on ty\n');
action.rx = false; action.ry = false; action.rz = false;
action.tx = false; action.ty = true; action.tz = false;
[ Ptrue.ty, Phat.ty ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% tz
fprintf('\n Working on tz\n');
action.rx = false; action.ry = false; action.rz = false;
action.tx = false; action.ty = false; action.tz = true;
[ Ptrue.tz, Phat.tz ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% rxy
fprintf('\n Working on rxy\n');
action.rx = true; action.ry = true; action.rz = false;
action.tx = false; action.ty = false; action.tz = false;
[ Ptrue.rxy, Phat.rxy ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% rxyz
fprintf('\n Working on rxyz\n');
action.rx = true; action.ry = true; action.rz = true;
action.tx = false; action.ty = false; action.tz = false;
[ Ptrue.rxyz, Phat.rxyz ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% txy
fprintf('\n Working on txy\n');
action.rx = false; action.ry = true; action.rz = false;
action.tx = true; action.ty = true; action.tz = false;
[ Ptrue.txy, Phat.txy ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% txyz
fprintf('\n Working on txyz\n');
action.rx = true; action.ry = true; action.rz = true;
action.tx = false; action.ty = false; action.tz = false;
[ Ptrue.txyz, Phat.txyz ] = ...
    determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );
% rxyztxyz
% fprintf('\n Working on rxyztxyz\n');
% action.rx = true; action.ry = true; action.rz = true;
% action.tx = true; action.ty = true; action.tz = true;
% [ Param.rxyztxyz.rx, Param.rxyztxyz.ry, Param.rxyztxyz.rz, Param.rxyztxyz.tx, Param.rxyztxyz.ty, Param.rxyztxyz.tz ] = ...
%     determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action );



for i = 1:n_S;
    for j = i:n_S;
        fprintf('\n\tMatching %g to %g\n',j,i);
        [ a,b, c, d ] = toy_registrationfunction(P_LCS{i}',P_LCS{j}',P_rad{i},P_rad{j});
        if ~isempty(a)&&~isempty(b)&&~isempty(c)&&~isempty(d);
            match_R{i,j} = a;
            match_t{i,j} = b;
            match_i{i,j} = c;
            match_j{i,j} = d;
        end
    end
end








