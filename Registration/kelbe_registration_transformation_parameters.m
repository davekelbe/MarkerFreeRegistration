%% Validation: transformation parameters 
% Clear without breakpoints 
tmp = dbstatus;
save('tmp.mat','tmp');
clear all; close all; clc;
load('tmp.mat');
dbstop(tmp);
clear tmp;
delete('tmp.mat');
% Default settings 
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = true;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_temp = 'D:\Users\djk2312\Documents\MATLAB\temp\';
load('tempvalidation.mat');
% Initialize program-level options and paths
% Outputs
% G_npath               number of nodes in each path
% G_rx                  array of x angles [degrees]
% G_ry                  array of y angles [degrees]
% G_rz                  array of z angles [degrees]
% G_tx                  array of x translations [m]
% G_ty                  array of y translations [m]
% G_tz                  array of z translations [m]
% info                  program level information
% ************************
% options               options
% paths                 paths
%% Rxyz vs # nodes 
options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_rxyz_vs_node(G_npath,G_rx, G_ry, G_rz,info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );


%% Txyz vs. # nodes 
options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_txyz_vs_node(G_npath,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% Rxyz vs path distance 

options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_rxyz_vs_node(G_npath,G_rx, G_ry, G_rz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% Txyz vs path distance 

options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_txyz_vs_node(G_npath,G_tx, G_ty, G_tz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );

%% Rxyz vs path distance 

options_matlabfig = true;
options_tikzfig = true;
save_tikz = true;
kelbe_reg_rxyz_vs_node(G_npath,G_rx, G_ry, G_rz, info_site,...
    path_tikz, options_matlabfig, options_tikzfig, save_tikz );


%% Compare tree locations 

G_tree_hat = cell(n_g,1);
G_tree_true = cell(n_g,1);
GR4_eye = cellfun(@double,repmat({[eye(3) zeros(3,1); zeros(1,3) 1]},n_g,1),'Un',0);
for s = 1:n_S;
    is_valid = (G_start==s);
    P4_temp = [P_LCS{s}; ones(1,size(P_LCS{s},2))];
    GR4_temp = G_R4(is_valid);
    Geye_temp = GR4_eye(is_valid);
    G_tree_hat(is_valid) = cellfun(@(x) x*P4_temp, GR4_temp, 'uniformoutput',false);
    G_tree_true(is_valid) =  cellfun(@(x) x*P4_temp, Geye_temp, 'uniformoutput',false);
end

n_inG = cellfun(@(x) size(x,2), G_tree_hat);
n_all = sum(n_inG);

all_tree_hat = zeros(n_all,3);
all_tree_true = zeros(n_all,3); 
all_npath = zeros(n_all,1);
all_dist = zeros(n_all,1);
all_minmatch = zeros(n_all,1);
start_ix = [0; cumsum(n_inG)]+1;
end_ix = circshift(start_ix, [-1,0])-1;
start_ix = start_ix(1:end-1);
end_ix = end_ix(1:end-1);

for g = 1:n_g;
    all_tree_hat(start_ix(g):end_ix(g),:) = G_tree_hat{g}(1:3,:)';
    all_tree_true(start_ix(g):end_ix(g),:) = G_tree_true{g}(1:3,:)';
    all_npath(start_ix(g):end_ix(g)) = repmat(G_npath(g),[n_inG(g),1]);
    all_dist(start_ix(g):end_ix(g)) = repmat(G_dist(g),[n_inG(g),1]);
    all_minmatch(start_ix(g):end_ix(g)) = repmat(G_min_nmatch(g),[n_inG(g),1]);
end

all_tree_exy = sqrt(sum((all_tree_hat(:,1:2) - all_tree_true(:,1:2)).^2,2)); 
all_tree_ex = all_tree_hat(:,1) - all_tree_true(:,1); 
all_tree_ey = all_tree_hat(:,2) - all_tree_true(:,2); 
all_tree_ez = all_tree_hat(:,3) - all_tree_true(:,3); 

%% Tree location vs. number of nodes 

%{
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
ylabel('Error in xy');
%title('Tree position accuracy after registration');
filename = sprintf('Error-tloc_vs_nodes_s%02.0f',info_site);
filepath_save = sprintf('%s%s.eps',path_save,filename);
%saveas(gcf,filepath_save,'psc2')
%}

% Make tikz figure 
%{
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
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%% Tree location vs path distance 

%{
clear legend_str
figure;
hold on
plot(all_dist,all_tree_ex,'ok','markerfacecolor','r');
plot(all_dist,all_tree_ey,'ok','markerfacecolor','g');
%plot(all_dist,all_tree_ez,'ok','markerfacecolor','b');
axisval = axis;
legend_str{1} = 'x';
legend_str{2} = 'y';
%legend_str{3} = 'z';
legend(legend_str);
xlabel('Total path distance');
ylabel('Error');
%title('Tree position accuracy after registration');
filename = sprintf('Error-tloc_vs_dist_s%02.0f',info_site);
filepath_save = sprintf('%s%s.eps',path_save,filename);
saveas(gcf,filepath_save,'psc2');
%}

% Make tikz figure 
%{
dist_axes = 10:10:90;
n_d = numel(dist_axes) - 1;
data_1 = cell(n_d,1);
data_2 = cell(n_d,1);
data_3 = cell(n_d,1);
labels = cell(n_d,1);
ctr = 1;
for d = 1:n_d;
    is_valid = (all_dist>=dist_axes(d)) & (all_dist<dist_axes(d+1));
    if sum(is_valid);
        data_1{ctr} = all_tree_ex(is_valid);
        data_2{ctr} = all_tree_ey(is_valid);
        data_3{ctr} = all_tree_ez(is_valid);
        labels{ctr} = sprintf('%g-%g',dist_axes(d), dist_axes(d+1));
        ctr = ctr + 1;
    end
end

tikzxlabel = 'Distance in meters';
tikzylabel = 'Error in $m$';
legendlabels = {'\gls{tx}', '\gls{ty}', '\gls{tz}'};

ymin = min([all_tree_ex; all_tree_ey; all_tree_ez]);
ymax = max([all_tree_ex; all_tree_ey; all_tree_ez]);
ylimval = [ymin ymax];

filename = sprintf('Error-tloc_vs_dist_s%02.0f',info_site);
filepath_tikz = sprintf('%s%s.tex',path_tikz,filename);
fid = fopen(filepath_tikz,'w+');
makeTikzBoxplot3( data_1, data_2, data_3,'fid', fid, ...
    'labels', labels,'xlabel', tikzxlabel,'ylabel', tikzylabel,...
    'legendlabels', legendlabels,...
    'boxsep',2.25,...
    'ylim', ylimval);
fclose all;
%}
clear Gminpath Gmaxpath tick legend_str filename filepath_save
clear data1 data2 data3 data_isvalid p 
clear tikzxlabel tikzylabel legendlabels labels i ymin myax ylimval 
clear filepath_tikz fid 
%%
%
clear legend_str
figure;
hold on
for g = 1:n_g;
plot(G_tree_true{g}(1,:),G_tree_hat{g}(1,:),'+','markeredgecolor',[.5 .5 .5]);
plot(G_tree_true{g}(2,:),G_tree_hat{g}(2,:),'x','markeredgecolor',[.8 .8 .8]);
%plot(G_tree_true{g}(3,:),G_tree_hat{g}(3,:),'+','markeredgecolor','b');
end
axisval = axis;
plot(axisval(1:2),axisval(3:4),'-k');
legend_str{1} = 'x';
legend_str{2} = 'y';
%legend_str{3} = 'z';
legend(legend_str);
xlabel('Truth location');
ylabel('Estimated location');
%title('Tree position accuracy after registration');
filepath_save = sprintf('%sError-treeloc_s%02.0f.eps',path_save, info_site);
%saveas(gcf,filepath_save,'psc2')
%}

    for p = 1:G_npath(g);
        path = G_path{g}(p,:);
        i = path(1);
        j = path(2);
        loop_P_true{g}{p} = P_LCS{i}(:,match_i_sub2{i,j});
        loop_P_hat{g}{p} = loop_R{g}{p}*P_LCS{i}(:,match_i_sub2{i,j}) +...
            repmat(loop_t{g}{p},[1,numel(match_i_sub2{i,j})]);
    end


all_true = [];
all_hat = [];
all_g1 = [];
all_npath = [];
for g = 1:n_S;
    for p = 1:G_npath(g);
       all_true = [all_true loop_P_true{g}{p}];
       all_hat = [all_hat loop_P_hat{g}{p}];
    end
end
all_error = sqrt((all_true - all_hat).^2); 

fprintf('\nRMSE x = %3.3f m\n', mean(all_error(1,:)))
fprintf('RMSE y = %3.3f m\n', mean(all_error(2,:)))
fprintf('RMSE z = %3.3f m\n', mean(all_error(3,:)))

%
figure;
hold on
plot(all_true(1,:),all_hat(1,:),'ok','markerfacecolor','r');
plot(all_true(2,:),all_hat(2,:),'ok','markerfacecolor','g');
plot(all_true(3,:),all_hat(3,:),'ok','markerfacecolor','b');
axisval = axis;
plot(axisval(1:2),axisval(3:4),'-k');
legend_str{1} = 'x';
legend_str{2} = 'y';
legend_str{3} = 'z';
legend(legend_str);
xlabel('Truth location');
ylabel('Estimated location');
title('Tree position accuracy after registration');
filepath_save = sprintf('%sError-treeloc_s%02.0f.eps',path_save, info_site);
saveas(gcf,filepath_save,'psc2')
%}
%% Validation: Point cloud error 
% Initial transformation of PLY
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
clear match_teff_arr vertex data_x2 data_y2 data_z2 color_ply xyz2t
%}

%% Validation: Error field 

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;
space =2;
xv = xmin:space:xmax;
yv = ymin:space:ymax;
zv = zmin:space:zmax;
[X,Y,Z] = meshgrid(xv, yv, zv);
[nx,ny,nz] = size(X);
test_rx = 0;
test_ry = 3;
test_rz = 1;
test_tx = .1;
test_ty = 0.5;
test_tz = 0.0;
test_R = compose_rotation(deg2rad(test_rx),deg2rad(test_ry),deg2rad(test_rz));
test_t = [test_tx test_ty test_tz]';
XYZ = [X(:),Y(:),Z(:)]';
n_pts = size(XYZ,2);
XYZhat = test_R*XYZ + repmat(test_t,[1, n_pts]);
error_XYZ = XYZ - XYZhat;
error = abs(mean((XYZ - XYZhat),1));

figure;
hold on
scatter3(XYZ(1,:),XYZ(2,:),XYZ(3,:),10,'r','filled');
scatter3(XYZhat(1,:),XYZhat(2,:),XYZhat(3,:),10,'b','filled');
axis equal
legend('True','Estimated');

figure;
hold on
scatter3(XYZ(1,:),XYZ(2,:),XYZ(3,:),10,error,'filled');
axis equal
colorbar
emin = min(error);
emax = max(error); 
cmin = emin -.1;
cmax = emax + .1;
caxis([cmin, cmax]);

figure;
xbin = linspace(min(error),max(error),100);
count = hist(error, xbin);
plot(xbin,count,'-x');

figure; 
scale = 2;
quiver3(XYZ(1,:),XYZ(2,:),XYZ(3,:),error_XYZ(1,:), error_XYZ(2,:), error_XYZ(3,:),...
    scale);


