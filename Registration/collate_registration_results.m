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
all_site = [];
all_plotI = [];
all_plotJ = [];
all_nit = [];
all_RMSE = [];
all12_rx = [];
all12_ry = [];
all12_rz = [];
all12_tx = [];
all12_ty = [];
all12_tz = [];
all_tx = [];
all_ty = [];
all_tz = [];
all2_tx = [];
all2_ty = [];
all2_tz = [];
all_true1 = [];
all_true2 = [];
n_plot = zeros(n_site,1);

exclude_site = [];
for s = 1:11;%n_site;
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
    load(filepath_nit);
    filepath_match12_RMSE = sprintf('%s%s',path_mat, 'match12_RMSE.mat');
    load(filepath_match12_RMSE);
    filepath_match12_R = sprintf('%s%s',path_mat, 'match12_R.mat');
    load(filepath_match12_R);
    filepath_match12_t = sprintf('%s%s',path_mat, 'match12_t.mat');
    load(filepath_match12_t);
    filepath_match1_t = sprintf('%s%s',path_mat, 'match1_t.mat');
    load(filepath_match1_t);
    filepath_match2_t = sprintf('%s%s',path_mat, 'match2_t.mat');
    load(filepath_match2_t);
    filepath_match_plotI = sprintf('%s%s',path_mat, 'match_plotI.mat');
    load(filepath_match_plotI);
    filepath_match_plotJ = sprintf('%s%s',path_mat, 'match_plotJ.mat');
    load(filepath_match_plotJ);
    match12_rx = 180*cellfun(@decompose_rotation_rx,match12_R)/pi;
    match12_ry = 180*cellfun(@decompose_rotation_ry,match12_R)/pi;
    match12_rz = 180*cellfun(@decompose_rotation_rz,match12_R)/pi;
    match12_tx = cellfun(@decompose_translation_tx,match12_t);
    match12_ty = cellfun(@decompose_translation_ty,match12_t);
    match12_tz = cellfun(@decompose_translation_tz,match12_t);
    match1_tx = cellfun(@decompose_translation_tx,match1_t);
    match1_ty = cellfun(@decompose_translation_ty,match1_t);
    match1_tz = cellfun(@decompose_translation_tz,match1_t);
    match2_tx = cellfun(@decompose_translation_tx,match2_t);
    match2_ty = cellfun(@decompose_translation_ty,match2_t);
    match2_tz = cellfun(@decompose_translation_tz,match2_t);
    
    n_plot(s) = size(match_nit,1);
    
    true1 = zeros(n_plot(s),n_plot(s));
    true2 = zeros(n_plot(s),n_plot(s));
    for p = 1:n_plot(s);
        filepath_true1 = sprintf('%s%02.0f%smat%s%s%s',path_site,p,info_slash,info_slash, 'match1_tf.mat');
        filepath_true2 = sprintf('%s%02.0f%smat%s%s%s',path_site,p,info_slash,info_slash, 'match2_tf.mat');
        load(filepath_true1);
        load(filepath_true2);
        true1(p,:) = match1_tf;
        true2(p,:) = match2_tf;
    end
    
    all_true1 = [all_true1; true1(:)];
    all_true2 = [all_true2; true2(:)];
    
    %{
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    [x_sph, y_sph, z_sph] = sphere;
    figure
    hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    % Plot all pairwise matches from j
    i = 1;
    for j = 1:n_plot(s);
        % if ~isempty0(i,j);
        x_axist = match12_R{i,j}*x_axis + repmat(match12_t{i,j},[1,2]);
        y_axist = match12_R{i,j}*y_axis + repmat(match12_t{i,j},[1,2]);
        z_axist = match12_R{i,j}*z_axis + repmat(match12_t{i,j},[1,2]);
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
        %   end
    end
    %}
    
    
    all_site = [all_site; repmat(info_site, [n_plot(s).^2,1])];
    all_plotI = [all_plotI; match_plotI(:)];
    all_plotJ = [all_plotJ; match_plotJ(:)];
    all_nit = [all_nit; match_nit(:)];
    all_RMSE = [all_RMSE; match12_RMSE(:)];
    all12_rx = [all12_rx; match12_rx(:)];
    all12_ry = [all12_ry; match12_ry(:)];
    all12_rz = [all12_rz; match12_rz(:)];
    all12_tx = [all12_tx; match12_tx(:)];
    all12_ty = [all12_ty; match12_ty(:)];
    all12_tz = [all12_tz; match12_tz(:)];
    all_tx = [all_tx; match1_tx(:)];
    all_ty = [all_ty; match1_ty(:)];
    all_tz = [all_tz; match1_tz(:)];
    all2_tx = [all2_tx; match2_tx(:)];
    all2_ty = [all2_ty; match2_ty(:)];
    all2_tz = [all2_tz; match2_tz(:)];
end

all_t = sqrt(all_tx.^2 + all_ty.^2 + all_tz.^2);
all2_t = sqrt(all2_tx.^2 + all2_ty.^2 + all2_tz.^2);

clear filename_csv fid C s info_site path_site path_mat filepath_nit
clear filepath_match12_RMSE filepath_match12_R filepath_match12_t filepath_match1_t
clear filepath_match_plotI filepath_match_plotJ
clear match12_R match12_RMSE match12_t match1_t match_nit match_plotI match_plotJ
clear match12_rx match12_ry match12_rz match12_tx match12_ty match12_tz
clear match1_tx match1_ty match1_tz

path_reg_results = 'D:\Users\djk2312\Documents\Harvard\reg_results\';
filepath_reg_results = sprintf('%s%s',path_reg_results,'collate_results.mat');
%save(filepath_reg_results);

% Histogram of RMSE values
all_true1nan = all_true1;
all_true1nan(isnan(all_true1nan))= 0;
all_true1nan = logical(all_true1nan);
all_true2nan = all_true2;
all_true2nan(isnan(all_true2nan))= 0;
all_true2nan = logical(all_true2nan);

binranges = 0:1:max(all_RMSE);
bincountsTT = histc(all_RMSE(all_true1nan&all_true2nan), binranges);
bincountsFF = histc(all_RMSE(~all_true1nan&~all_true2nan), binranges);
bincountsTF = histc(all_RMSE(all_true1nan&~all_true2nan), binranges);
bincountsFT = histc(all_RMSE(~all_true1nan&all_true2nan), binranges);
bincountsF = histc(all_RMSE(~all_true1nan|~all_true2nan), binranges);

bincountsTT = uint16(bincountsTT);
bincountsFF = uint16(bincountsFF);
bincountsTF = uint16(bincountsTF);
bincountsFT = uint16(bincountsFT);
bincountsF = uint16(bincountsF);

% linear just histogram
%{
figure;
clear legend_str
hold on;
plot(binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
%semilogy(binranges, bincountsTT, '-k', 'linewidth', 2);
%semilogy(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
%semilogy(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
%semilogy(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
%plot(binranges, bincountsF, '-', 'color', [.8 .8 .8], 'linewidth', 1.4);

ylabel('Count');
xlabel('RMSE [m]');
legend(legend_str)
filepath_countvsrmse = 'Z:\Desktop\countvsrmse-crop.tex';
%matlab2tikz(filepath_countvsrmse)
%}

% log log just histogram
%{
figure;
clear legend_str
hold on;
loglog(binranges, bincountsTT, '-k', 'linewidth', 2);
legend_str{1} = 'TT';
loglog(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
loglog(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
loglog(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
%semilogy(binranges, bincountsTT, '-k', 'linewidth', 2);
%semilogy(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
%semilogy(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
%semilogy(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
%plot(binranges, bincountsF, '-', 'color', [.8 .8 .8], 'linewidth', 1.4);

ylabel('Count');
xlabel('RMSE [m]');
legend(legend_str)
%}


%% ROC curve log log 
%{
[rocx, rocy, roct] = perfcurve(all_true1nan&all_true2nan, all_RMSE,false);
lw = 1.5;
figure;
subplot(2,2,1);
semilogy(roct, rocy, '-k', 'linewidth',lw)
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
loglog(rocx, rocy, '-k', 'linewidth',lw)
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
loglog(rocx, roct, '-k', 'linewidth',lw)
xlabel('Incorrect Detection Rate');
ylabel('RMSE threshold [m]');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
subplot(2,2,1);
%ylim([.0003, 1.5]);
%axis(axis0);
subplot(2,2,2);
%ylim([.0003, 1.5])
%axis(axis0);
%axis(axis0);
clear legend_str
subplot(2,2,3);
hold on
plot(binranges, bincountsTT, '-k', 'linewidth',lw);
legend_str{1} = 'TT';
plot(binranges, bincountsF, '--', 'color', [0 0 0], 'linewidth',lw);
legend_str{2} = 'F';
plot(binranges, bincountsFF, '-', 'color', [.6 .6 .6], 'linewidth',lw);
legend_str{3} = 'FF';
plot(binranges, bincountsTF, ':', 'color', [.6 .6 .6], 'linewidth',lw);
legend_str{4} = 'TF';
plot(binranges, bincountsFT, '-.', 'color', [.6 .6 .6], 'linewidth',lw);
legend_str{5} = 'FT';
%semilogy(binranges, bincountsTT, '-k', 'linewidth', 2);
%semilogy(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
%semilogy(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
%semilogy(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
%plot(binranges, bincountsF, '-', 'color', [.8 .8 .8], 'linewidth', 1.4);
ylabel('Count');
xlabel('RMSE [m]');
legend(legend_str)
filepath_countvsrmse = 'Z:\Desktop\countvsrmselog.tex';
matlab2tikz(filepath_countvsrmse)
%}

% ROC curve test
%{
mu = 5;
sigma = 3;
m = 1;
n = 500;
data1 = normrnd(mu,sigma,m,n);
label1 = repmat(true,[m,n]);
mu = 10;
sigma = 3;
m = 1;
n = 500;
data2 = normrnd(mu,sigma,m,n);
label2 = repmat(false,[m,n]);
data = [data1 data2];
label = [label1 label2];
% Histogram
binranges = min(data):1:max(data);
bincountsTT = histc(data(label==true), binranges);
bincountsFF = histc(data(label==false), binranges);
[rocx, rocy, roct] = perfcurve(label, data,false);
figure;
subplot(2,2,1);
plot(roct, rocy, '-k')
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
plot(rocx, rocy, '-k')
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(roct, rocx, '-k', 'linewidth',1.5)
xlabel('RMSE threshold [m]');
ylabel('Incorrect Detection Rate');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
clear legend_str
subplot(2,2,3);
hold on;
plot(binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
%}

% ROC curve linear 
%{
labels = all_true1nan&all_true2nan;
%figure;
%scatter(all_RMSE, labels, 10, labels);
[rocx, rocy, roct] = perfcurve(labels, all_RMSE,true);
figure;
subplot(2,2,1);
plot(roct, rocy, '-k')
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
plot(rocx, rocy, '-k')
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(roct, rocx, '-k', 'linewidth',1.5)
xlabel('RMSE threshold [m]');
ylabel('Incorrect Detection Rate');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
clear legend_str
subplot(2,2,3);
hold on;
plot(binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
legend(legend_str);
%}

% Fixed normal perfcurve 
%{
lw = 1.5;
labels = all_true1nan&all_true2nan;
%figure;
%scatter(all_RMSE, labels, 10, labels);
[rocx, rocy, roct] = perfcurve(labels, all_RMSE,false);
rocx = 1-rocx;
rocy = 1-rocy;
figure;
subplot(2,2,1);
plot(roct, rocx, '-k', 'linewidth', lw)
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
plot(rocy, rocx, '-k', 'linewidth', lw)
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(rocy, roct, '-k', 'linewidth',lw)
ylabel('RMSE threshold [m]');
xlabel('Incorrect Detection Rate');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
clear legend_str
subplot(2,2,3);
xlabel('RMSE threshold [m]');
ylabel('Count');
hold on;
plot(binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
legend(legend_str);
filepath_countvsrmse = 'Z:\Desktop\countvsrmse.tex';
matlab2tikz(filepath_countvsrmse)
%}

% Fixed normal perfcure 4 separate plots  *IN PAPER
%{
lw = 2.5;
fs = 18;
labels = all_true1nan&all_true2nan;
%figure;
%scatter(all_RMSE, labels, 10, labels);
[rocx, rocy, roct] = perfcurve(labels, all_RMSE,false);
rocx = 1-rocx;
rocy = 1-rocy;
f1 = figure;
plot(roct, rocx, '-k', 'linewidth', lw)
xlabel('RMSE threshold [m]','fontsize', fs);
ylabel('Correct Detection Rate','fontsize', fs);
set(gcs, 'fontsize', fs);
axis1 = axis;
filepath_cdvsrmse = 'Z:\Desktop\cdvsrmse.tex';
matlab2tikz(filepath_cdvsrmse)
f2 = figure; 
plot(rocy, rocx, '-k', 'linewidth', lw)
xlabel('Incorrect Detection Rate','fontsize', fs);
ylabel('Correct Detection Rate','fontsize', fs);
set(gcs, 'fontsize', fs);
axis2 = axis;
filepath_cdvsid = 'Z:\Desktop\cdvsid.tex';
matlab2tikz(filepath_cdvsid)
f3 = figure;
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(roct, rocy, '-k', 'linewidth',lw)
xlabel('RMSE threshold [m],','fontsize', fs);
ylabel('Incorrect Detection Rate','fontsize', fs);
set(gcs, 'fontsize', fs);
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
filepath_idvsrmse = 'Z:\Desktop\idvsrmse.tex';
matlab2tikz(filepath_idvsrmse)
clear legend_str
f4 = figure;
xlabel('RMSE threshold [m]','fontsize', fs);
ylabel('Count','fontsize', fs);
hold on;
plot(binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
legend(legend_str);
set(gcs, 'fontsize', fs);

filepath_countvsrmse = 'Z:\Desktop\countvsrmse.tex';
matlab2tikz(filepath_countvsrmse)
%

% Fixed log perfcurve
%{
labels = all_true1nan&all_true2nan;
%figure;
%scatter(all_RMSE, labels, 10, labels);
[rocx, rocy, roct] = perfcurve(labels, all_RMSE,true);
figure;
subplot(2,2,1);
semilogy(roct, rocy, '-k')
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
plot(rocy, rocx, '-k')
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(roct, rocx, '-k', 'linewidth',1.5)
xlabel('RMSE threshold [m]');
ylabel('Incorrect Detection Rate');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
clear legend_str
subplot(2,2,3);
hold on;
plot(binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
legend(legend_str);
%}

% Reverse perfcurve
%{
labels = all_true1nan&all_true2nan;
%figure;
%scatter(all_RMSE, labels, 10, labels);
[rocx, rocy, roct] = perfcurve(labels, -all_RMSE,true);
figure;
subplot(2,2,1);
plot(roct, rocy, '-k')
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
plot(rocx, rocy, '-k')
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(roct, rocx, '-k', 'linewidth',1.5)
xlabel('RMSE threshold [m]');
ylabel('Incorrect Detection Rate');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
clear legend_str
subplot(2,2,3);
hold on;
plot(-binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(-binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(-binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(-binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
legend(legend_str);
%}

% Fixed reverse perfcurve
%{
labels = all_true1nan&all_true2nan;
%figure;
%scatter(all_RMSE, labels, 10, labels);
[rocx, rocy, roct] = perfcurve(labels, -all_RMSE,true);
figure;
subplot(2,2,1);
plot(-roct, rocy, '-k')
xlabel('RMSE threshold [m]');
ylabel('Correct Detection Rate');
axis1 = axis;
subplot(2,2,2);
plot(rocx, rocy, '-k')
xlabel('Incorrect Detection Rate');
ylabel('Correct Detection Rate');
axis2 = axis;
subplot(2,2,4)
%plot(rocx, roct, '-k', 'linewidth',1.5)
%xlabel('Incorrect Detection Rate');
%ylabel('RMSE threshold [m]');
plot(roct, rocx, '-k', 'linewidth',1.5)
xlabel('RMSE threshold [m]');
ylabel('Incorrect Detection Rate');
axis3 = axis;
axis0 = [min([axis1(1),axis2(1),axis3(1)]), max([axis1(2),axis2(2),axis3(2)]),...
    min([axis1(3),axis2(3),axis3(3)]), max([axis1(4),axis2(4),axis3(4)])];
%subplot(2,2,1);
%axis(axis0);
%subplot(2,2,2);
%axis(axis0);
%subplot(2,2,3);
%axis(axis0);
clear legend_str
subplot(2,2,3);
hold on;
plot(-binranges, bincountsTT', '-k', 'linewidth', 2);
legend_str{1} = 'TT';
plot(-binranges, bincountsFF, '--', 'color', [.3 .3 .3], 'linewidth', 1.4);
legend_str{2} = 'FF';
plot(-binranges, bincountsTF, ':', 'color', [.5 .5 .5], 'linewidth', 1.4);
legend_str{3} = 'TF';
plot(-binranges, bincountsFT, '-.', 'color', [.7 .7 .7], 'linewidth', 1.4);
legend_str{4} = 'FT';
legend(legend_str);
%}
% Manual perfcurve

% RANSAC error vs. number of iterations
%{
figure;
scatter(all_nit(all_true1nan&all_true2nan), all_RMSE(all_true1nan&all_true2nan), '+b');
hold on
scatter(all_nit(~(all_true1nan&all_true2nan)), all_RMSE(~(all_true1nan&all_true2nan)), '*r');
xlabel('Number of iterations');
ylabel('RMSE [m]');
%}

% RMSE vs sensor distance scatter
%{
figure;
scatter(all_t(all_true1nan&all_true2nan), all_RMSE(all_true1nan&all_true2nan), '+b');
hold on
scatter(all_t(~(all_true1nan&all_true2nan)), all_RMSE(~(all_true1nan&all_true2nan)), '*r');
xlabel('Distance between sensors');
ylabel('RMSE [m]')
%}

% RMSE vs sensor distance of true scatter
%{
figure;
scatter(all_t(all_true1nan&all_true2nan), all_RMSE(all_true1nan&all_true2nan), '+b');
xlabel('Distance between sensors');
ylabel('RMSE [m]')
%}

% Box plot 
%{
%makeTikzBoxplot
%data = [all_t(all_true1nan&all_true2nan) all_RMSE(all_true1nan&all_true2nan)]'; 
%cjrsBoxPlot(data)
figure
all_tT = all_t(all_true1nan&all_true2nan);
all_RMSET = all_RMSE(all_true1nan&all_true2nan);
int = 5;
xrange = 0:int:max(all_tT)+int;
all_ownership = zeros(size(all_tT));
label = cell(numel(xrange)-1);
RMSE_mean = zeros(1,numel(xrange)-1);
for i = 1:numel(xrange)-1;
    ix_int = (all_tT>= xrange(i) & all_tT< xrange(i+1));
    all_ownership(ix_int) = i;
    RMSE_mean(i) = nanmean(all_RMSET(ix_int));
    label{i} = sprintf('%2.0f - %2.0f', xrange(i), xrange(i+1));
end
boxplot(all_RMSET, all_ownership,'widths',[.9 .9 .9 .9 .9]);
hold on
plot(1:numel(xrange)-1, RMSE_mean, '*k');
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
set(gca, 'xticklabel',label );
xlabel('Distance between sensors [m]');
ylabel('RMSE [m]')
% regression
is_valid = ~isnan(all_RMSET);
p = polyfit(all_tT(is_valid), all_RMSET(is_valid),1);
xint = 1:numel(xrange)-1;
xmean = (xrange + circshift(xrange,[0,1]))/2;
xmean = xmean(2:end);
yfit = polyval(p,xmean);
plot(xint,yfit,'--', 'color', [.5 .5 .5], 'linewidth', lw)
clear legend_str
legend_str{1} = 'Mean';
legend_str{2} = sprintf('Fit: RMSE = %3.2f*Range + %3.2f', p(1), p(2));
legend(legend_str);
filepath_rmsevsdist = 'Z:\Desktop\rmsevsdist.tex';
matlab2tikz(filepath_rmsevsdist)
%}

% Count vs sensor distance 
%{
binranges = 0:1:max(all_t);
bincountsTT = histc(all_RMSE(all_true1nan&all_true2nan), binranges);
bincountsF = histc(all_RMSE(~(all_true1nan&all_true2nan)), binranges);
figure
hold on
plot(binranges, bincountsTT, '-', 'color', [.4 .4 .4], 'linewidth', 2);
plot(binranges, bincountsTT, '--k', 'linewidth', 2);
xlabel('Distance between sensors');
ylabel('Count')
%}

% range to each scanner location 
%
n_scan = 25;
vec_plot = 1:25;
vec_x = [10:-5:-10, -10:5:10, 10:-5:-10, -10:5:10, 10:-5:-10];
vec_y = [repmat(-10, [1,5]),repmat(-5, [1,5]),zeros(1,5),repmat(5, [1,5]),repmat(10, [1,5])];
dist_grid = zeros(n_scan, n_scan);
for i = 1:n_scan;
    dist_grid(i,:) = sqrt((vec_x-vec_x(i)).^2 + (vec_y-vec_y(i)).^2);
end
figure;
hold on
axis([-12 12 -12 12]);
info_plot = 1;
scatter(vec_x, vec_y, 500, dist_grid(info_plot,:), 'filled');
for i = 1:n_scan;
    text(vec_x(i), vec_y(i),sprintf('%g', vec_plot(i)), 'fontsize', 14);
end
axis equal
%}

% Percent detection *IN PAPER
%
% Determine expected number of scans
n_scan = 25;
vec_plot = 1:25;
vec_x = [10:-5:-10, -10:5:10, 10:-5:-10, -10:5:10, 10:-5:-10];
vec_y = [repmat(-10, [1,5]),repmat(-5, [1,5]),zeros(1,5),repmat(5, [1,5]),repmat(10, [1,5])];
dist_grid = zeros(n_scan, n_scan);
for i = 1:n_scan;
    dist_grid(i,:) = sqrt((vec_x-vec_x(i)).^2 + (vec_y-vec_y(i)).^2);
end
t_RMSE = 1;
m_det = zeros(n_scan, n_scan);
m_detTF = zeros(n_scan, n_scan);
m_detFT = zeros(n_scan, n_scan);
for i = 1:n_scan;
    for j = 1:n_scan;
        is_valid = (all_plotI ==i) & (all_plotJ == j) & (all_RMSE < t_RMSE);
        is_validTF = (all_plotI ==i) & (all_plotJ == j) & (all_true1nan);
        is_validFT = (all_plotI ==i) & (all_plotJ == j) & (all_true2nan);
        m_det(i,j) = sum(is_valid);
        m_detTF(i,j) = sum(is_validTF);
        m_detFT(i,j) = sum(is_validFT);
    end
end

% collate into unique values
uniq_dist = unique(dist_grid(:));
n_uniq = numel(uniq_dist);
uniq_det = zeros(n_uniq,1);
uniq_detTF = zeros(n_uniq,1);
uniq_detFT = zeros(n_uniq,1);
for u = 1:n_uniq;
    is_u = (dist_grid==uniq_dist(u));
    uniq_det(u) = mean(m_det(is_u));
    uniq_detTF(u) = mean(m_detTF(is_u));
    uniq_detFT(u) = mean(m_detFT(is_u));
end

% Percent detection
%{
figure;
plot(uniq_dist(:), uniq_det(:)./max(m_det(:)), '-ok');
xlabel('Range [m]');
ylabel('Percent Detection');
ylim([0 1]);
%}

% Normalize by integrated circle overlap 
rmax = 16.0377;
res = .1;
npix = 2*rmax/res;
imageSizeX = ceil(npix);
imageSizeY = ceil(npix);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = floor(npix/2);
centerY = floor(npix/2);
radius = floor(rmax./res);
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
%image(circlePixels) ;
%colormap([0 0 0; 1 1 1]);
%title('Binary image of a circle');

% Pad 
I1 = [circlePixels zeros(imageSizeX, imageSizeY)];
I2 = [zeros(imageSizeX, imageSizeY) circlePixels ];
integral = zeros(imageSizeX,1);
for i = 1:2*imageSizeX;
    integral(i) = sum(sum(circshift(I2,[0,-i]) & I1));
end

pixel_step = 1:2*imageSizeX;
m_step = pixel_step*res;
integral_norm = integral./max(integral);
figure; plot(m_step, integral_norm);
[~,ix_max] = max(integral_norm);

m_int = integral_norm(ix_max:end);
m_range = (m_step(ix_max:end)- m_step(ix_max))';
figure; plot(m_range, m_int);

uniq_det_norm = uniq_det(:)./max(m_det(:));
uniq_detTF_norm = uniq_detTF(:)./max(m_detTF(:));
uniq_detFT_norm = uniq_detFT(:)./max(m_detFT(:));

% Percent detection
figure;
plot(uniq_dist, uniq_det_norm, '-ok');
m_int_interp = interp1q(m_range, m_int, uniq_dist);
normnorm_det = uniq_det_norm./m_int_interp; 
hold on
plot(uniq_dist, normnorm_det, '-*k');
plot(m_range, m_int, '--', 'color', [.5 .5 .5]);
xlabel('Range [m]');
ylabel('Percent Detection');
ylim([0 1]);
legend_str{1} = 'Uncorrected';
legend_str{2} = 'Corrected for Reduced Area of Overlap';
legend_str{3} = 'Normalized Area of Overlap for Circle, r = 16 m';
legend(legend_str);

figure;
lw = 1.2;
ms = 10;
hold on
plot(uniq_dist, 100*uniq_detTF_norm, '--*','color', [.5 .5 .5], 'linewidth', lw, 'markersize', ms);
plot(uniq_dist, 100*uniq_detFT_norm, '--o','color', [.5 .5 .5], 'linewidth', lw, 'markersize', ms);
m_int_interp = interp1q(m_range, m_int, uniq_dist);
normnorm_detTF = uniq_detTF_norm./m_int_interp; 
normnorm_detFT = uniq_detFT_norm./m_int_interp; 
plot(uniq_dist, 100*normnorm_detTF, '-*k', 'linewidth', lw, 'markersize', ms);
plot(uniq_dist, 100*normnorm_detFT, '-ok', 'linewidth', lw, 'markersize', ms);
plotyy([],[],m_range, 100*m_int);
xlabel('Range [m]');
ylabel('Percent Detection [%]');
ylim([0 100]);
legend_str{1} = 'Uncorrected Forward';
legend_str{2} = 'Uncorrected Reverse';
legend_str{3} = 'Corrected Forward';
legend_str{4} = 'Corrected Reverse';
legend_str{5} = 'Normalized Area of Circle Overlap';
legend(legend_str);
filepath_pdvsrange = 'Z:\Desktop\pdvsrange.tex';
matlab2tikz(filepath_pdvsrange)
%}

%%
isnnan = ~isnan(all12_tx);
all_true1nan = all_true1;
all_true1nan(isnan(all_true1)) = 0;
all_true1nan = logical(all_true1nan);
all_true2nan = all_true2;
all_true2nan(isnan(all_true2)) = 0;
all_true2nan = logical(all_true2nan);

% txyz vs site scatter
%{
isnnan = ~isnan(all12_tx);
all_true1nan = all_true1;
all_true1nan(isnan(all_true1)) = 0;
all_true1nan = logical(all_true1nan);
all_true2nan = all_true2;
all_true2nan(isnan(all_true2)) = 0;
all_true2nan = logical(all_true2nan);
all12_txyz = sqrt(all12_tx.^2 + all12_ty.^2 + all12_tz.^2);
ixunique = 1:numel(site_unique);
all_ix = all_site;
for s = 1:numel(site_unique);
    all_ix(all_site==site_unique(s)) = ixunique(s);
end
figure;
scatter(all_ix(isnnan&all_true1nan&all_true2nan),...
   all12_txyz(isnnan&all_true1nan&all_true2nan), 10,...
   all_RMSE(isnnan&all_true1nan&all_true2nan), 'filled');
%}

% txyz vs site boxplot 
%{
figure;
boxplot(all12_txyz(isnnan&all_true1nan&all_true2nan),...
    all_ix(isnnan&all_true1nan&all_true2nan),'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('Translation error [m]')
%}

% rxyz vs site scatter
%{
all12_rxyz = sqrt(all12_rx.^2 + all12_ry.^2 + all12_rz.^2);
figure;
scatter(all_site(isnnan&all_true1nan&all_true2nan),...
   all12_rxyz(isnnan&all_true1nan&all_true2nan), 10,...
   all_RMSE(isnnan&all_true1nan&all_true2nan),'filled');
%}

% rxyz vs site boxplot
%{
figure;
boxplot(all12_rxyz(isnnan&all_true1nan&all_true2nan),...
    all_ix(isnnan&all_true1nan&all_true2nan),'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('Rotation error [deg]')
%}

% broken rx, ry, rz vs site boxplot
%{
temp = [all12_rx(isnnan&all_true1nan&all_true2nan) ...
    all12_ry(isnnan&all_true1nan&all_true2nan) ... 
    all12_rz(isnnan&all_true1nan&all_true2nan)];
temp2 = all_ix(isnnan&all_true1nan&all_true2nan);
figure;
boxplot(temp,...
    temp2,'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('Rotation error [deg]')
%}

% rx vs site boxplot
%{
f1= figure;
boxplot(all12_rx(isnnan&all_true1nan&all_true2nan),...
    all_ix(isnnan&all_true1nan&all_true2nan),'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('Rotation x error [deg]')
a1 = gca;
axis1 = axis;
%}

% ry vs site boxplot
%{
f2 = figure;
boxplot(all12_ry(isnnan&all_true1nan&all_true2nan),...
    all_ix(isnnan&all_true1nan&all_true2nan),'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('Rotation y error [deg]')
a2 = gca;
axis2 = axis;
%}

% rz vs site boxplot
%{
f3 = figure;
boxplot(all12_rz(isnnan&all_true1nan&all_true2nan),...
    all_ix(isnnan&all_true1nan&all_true2nan),'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('Rotation z error [deg]')
axis3 = axis;
a3 = gca;
axis0 = [min([axis1(1) axis2(1) axis3(1)]) ...
    max([axis1(2) axis2(2) axis3(2)]) ...
    min([axis1(3) axis2(3) axis3(3)]) ...
    max([axis1(4) axis2(4) axis3(4)])];
figure(f1);
axis(axis0);
figure(f2);
axis(axis0);
figure(f3);
axis(axis0);
%}

% RMSE vs site boxplot  *IN PAPER
%{
figure;
boxplot(all_RMSE(isnnan&all_true1nan&all_true2nan),...
    all_ix(isnnan&all_true1nan&all_true2nan),'widths',[.9 .9 .9 .9 .9]);
%plot(1:numel(xrange), RMSE_mean, '*k');
%plot([0 numel(xrange)+1],[.3464 .3464], ':', 'color', [.5 .5 .5], 'linewidth', lw)
%plot([1 numel(xrange)],1.6*[.2 .2], ':', 'color', [.5 .5 .5], 'linewidth', lw)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
h = findobj('Tag','Box');
set(h,'Color',[0,0,0]);
h = findobj(gca,'Tag','Outliers');
set(h,'MarkerEdgeColor',[.5,.5,.5]);
h = findobj('Tag','Median');
set(h,'Color',[.5,.5,.5]);
set(findobj(gca,'type','line'),'linew',1.5)
xlabel('Site');
set(gca, 'xticklabel',  site_unique)
ylabel('RMSE error [m]')
filepath_rmsevssite = 'Z:\Desktop\rmsevssite.tex';
matlab2tikz(filepath_rmsevssite)
%}

figure;
scatter(all12_txyz(isnnan&all_true1nan&all_true2nan),...
   all_RMSE(isnnan&all_true1nan&all_true2nan), 10, 'k', 'filled');
xlabel('Translation error [m]');
ylabel('RMSE [m]');
figure;
hold on
rall = [abs(all12_rx(isnnan&all_true1nan&all_true2nan)),...
   abs(all12_ry(isnnan&all_true1nan&all_true2nan)),...
   abs(all12_rz(isnnan&all_true1nan&all_true2nan))];
rmax = max(rall, [], 2);
scatter(rmax, all_RMSE(isnnan&all_true1nan&all_true2nan), 10, 'k', 'filled');
xlabel('Maximum rotation error [deg]');
ylabel('RMSE [m]');

fprintf('Mean tx = %3.3f m\n', mean(abs(all12_tx(isnnan&all_true1nan&all_true2nan))));
fprintf('Mean ty = %3.3f m\n', mean(abs(all12_ty(isnnan&all_true1nan&all_true2nan))));
fprintf('Mean tz = %3.3f m\n', mean(abs(all12_tz(isnnan&all_true1nan&all_true2nan))));
fprintf('Mean rx = %3.3f deg\n', mean(abs(all12_rx(isnnan&all_true1nan&all_true2nan))));
fprintf('Mean ry = %3.3f deg\n', mean(abs(all12_ry(isnnan&all_true1nan&all_true2nan))));
fprintf('Mean rz = %3.3f deg\n', mean(abs(all12_rz(isnnan&all_true1nan&all_true2nan))));
fprintf('Mean txyz = %3.3f m\n', mean(abs(all12_txyz(isnnan&all_true1nan&all_true2nan))));
