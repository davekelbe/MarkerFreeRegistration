% batch_graph_analyses
options_showfig = true;
%{
% GraphAnalRMSE2
%aux.info_exp = 'GraphAnalRMSE2';
%sigma_mean = 0:.025:.25;
%n_trial = 10;
% GraphAnalRMSE3
%aux.info_exp = 'GraphAnalRMSE3';
%sigma_mean = .025:.025:.05;
%n_trial = 1;
% GraphAnalRMSE3
%aux.info_exp = 'GraphAnalRMSE4';
%sigma_mean = .025:.025:.05;
%n_trial = 1;
% GraphAnalRMSE5
%aux.info_exp = 'GraphAnalRMSE5';
%sigma_mean = .025:.025:.025;
%n_trial = 1;
%aux.info_exp = 'GraphAnalRMSE6';
%sigma_mean = 0:.001:.05;
%n_trial = 1;
%aux.info_exp = 'GraphAnalRMSE7';
%sigma_mean = 0:.001:.005;
%n_trial = 1;
%aux.info_exp = 'GraphAnalRMSE8';
%sigma_mean = 0:.01:.05;
%n_trial = 1;
%aux.info_exp = 'GraphAnalRMSE9';
%sigma_mean = 0:.01:.05;
%n_trial = 1;
%aux.info_exp = 'GraphAnalRMSE10';
%sigma_mean = 0:.001:.05;
%n_trial = 1;
%aux.info_exp = 'GraphAnalRMSE11';
%sigma_mean = 0:.001:.1;
%n_trial = 1;
aux.info_exp = 'GraphAnalRMSE12';
sigma_mean = 0:.001:.1;
n_trial = 1;
aux.info_exp = 'GraphAnalRMSE12Copy';
sigma_mean = 0:.001:.1;
n_trial = 1;
%}
aux.info_exp = 'GraphAnalRMSE13';
sigma_mean = .1:.1:.5;
n_trial = 30;
n_step = numel(sigma_mean);
for step = 3:3;%n_step;
    fprintf('\n Step %g of %g\n',step, n_step);
    aux.info_suffix = sprintf('step-%02.0f',step);
    aux.S_n = sigma_mean(step);
    for trial = 1:n_trial;
        fprintf('\n Trial %g of %g\n',trial, n_trial);
        aux.trial = sprintf('%02.0f',trial);
        graph_analysis_rmsefun_generic( aux )
    end
end
return
path_fig = sprintf('%s%s%s','Z:\Desktop\KelbeGraph\',aux.info_exp,'\');
if ~exist(path_fig, 'dir');
    mkdir(path_fig);
end
path_tikz = sprintf('%s%s%s','Z:\Documents\Research\Documents\KelbeThesis\thesis\');

%% Load variables and collate
n_S = 25;
all_P_RMSEin = nan(n_step,n_trial,n_S);
all_RMSE_pair = nan(n_step,n_trial,n_S);
all_RMSE_pairs = nan(n_step,n_trial,n_S,n_S);
all_RMSE_MST = nan(n_step,n_trial,n_S);
all_RMSE_sum_MST = nan(n_step,n_trial,n_S);
all_RMSE_geo_MST = nan(n_step,n_trial,n_S);
all_RMSE_mult_MST = nan(n_step,n_trial,n_S);
all_RMSE_max_MST = nan(n_step,n_trial,n_S);
all_RMSEh_sum_MST = nan(n_step,n_trial,n_S);
all_RMSEh_geo_MST = nan(n_step,n_trial,n_S);
all_RMSEh_mult_MST = nan(n_step,n_trial,n_S);
all_RMSEh_max_MST = nan(n_step,n_trial,n_S);
all_RMSE_path_MST = nan(n_step,n_trial,n_S);
all_RMSE_WMF = nan(n_step,n_trial,n_S);
all_RMSE_SVD = nan(n_step,n_trial,n_S);
all_R_pair = cell(n_step,n_trial);
all_t_pair = cell(n_step,n_trial);
all_R12_pair = cell(n_step,n_trial);
all_t12_pair = cell(n_step,n_trial);
all_R_pairs = cell(n_step,n_trial);
all_t_pairs = cell(n_step,n_trial);
all_R_MST = cell(n_step,n_trial);
all_t_MST = cell(n_step,n_trial);
all_R_WMF = cell(n_step,n_trial);
all_t_WMF = cell(n_step,n_trial);
all_R_SVD = cell(n_step,n_trial);
all_t_SVD = cell(n_step,n_trial);
all_P_LCS = cell(n_step,n_trial);
all_P_LCSn = cell(n_step,n_trial);
all_P_radn = cell(n_step,n_trial);
all_P_ix = cell(n_step,n_trial);
all_P_WCSn = cell(n_step,n_trial);
all_P_WCS = cell(n_step,n_trial);
all_S_R = cell(n_step,n_trial);
all_S_t = cell(n_step,n_trial);
all_S_range = zeros(n_step,n_trial,n_S);
all_w = zeros(n_step,n_trial);
all_T_x = cell(n_step,n_trial);
all_T_y = cell(n_step,n_trial);
all_T_z = cell(n_step,n_trial);
all_T_r = cell(n_step,n_trial);
all_step = zeros(n_step,n_trial,n_S);
all_step_pairs = zeros(n_step,n_trial,n_S,n_S);
for step = 1:n_step;
    aux.info_suffix = sprintf('step-%02.0f',step);
    for trial = 1:n_trial
        all_step(step,trial,:) = repmat(step, [1,n_S]);
        all_step_pairs(step,trial,:,:) = repmat(step, [n_S,n_S]);
        aux.trial = sprintf('%02.0f',trial);
        path_mat = sprintf('%s%s%s%s%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix,'\',aux.trial,'\');
        % Truth values
        filepath_T_x = sprintf('%s%s',path_mat, 'T_x.mat');
        load(filepath_T_x);
        all_T_x{step,trial}= T_x;
        filepath_T_y = sprintf('%s%s',path_mat, 'T_y.mat');
        load(filepath_T_y);
        all_T_y{step,trial}= T_y;
        filepath_T_z = sprintf('%s%s',path_mat, 'T_z.mat');
        load(filepath_T_z);
        all_T_z{step,trial}= T_z;
        filepath_T_r = sprintf('%s%s',path_mat, 'T_r.mat');
        load(filepath_T_r);
        all_T_r{step,trial}= T_r;
        % RMSE in
        filepath_P_RMSEin = sprintf('%s%s',path_mat, 'P_RMSEin.mat');
        load(filepath_P_RMSEin);
        all_P_RMSEin(step,trial,:)= P_RMSEin;
        % P_LCS
        filepath_P_LCS = sprintf('%s%s',path_mat, 'P_LCS.mat');
        load(filepath_P_LCS);
        all_P_LCS{step,trial}= P_LCS;
        % P_LCSn
        filepath_P_LCSn = sprintf('%s%s',path_mat, 'P_LCSn.mat');
        load(filepath_P_LCSn);
        all_P_LCSn{step,trial}= P_LCSn;
        % P_radn
        filepath_P_radn = sprintf('%s%s',path_mat, 'P_radn.mat');
        load(filepath_P_radn);
        all_P_radn{step,trial}= P_radn;
        % P_ix
        filepath_P_ix = sprintf('%s%s',path_mat, 'P_ix.mat');
        load(filepath_P_ix);
        all_P_ix{step,trial}= P_ix;
        % P_WCS
        filepath_P_WCS = sprintf('%s%s',path_mat, 'P_WCS.mat');
        load(filepath_P_WCS);
        all_P_WCS{step,trial}= P_WCS;
        % P_WCSn
        filepath_P_WCSn = sprintf('%s%s',path_mat, 'P_WCSn.mat');
        load(filepath_P_WCSn);
        all_P_WCSn{step,trial}= P_WCSn;
        % Reference node
        filepath_w = sprintf('%s%s',path_mat, 'w.mat');
        load(filepath_w);
        all_w(step,trial)= w;
        % Sensor poses
        filepath_S_R = sprintf('%s%s',path_mat, 'S_R.mat');
        load(filepath_S_R);
        all_S_R{step,trial} = S_R;
        filepath_S_t = sprintf('%s%s',path_mat, 'S_t.mat');
        load(filepath_S_t);
        all_S_t{step,trial} = S_t;
        % Distance from each j to w
        S_range = cellfun(@(x) sqrt(sum((x-S_t{w}).^2)), S_t);
        all_S_range(step,trial,:) = S_range;
        % Pair
        filepath_match_R = sprintf('%s%s',path_mat, 'match_R.mat');
        load(filepath_match_R);
        all_R_pair{step,trial} = match_R(w,:)';
        filepath_match_t = sprintf('%s%s',path_mat, 'match_t.mat');
        load(filepath_match_t);
        all_t_pair{step,trial} = match_t(w,:)';
        filepath_match12_RMSE = sprintf('%s%s',path_mat, 'match12_RMSE.mat');
        load(filepath_match12_RMSE);
        all_RMSE_pair(step,trial,:) = match12_RMSE(w,:)';
        filepath_match12_R = sprintf('%s%s',path_mat, 'match12_R.mat');
        load(filepath_match12_R);
        all_R12_pair{step,trial} = match12_R(w,:)';
        filepath_match12_t = sprintf('%s%s',path_mat, 'match12_t.mat');
        load(filepath_match12_t);
        all_t12_pair{step,trial} = match12_t(w,:)';
        % Pairs
        all_R_pairs{step,trial} = match_R;
        all_t_pairs{step,trial} = match_t;
        all_RMSE_pairs(step,trial,:,:) = match12_RMSE;
        % MST
        filepath_G_R_MST = sprintf('%s%s',path_mat, 'G_R_MST.mat');
        load(filepath_G_R_MST);
        all_R_MST{step,trial} = G_R_MST;
        filepath_G_t_MST = sprintf('%s%s',path_mat, 'G_t_MST.mat');
        load(filepath_G_t_MST);
        all_t_MST{step,trial} = G_t_MST;
        filepath_G_RMSE_MST = sprintf('%s%s',path_mat, 'G_RMSE_MST.mat');
        load(filepath_G_RMSE_MST);
        filepath_G_RMSE_sum_MST = sprintf('%s%s',path_mat, 'G_RMSE_sum_MST.mat');
        load(filepath_G_RMSE_sum_MST);
        all_RMSE_sum_MST(step,trial,:)= G_RMSE_sum_MST;
        filepath_G_RMSE_geo_MST = sprintf('%s%s',path_mat, 'G_RMSE_geo_MST.mat');
        load(filepath_G_RMSE_geo_MST);
        all_RMSE_geo_MST(step,trial,:)= G_RMSE_geo_MST;
        filepath_G_RMSE_mult_MST = sprintf('%s%s',path_mat, 'G_RMSE_mult_MST.mat');
        load(filepath_G_RMSE_mult_MST);
        all_RMSE_mult_MST(step,trial,:)= G_RMSE_mult_MST;
        filepath_G_RMSE_max_MST = sprintf('%s%s',path_mat, 'G_RMSE_max_MST.mat');
        load(filepath_G_RMSE_max_MST);
        all_RMSE_max_MST(step,trial,:)= G_RMSE_max_MST;
        filepath_G_RMSE_path_MST = sprintf('%s%s',path_mat, 'G_RMSE_path_MST.mat');
        load(filepath_G_RMSE_path_MST);
        all_RMSE_path_MST(step,trial,:)= G_RMSE_path_MST;
        % half 
        filepath_G_RMSEh_sum_MST = sprintf('%s%s',path_mat, 'G_RMSEh_sum_MST.mat');
        load(filepath_G_RMSEh_sum_MST);
        all_RMSEh_sum_MST(step,trial,:)= G_RMSEh_sum_MST;
        filepath_G_RMSEh_geo_MST = sprintf('%s%s',path_mat, 'G_RMSEh_geo_MST.mat');
        load(filepath_G_RMSEh_geo_MST);
        all_RMSEh_geo_MST(step,trial,:)= G_RMSEh_geo_MST;
        filepath_G_RMSEh_mult_MST = sprintf('%s%s',path_mat, 'G_RMSEh_mult_MST.mat');
        load(filepath_G_RMSEh_mult_MST);
        all_RMSEh_mult_MST(step,trial,:)= G_RMSEh_mult_MST;
        filepath_G_RMSEh_max_MST = sprintf('%s%s',path_mat, 'G_RMSEh_max_MST.mat');
        load(filepath_G_RMSEh_max_MST);
        all_RMSEh_max_MST(step,trial,:)= G_RMSEh_max_MST;
        % WMF
        filepath_G_R_WMF = sprintf('%s%s',path_mat, 'G_R_WMF.mat');
        load(filepath_G_R_WMF);
        all_R_WMF{step,trial} = G_R_WMF;
        filepath_G_t_WMF = sprintf('%s%s',path_mat, 'G_t_WMF.mat');
        load(filepath_G_t_WMF);
        all_t_WMF{step,trial} = G_t_WMF;
        filepath_G_RMSE_WMF = sprintf('%s%s',path_mat, 'G_RMSE_WMF.mat');
        load(filepath_G_RMSE_WMF);
        all_RMSE_WMF(step,trial,:)= G_RMSE_WMF;
        % SVD
        filepath_G_R_SVD= sprintf('%s%s',path_mat, 'G_R_SVD.mat');
        load(filepath_G_R_SVD);
        all_R_SVD{step,trial} = G_R_SVD;
        filepath_G_t_SVD = sprintf('%s%s',path_mat, 'G_t_SVD.mat');
        load(filepath_G_t_SVD);
        all_t_SVD{step,trial} = G_t_SVD;
        filepath_G_RMSE_SVD = sprintf('%s%s',path_mat, 'G_RMSE_SVD.mat');
        load(filepath_G_RMSE_SVD);
        all_RMSE_SVD(step,trial,:)= G_RMSE_SVD;
    end
end 

clear filepath* step trial w
%% Transform points into different coordinate systems
isvalid_pair= zeros(n_step,n_trial,n_S);
isvalid_pairs= zeros(n_step,n_trial,n_S,n_S);
isvalid_MST= zeros(n_step,n_trial,n_S);
isvalid_12= zeros(n_step,n_trial,n_S);
isvalid_WMF= zeros(n_step,n_trial,n_S);
isvalid_SVD= zeros(n_step,n_trial,n_S);
for step = 1:n_step;
    for trial = 1:n_trial;
        isvalid_pair(step,trial,:)= ~cellfun(@isempty, all_R_pair{step,trial});
        isvalid_pairs(step,trial,:,:)= ~cellfun(@isempty, all_R_pairs{step,trial});
        isvalid_12(step,trial,:)= ~cellfun(@isempty, all_R12_pair{step,trial});
        isvalid_MST(step,trial,:)= ~cellfun(@isempty, all_R_MST{step,trial});
        isvalid_WMF(step,trial,:)= ~cellfun(@isempty, all_R_WMF{step,trial});
        isvalid_SVD(step,trial,:)= ~cellfun(@isempty, all_R_SVD{step,trial});
    end
end

all_P_LCSnw_pair = cell(n_step, n_trial, n_S);
all_P_LCSni_pairs = cell(n_step, n_trial, n_S,n_S);
all_P_LCSnj12 = cell(n_step, n_trial, n_S);
all_P_LCSnW12 = cell(n_step, n_trial, n_S);
all_P_LCSnw_MST = cell(n_step, n_trial, n_S);
all_P_LCSnw_WMF = cell(n_step, n_trial, n_S);
all_P_LCSnw_SVD = cell(n_step, n_trial, n_S);
all_P_LCSnW_pair = cell(n_step, n_trial, n_S);
all_P_LCSnW_MST = cell(n_step, n_trial, n_S);
all_P_LCSnW_WMF = cell(n_step, n_trial, n_S);
all_P_LCSnW_SVD = cell(n_step, n_trial, n_S);
all_P_WCSnw = cell(n_step, n_trial, n_S);
all_P_WCSni = cell(n_step, n_trial, n_S,n_S);
all_P_WCSi = cell(n_step, n_trial, n_S,n_S);
all_P_ntree = zeros(n_step, n_trial,n_S);
for step = 1:n_step;
    for trial = 1:n_trial
        for j = 1:n_S;
            [~,P_ntree] = cellfun(@size, all_P_LCSn{step, trial});
            all_P_ntree(step,trial,:) = P_ntree;
            w = all_w(step,trial);
            if isvalid_12(step,trial,j);
                all_P_LCSnj12{step,trial,j} = all_R12_pair{step, trial}{j}*all_P_LCSn{step,trial}{j} + ...
                    repmat(all_t12_pair{step, trial}{j}, [1,P_ntree(j)]);
                all_P_LCSnW12{step,trial,j}= all_S_R{step, trial}{j}*all_R12_pair{step, trial}{j}*...
                    all_P_LCSn{step,trial}{j} + ...
                    repmat(...
                    all_S_R{step, trial}{j}*all_t12_pair{step, trial}{j} + ...
                    all_S_t{step, trial}{j},...
                    [1 P_ntree(j)]);
            end
            if isvalid_MST(step,trial,j);
                all_P_LCSnw_MST{step,trial,j} = all_R_MST{step, trial}{j}*all_P_LCSn{step,trial}{j} + ...
                    repmat(all_t_MST{step, trial}{j}, [1,P_ntree(j)]);
                all_P_LCSnW_MST{step,trial,j}= all_S_R{step, trial}{w}*all_R_MST{step, trial}{j}*...
                    all_P_LCSn{step,trial}{j} + ...
                    repmat(...
                    all_S_R{step, trial}{w}*all_t_MST{step, trial}{j} + ...
                    all_S_t{step, trial}{w},...
                    [1 P_ntree(j)]);
            end
            if isvalid_WMF(step,trial,j);
                all_P_LCSnw_WMF{step,trial,j} = all_R_WMF{step, trial}{j}*all_P_LCSn{step,trial}{j} + ...
                    repmat(all_t_WMF{step, trial}{j}, [1,P_ntree(j)]);
                all_P_LCSnW_WMF{step,trial,j}= all_S_R{step, trial}{w}*all_R_WMF{step, trial}{j}*...
                    all_P_LCSn{step,trial}{j} + ...
                    repmat(...
                    all_S_R{step, trial}{w}*all_t_WMF{step, trial}{j} + ...
                    all_S_t{step, trial}{w},...
                    [1 P_ntree(j)]);
            end
            if isvalid_SVD(step,trial,j);
                all_P_LCSnw_SVD{step,trial,j} = all_R_SVD{step, trial}{j}*all_P_LCSn{step,trial}{j} + ...
                    repmat(all_t_SVD{step, trial}{j}, [1,P_ntree(j)]);
                all_P_LCSnW_SVD{step,trial,j}= all_S_R{step, trial}{w}*all_R_SVD{step, trial}{j}*...
                    all_P_LCSn{step,trial}{j} + ...
                    repmat(...
                    all_S_R{step, trial}{w}*all_t_SVD{step, trial}{j} + ...
                    all_S_t{step, trial}{w},...
                    [1 P_ntree(j)]);
            end
            % Pairwise
            if isvalid_pair(step,trial,j)
                all_P_LCSnw_pair{step,trial,j} = all_R_pair{step, trial}{j}*...
                    all_P_LCSn{step,trial}{j} + ...
                    repmat(...
                    all_t_pair{step, trial}{j},[1 P_ntree(j)]);
                all_P_LCSnW_pair{step,trial,j} = all_S_R{step, trial}{w}*all_R_pair{step, trial}{j}*...
                    all_P_LCSn{step,trial}{j} + ...
                    repmat(...
                    all_S_R{step, trial}{w}*all_t_pair{step, trial}{j} + ...
                    all_S_t{step, trial}{w},...
                    [1 P_ntree(j)]);
            end
            % WCS into w
            all_P_WCSnw{step,trial,j} = all_S_R{step, trial}{w}'*all_P_WCSn{step,trial}{j}-...
                repmat(all_S_R{step, trial}{w}'*all_S_t{step, trial}{w}, [1, P_ntree(j)]);
            % WCS into i
            for i = 1:n_S;
                all_P_WCSni{step,trial,i,j} = all_S_R{step, trial}{i}'*all_P_WCSn{step,trial}{j}-...
                    repmat(all_S_R{step, trial}{i}'*all_S_t{step, trial}{i}, [1, P_ntree(j)]);
                all_P_WCSi{step,trial,i,j} = all_S_R{step, trial}{i}'*all_P_WCS{step,trial}{j}-...
                    repmat(all_S_R{step, trial}{i}'*all_S_t{step, trial}{i}, [1, P_ntree(j)]);
            end
            
            % Pairs
            for i = 1:n_S;
                if isvalid_pairs(step,trial,i,j)
                    all_P_LCSni_pairs{step,trial,i,j} = all_R_pairs{step, trial}{i,j}*...
                        all_P_LCSn{step,trial}{j} + ...
                        repmat(...
                        all_t_pairs{step, trial}{i,j},[1 P_ntree(j)]);
                end
            end
        end
    end
end
%% Sanity check
%{
step = 1;
trial = 1;
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
figure;
hold on
axis equal
for j =1:n_S;
    % MST
    if isvalid_MST(step,trial,j)
        x_axist = all_R_MST{step, trial}{j}*x_axis + repmat(all_t_MST{step, trial}{j},[1,2]);
        y_axist = all_R_MST{step, trial}{j}*y_axis + repmat(all_t_MST{step, trial}{j},[1,2]);
        z_axist = all_R_MST{step, trial}{j}*z_axis + repmat(all_t_MST{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'r', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'r', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'r', 'linewidth',2)
    end
    % WMF
    if isvalid_WMF(step,trial,j)
        x_axist = all_R_WMF{step, trial}{j}*x_axis + repmat(all_t_WMF{step, trial}{j},[1,2]);
        y_axist = all_R_WMF{step, trial}{j}*y_axis + repmat(all_t_WMF{step, trial}{j},[1,2]);
        z_axist = all_R_WMF{step, trial}{j}*z_axis + repmat(all_t_WMF{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'g', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'g', 'linewidth',2)
    end
    % SVD
    if isvalid_SVD(step,trial,j)
        x_axist = all_R_SVD{step, trial}{j}*x_axis + repmat(all_t_SVD{step, trial}{j},[1,2]);
        y_axist = all_R_SVD{step, trial}{j}*y_axis + repmat(all_t_SVD{step, trial}{j},[1,2]);
        z_axist = all_R_SVD{step, trial}{j}*z_axis + repmat(all_t_SVD{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'b', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'b', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'b', 'linewidth',2)
    end
    % Pair
    if isvalid_pair(step,trial,j)
        x_axist = all_R_pair{step, trial}{j}*x_axis + repmat(all_t_pair{step, trial}{j},[1,2]);
        y_axist = all_R_pair{step, trial}{j}*y_axis + repmat(all_t_pair{step, trial}{j},[1,2]);
        z_axist = all_R_pair{step, trial}{j}*z_axis + repmat(all_t_pair{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'c', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'c', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'c', 'linewidth',2)
    end
end
%}
%% Tree-wise error 
%{
true_tRMSE_MST = cell(n_step,n_trial,n_S);
all_tRMSE_MST = cell(n_step,n_trial,n_S);
for step = 1:n_step;
    for trial = 1:n_trial;
        for j = 1:n_S;
            if ~isempty(all_P_LCSnw_MST{step,trial,j});
            all_tRMSE_MST{step, trial, j} = sum((all_P_LCSnj12{step,trial,j} - all_P_LCSn{step,trial}{j}).^2);        
            end            
            if ~isempty(all_P_LCSnw_MST{step,trial,j});
            true_tRMSE_MST{step, trial, j} = sum((all_P_LCSnw_MST{step,trial,j} - all_P_WCSnw{step,trial,j}).^2);        
            end      
        end
    end
end
%}
%% Compute true RMSE values

true_RMSE_MST = nan(n_step, n_trial, n_S);
true_RMSE_WMF = nan(n_step, n_trial, n_S);
true_RMSE_SVD = nan(n_step, n_trial, n_S);
true_RMSE_pair = nan(n_step, n_trial, n_S);
true_RMSE_pair2 = nan(n_step, n_trial, n_S);
true_RMSEn_pairs = nan(n_step, n_trial, n_S,n_S);
true_RMSE_pairs = nan(n_step, n_trial, n_S,n_S);
all_P_ntree = nan(n_step, n_trial,n_S);
for step = 1:n_step;
    for trial = 1:n_trial
        for j = 1:n_S;
            P_ntree = all_P_ntree(step,trial,:);
            w = all_w(step,trial);
            % Sanity check
            %{
            % WCS_j vs T_x (truth)
            %{
            figure;
            hold on
            for t = 1:numel(all_T_x{step,trial});
                h1 = filledCircle([all_T_x{step,trial}(t);...
                    all_T_y{step,trial}(t)]',all_T_r{step,trial}(t),1000,'b');
                set(h1, 'facealpha', 0.5);
            end
            for t = 1:numel(P_rad{j});
                h1 = filledCircle([all_P_WCSn{step,trial}{j}(1,t);...
                    all_P_WCSn{step,trial}{j}(2,t)]',all_P_radn{step,trial}{j}(t),1000,'r');
                set(h1, 'facealpha', 0.5);
            end
            %}
            % P_LCS w reference node vs T_x
            %{
            temp = all_S_R{step, trial}{w}*all_P_LCSn{step,trial}{w} + ...
                repmat(all_S_t{step, trial}{w}, [1,P_ntree(w)]);
            figure;
            hold on
            for t = 1:numel(all_T_x{step,trial});
                h1 = filledCircle([all_T_x{step,trial}(t);...
                    all_T_y{step,trial}(t)]',all_T_r{step,trial}(t),1000,'b');
                set(h1, 'facealpha', 0.5);
            end
            for t = 1:P_ntree(w);
                h1 = filledCircle(temp(1:2,t),all_P_radn{step,trial}{w}(t),1000,'r');
                set(h1, 'facealpha', 0.5);
            end
            %}
            % P_LCS j into P_LCS w
            %{
            temp = all_R_MST{step, trial}{j}*all_P_LCSn{step,trial}{j} + ...
                repmat(all_t_MST{step, trial}{j}, [1,P_ntree(j)]);
            temp_all = [all_T_x{step,trial} all_T_y{step,trial} all_T_z{step,trial}]';
            temp2 = all_S_R{step, trial}{w}'*temp_all-...
                repmat(all_S_R{step, trial}{w}'*all_S_t{step, trial}{w}, [1, size(temp_all,2)]);
            figure;
            hold on
            for t = 1:P_ntree(w);
                h1 = filledCircle(all_P_LCSn{step,trial}{w}(1:2,t),...
                    all_P_radn{step,trial}{w}(t),1000,'b');
                set(h1, 'facealpha', 0.5);
            end
            for t = 1:P_ntree(j);
                h1 = filledCircle(temp(1:2,t),all_P_radn{step,trial}{j}(t),1000,'r');
                set(h1, 'facealpha', 0.5);
            end
            for t = 1:numel(all_T_r{step,trial});
                h1 = filledCircle(temp2(1:2,t),all_T_r{step,trial}(t),1000,'g');
                set(h1, 'facealpha', 0.2);
            end
            %}
            % P_LCS into WCS via w node
            %{
            temp1 = all_R_MST{step, trial}{j}*all_P_LCSn{step,trial}{j} + ...
                repmat(all_t_MST{step, trial}{j}, [1,P_ntree(j)]);
            temp2 = all_S_R{step, trial}{w}*temp1 + ...
                repmat(all_S_t{step, trial}{w}, [1,P_ntree(w)]);
            %}
            % Direct route:
            %{
            PLCSW = all_S_R{step, trial}{w}*all_R_MST{step, trial}{j}*...
                all_P_LCSn{step,trial}{j} + ...
                repmat(...
                all_S_R{step, trial}{w}*all_t_MST{step, trial}{j} + ...
                all_S_t{step, trial}{w},...
                [1 P_ntree(j)]);
            % Figure
            figure;
            hold on
            for t = 1:numel(all_T_x{step,trial});
                h1 = filledCircle([all_T_x{step,trial}(t);...
                    all_T_y{step,trial}(t)]',all_T_r{step,trial}(t),1000,'b');
                set(h1, 'facealpha', 0.5);
            end
            for t = 1:P_ntree(j);
                h1 = filledCircle(temp2(1:2,t),all_P_radn{step,trial}{j}(t),1000,'r');
                set(h1, 'facealpha', 0.5);
            end
            for t = 1:P_ntree(j);
                h1 = filledCircle(PLCSW(1:2,t),all_P_radn{step,trial}{j}(t),1000,'g');
                set(h1, 'facealpha', 0.2);
            end
            %}
            %}
            % MST
            if isvalid_MST(step,trial,j);
                true_RMSE_MST(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_MST{step,trial,j}).^2)));
            end
            % WMF
            if isvalid_WMF(step,trial,j);
                true_RMSE_WMF(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_WMF{step,trial,j}).^2)));
            end
            % SVD
            if isvalid_SVD(step,trial,j);
                true_RMSE_SVD(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_SVD{step,trial,j}).^2)));
            end
            % True pair
            if isvalid_pair(step,trial,j);
                true_RMSE_pair(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSnw{step,trial,j} -...
                    all_P_LCSnw_pair{step,trial,j}).^2)));
                true_RMSE_pair2(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_pair{step,trial,j}).^2)));
            end
            % True pairs
            for i = 1:n_S;
                if isvalid_pairs(step,trial,i,j);
                    true_RMSEn_pairs(step,trial,i,j) = ...
                        sqrt(nanmean(sum((all_P_WCSni{step,trial,i,j} -...
                        all_P_LCSni_pairs{step,trial,i,j}).^2)));
                    true_RMSE_pairs(step,trial,i,j) = ...
                        sqrt(nanmean(sum((all_P_WCSi{step,trial,i,j} -...
                        all_P_LCSni_pairs{step,trial,i,j}).^2)));
                end
            end
        end
        
    end
end
%% Sanity check
%{
err_vec = reshape(true_RMSE_pair, [n_step*n_trial*n_S,1]);
err_color_vec = vec2cmap(err_vec, 'jet' );
err_color = reshape(err_color_vec, [n_step,n_trial,n_S,3]);
[x_sph, y_sph, z_sph] = sphere;
step = 60;
trial = 1;
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
figure;
hold on
axis equal
for j =1:n_S;
    % MST
    if isvalid_MST(step,trial,j)
        x_axist = all_R_MST{step, trial}{j}*x_axis + repmat(all_t_MST{step, trial}{j},[1,2]);
        y_axist = all_R_MST{step, trial}{j}*y_axis + repmat(all_t_MST{step, trial}{j},[1,2]);
        z_axist = all_R_MST{step, trial}{j}*z_axis + repmat(all_t_MST{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'r', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'r', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'r', 'linewidth',2)
    end
    % WMF
    if isvalid_WMF(step,trial,j)
        x_axist = all_R_WMF{step, trial}{j}*x_axis + repmat(all_t_WMF{step, trial}{j},[1,2]);
        y_axist = all_R_WMF{step, trial}{j}*y_axis + repmat(all_t_WMF{step, trial}{j},[1,2]);
        z_axist = all_R_WMF{step, trial}{j}*z_axis + repmat(all_t_WMF{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'g', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'g', 'linewidth',2)
    end
    % SVD
    if isvalid_SVD(step,trial,j)
        x_axist = all_R_SVD{step, trial}{j}*x_axis + repmat(all_t_SVD{step, trial}{j},[1,2]);
        y_axist = all_R_SVD{step, trial}{j}*y_axis + repmat(all_t_SVD{step, trial}{j},[1,2]);
        z_axist = all_R_SVD{step, trial}{j}*z_axis + repmat(all_t_SVD{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'b', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'b', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'b', 'linewidth',2)
    end
    % Pair
    if isvalid_pair(step,trial,j)
        x_axist = all_R_pair{step, trial}{j}*x_axis + repmat(all_t_pair{step, trial}{j},[1,2]);
        y_axist = all_R_pair{step, trial}{j}*y_axis + repmat(all_t_pair{step, trial}{j},[1,2]);
        z_axist = all_R_pair{step, trial}{j}*z_axis + repmat(all_t_pair{step, trial}{j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'color',...
            'c', 'linewidth',2);
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'color',...
            'c', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'color',...
            'c', 'linewidth',2)
        h = surf(x_sph + all_t_pair{step, trial}{j}(1),...
            y_sph + all_t_pair{step, trial}{j}(2), ...
            z_sph + all_t_pair{step, trial}{j}(3));
        set(h, 'facecolor',squeeze(err_color(step,trial,j,:)), 'facealpha', 0.2,...
            'EdgeColor','none','LineStyle','none')
    end
end
%}
%% Pair RMSE out vs RMSE true 
    figure;
    hold on
    isvalid = ~isnan(true_RMSE_pair) & ~isnan(all_RMSE_pair);
    scatter(true_RMSE_pair(isvalid),all_RMSE_pair(isvalid),5, all_step(isvalid), 'filled');
    xlabel('True Pairwise RMSE [m]');
    ylabel('Reported Pairwie RMSE [m]');
    %%
    figure;
    hold on
        isvalid = ~isnan(true_RMSEn_pairs) & ~isnan(all_RMSE_pairs);
    scatter(true_RMSEn_pairs(isvalid),all_RMSE_pairs(isvalid),5, 'r', 'filled');
    isvalid = ~isnan(true_RMSE_pairs) & ~isnan(all_RMSE_pairs);
    scatter(true_RMSE_pairs(isvalid),all_RMSE_pairs(isvalid),5, 'b', 'filled');

    xlabel('True Pairwise RMSE [m]');
    ylabel('Reported Pairwie RMSE [m]');
    %%
    figure;
    hold on
     xlabel('True Pairwise RMSE [m]');
    ylabel('Reported Pairwie RMSE [m]');
%% Pairs (all) RMSE out vs RMSE true
%
if options_showfig;
    fig_str = 'graph-anal-PairRMSEoutvsRMSEtrue';
    filepath_png = sprintf('%s%s.png',path_fig, fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz, fig_str);
    figure;
    hold on
    %title('Is reported pairwise RMSE accurate predictor of true RMSE?')
    isvalid = ~isnan(true_RMSEn_pairs) & ~isnan(all_RMSE_pairs) & ...
        true_RMSEn_pairs <1 & all_RMSE_pairs <1 ;
    plot(true_RMSEn_pairs(isvalid),all_RMSE_pairs(isvalid),'.');
    isvalid = ~isnan(true_RMSEn_pairs) & ~isnan(all_RMSE_pairs);
    scatter(true_RMSEn_pairs(isvalid),all_RMSE_pairs(isvalid),5, all_step_pairs(isvalid), 'filled');

    plot(true_RMSEn_pairs(isvalid),all_RMSE_pairs(isvalid),'.');
    % scatter(true_RMSE_pairs(isvalid),all_RMSE_pairs(isvalid),'+c');
    xlabel('True Pairwise RMSE [m]');
    ylabel('Reported Pairwie RMSE [m]');
    clear legend_str
    xmin = 0;
    xmax = 0.5;
    %axis([xmin xmax xmin xmax]);
    plot([xmin xmax],[xmin xmax],'-k');
    axis auto
    x = true_RMSEn_pairs(isvalid);
    x = x(:);
    y = all_RMSE_pairs(isvalid);
    y = y(:);
    p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsqSVD = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-m');
    legend_str{1} = 'Data';
    legend_str{2} = '1:1';
    legend_str{3} = sprintf('Fit (R^2 = %3.2f)',rsqSVD);
    plot(xmax, 0.1, '.w');
    legend_str{4} = sprintf('  y = %3.2fx + %3.2f',p(1), p(2));
    legend(legend_str);
        set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
        export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
end
%}
%% Different RMSE metrics Linear 
if false;%options_showfig;
    fig_str = 'graph-anal-differentmetrics_linear';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    figure;
    %title('Different RMSE metrics?')
    %scatter(true_RMSE_MST(:),all_RMSE_MST(:),'+k')
    t_cut = 20;
    is_valid_sum  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_sum_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_sum_MST<t_cut&all_RMSE_sum_MST>10^-10;
    is_valid_geo  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_geo_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_geo_MST<t_cut&all_RMSE_geo_MST>10^-10;
    is_valid_mult  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_mult_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_mult_MST<t_cut&all_RMSE_mult_MST>10^-10;
    is_valid_max  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_max_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_max_MST<t_cut&all_RMSE_max_MST>10^-10;
    is_valid_path  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_path_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_path_MST<t_cut&all_RMSE_path_MST>10^-10;
    xmin = 0;
    xmax = 1;
    clear legend_str;
    plot([10^-4 xmax],[10^-4 xmax],'-k');
    legend_ctr = 1;
    legend_str{legend_ctr} = '1:1';
    legend_ctr = legend_ctr + 1;
    hold on
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
    axis([xmin xmax xmin xmax]);
    xstep = 10^-4:.05:xmax+.05;
    % plot fits : Linear
    % sum
    x = true_RMSE_MST(is_valid_sum);
    y = all_RMSE_sum_MST(is_valid_sum);
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(x,y,'.r');
    legend_str{legend_ctr} = 'Sum';
    legend_ctr = legend_ctr + 1;
    plot(xstep, yfit, '-r');
    legend_str{legend_ctr} =  sprintf('Sum fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % geo
    x = true_RMSE_MST(is_valid_geo);
    y = all_RMSE_geo_MST(is_valid_geo);
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(x,y,'.c');
    legend_str{legend_ctr} = 'Geo';
    legend_ctr = legend_ctr + 1;
    plot(xstep, yfit, '-c');
    legend_str{legend_ctr} =  sprintf('Geo fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % mult
    %{
    x = true_RMSE_MST(is_valid_mult);
    y = all_RMSE_mult_MST(is_valid_mult);
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(x,y,'.g');
    legend_str{legend_ctr} = 'Mult';
    legend_ctr = legend_ctr + 1;
    plot(xstep, yfit, '-g');
    legend_str{legend_ctr} =  sprintf('Mult fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %}
    % max
    x = true_RMSE_MST(is_valid_max);
    y = all_RMSE_max_MST(is_valid_max);
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(x,y,'.b');
    legend_str{legend_ctr} = 'Max';
    legend_ctr = legend_ctr + 1;
    plot(xstep, yfit, '-b');
    legend_str{legend_ctr} =  sprintf('Max fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % path
    x = true_RMSE_MST(is_valid_path);
    y = all_RMSE_path_MST(is_valid_path);
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(x,y,'.m');
    legend_str{legend_ctr} = 'Path';
    legend_ctr = legend_ctr + 1;
    plot(xstep, yfit, '-m');
    legend_str{legend_ctr} =  sprintf('Path fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %}
    % Plot fits: Half 
    %{
    % halfs
    plot(true_RMSE_MST(is_validh_sum),all_RMSEh_sum_MST(is_validh_sum),'+r')
    plot(true_RMSE_MST(is_validh_geo),all_RMSEh_geo_MST(is_validh_geo),'+c')
    plot(true_RMSE_MST(is_validh_mult),all_RMSEh_mult_MST(is_validh_mult),'+g')
    plot(true_RMSE_MST(is_validh_max),all_RMSEh_max_MST(is_validh_max),'+b')
    legend_str{12} = 'Sum h';
    legend_str{13} = 'Geo h';
    legend_str{14} = 'Mult h';
    legend_str{15} = 'Max h';    
    % sum
    %
    x = true_RMSE_MST(is_validh_sum);
    y = all_RMSEh_sum_MST(is_validh_sum);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{16} = sprintf('Sum h fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r');
    % geo
    x = true_RMSE_MST(is_validh_geo);
    y = all_RMSEh_geo_MST(is_validh_geo);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{17} = sprintf('Geo fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-c');
    % mult
    x = true_RMSE_MST(is_validh_mult);
    y = all_RMSEh_mult_MST(is_validh_mult);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{18} = sprintf('Mult fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % max
    x = true_RMSE_MST(is_validh_mult);
    y = all_RMSEh_max_MST(is_validh_mult);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend_str{19} = sprintf('Max fit (R^2 = %3.2f)',rsq);
    %} 
    legend(legend_str);
    set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
            export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
end
%% Sanity check 
isvalid = all_RMSE_path_MST>.01&(all_step<30);
 x = true_RMSE_MST(isvalid);
 y = all_RMSE_path_MST(isvalid);
 c = all_step(isvalid);
  figure;
  hold on
 scatter(x,y,10,c,'filled');
 [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
 xstep = [0 .5];
    yfit = polyval(p,xstep);
    plot(xstep, yfit, '-k');
    clear legend_str
    legend_str =  {'Data', sprintf('Max fit (R^2 = %3.2f)',stats(1))};
    legend(legend_str);
    
%% Different RMSE metrics LogLog 
if false;%options_showfig;
    fig_str = 'graph-anal-differentmetrics';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    figure;
    %title('Different RMSE metrics?')
    %scatter(true_RMSE_MST(:),all_RMSE_MST(:),'+k')
    t_cut = 20;
    is_valid_sum  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_sum_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_sum_MST<t_cut&all_RMSE_sum_MST>10^-10;
    is_valid_geo  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_geo_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_geo_MST<t_cut&all_RMSE_geo_MST>10^-10;
    is_valid_mult  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_mult_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_mult_MST<t_cut&all_RMSE_mult_MST>10^-10;
    is_valid_max  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_max_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_max_MST<t_cut&all_RMSE_max_MST>10^-10;
    is_valid_path  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_path_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_path_MST<t_cut&all_RMSE_path_MST>10^-10;
    % half 
    all_RMSEh_path_MST = all_RMSE_path_MST./2;
    is_validh_sum  = ~isnan(true_RMSE_MST)&~isinf(all_RMSEh_sum_MST)&...
        true_RMSE_MST<t_cut&all_RMSEh_sum_MST<t_cut&all_RMSEh_sum_MST>10^-10;
    is_validh_geo  = ~isnan(true_RMSE_MST)&~isinf(all_RMSEh_geo_MST)&...
        true_RMSE_MST<t_cut&all_RMSEh_geo_MST<t_cut&all_RMSEh_geo_MST>10^-10;
    is_validh_mult  = ~isnan(true_RMSE_MST)&~isinf(all_RMSEh_mult_MST)&...
        true_RMSE_MST<t_cut&all_RMSEh_mult_MST<t_cut&all_RMSEh_mult_MST>10^-10;
    is_validh_max  = ~isnan(true_RMSE_MST)&~isinf(all_RMSEh_max_MST)&...
        true_RMSE_MST<t_cut&all_RMSEh_max_MST<t_cut&all_RMSEh_max_MST>10^-10;
    is_validh_path  = ~isnan(true_RMSE_MST)&~isinf(all_RMSEh_path_MST)&...
        true_RMSE_MST<t_cut&all_RMSEh_path_MST<t_cut&all_RMSEh_path_MST>10^-10;
    xmin = 0;
    xmax = 1;
    clear legend_str;
    loglog([10^-4 xmax],[10^-4 xmax],'-k');
    legend_ctr = 1;
    legend_str{legend_ctr} = '1:1';
    legend_ctr = legend_ctr + 1;
    hold on
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
    axis([xmin xmax xmin xmax]);
    xstep = 10^-4:.05:xmax+.05;
    % plot fits : Linear
    %{
    plot(true_RMSE_MST(is_valid_sum),all_RMSE_sum_MST(is_valid_sum),'.r')
    plot(true_RMSE_MST(is_valid_geo),all_RMSE_geo_MST(is_valid_geo),'.c')
    plot(true_RMSE_MST(is_valid_mult),all_RMSE_mult_MST(is_valid_mult),'.g')
    plot(true_RMSE_MST(is_valid_max),all_RMSE_max_MST(is_valid_max),'.b')
    plot(true_RMSE_MST(is_valid_path),all_RMSE_path_MST(is_valid_path),'.m')
    % sum
    x = true_RMSE_MST(is_valid_sum);
    y = all_RMSE_sum_MST(is_valid_sum);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{7} = sprintf('Sum fit (R^2 = %3.2f)',rsq);
    plot([xmin:.1:xmax],  polyval(p,[xmin:.1:xmax]), '-k');
    % geo
    x = true_RMSE_MST(is_valid_geo);
    y = all_RMSE_geo_MST(is_valid_geo);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{8} = sprintf('Geo fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-c');
    % mult
    x = true_RMSE_MST(is_valid_mult);
    y = all_RMSE_mult_MST(is_valid_mult);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{9} = sprintf('Mult fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % max
    x = true_RMSE_MST(is_valid_mult);
    y = all_RMSE_max_MST(is_valid_mult);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend_str{10} = sprintf('Max fit (R^2 = %3.2f)',rsq);
    % path
    x = true_RMSE_MST(is_valid_path);
    y = all_RMSE_path_MST(is_valid_path);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-m');
    legend_str{11} = sprintf('Path fit (R^2 = %3.2f)',rsq);
    %}
    % Plot fits: Half 
    %{
    % halfs
    plot(true_RMSE_MST(is_validh_sum),all_RMSEh_sum_MST(is_validh_sum),'+r')
    plot(true_RMSE_MST(is_validh_geo),all_RMSEh_geo_MST(is_validh_geo),'+c')
    plot(true_RMSE_MST(is_validh_mult),all_RMSEh_mult_MST(is_validh_mult),'+g')
    plot(true_RMSE_MST(is_validh_max),all_RMSEh_max_MST(is_validh_max),'+b')
    legend_str{12} = 'Sum h';
    legend_str{13} = 'Geo h';
    legend_str{14} = 'Mult h';
    legend_str{15} = 'Max h';    
    % sum
    %
    x = true_RMSE_MST(is_validh_sum);
    y = all_RMSEh_sum_MST(is_validh_sum);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{16} = sprintf('Sum h fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r');
    % geo
    x = true_RMSE_MST(is_validh_geo);
    y = all_RMSEh_geo_MST(is_validh_geo);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{17} = sprintf('Geo fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-c');
    % mult
    x = true_RMSE_MST(is_validh_mult);
    y = all_RMSEh_mult_MST(is_validh_mult);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{18} = sprintf('Mult fit (R^2 = %3.2f)',rsq);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % max
    x = true_RMSE_MST(is_validh_mult);
    y = all_RMSEh_max_MST(is_validh_mult);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend_str{19} = sprintf('Max fit (R^2 = %3.2f)',rsq);
    %} 
    % Plot fits: Log Log to overcome Heteroscedasticity 
    %{
    % sum
    x = true_RMSE_MST(is_valid_sum);
    y = all_RMSE_sum_MST(is_valid_sum);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.r');
    legend_str{legend_ctr} = 'Sum';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-r');
    legend_str{legend_ctr} =  sprintf('Sum fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % geo
    x = true_RMSE_MST(is_valid_geo);
    y = all_RMSE_geo_MST(is_valid_geo);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.c');
    loglog(xstep, yfit, '-c');
    legend_str{legend_ctr} =  sprintf('Geo fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % mult
    %{
    x = true_RMSE_MST(is_valid_mult);
    y = all_RMSE_mult_MST(is_valid_mult);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.g');
    legend_str{legend_ctr} = 'Mult';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-g');
    legend_str{legend_ctr} =  sprintf('Mult fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %}
    % max
    x = true_RMSE_MST(is_valid_max);
    y = all_RMSE_max_MST(is_valid_max);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.b');
    legend_str{legend_ctr} = 'Max';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-b');
    legend_str{legend_ctr} =  sprintf('Max fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % path
    x = true_RMSE_MST(is_valid_path);
    y = all_RMSE_path_MST(is_valid_path);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.m');
    legend_str{legend_ctr} = 'Path';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-m');
    legend_str{legend_ctr} =  sprintf('Path fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %}
    % Plot half fits: Log Log
    % sum
    x = true_RMSE_MST(is_validh_sum);
    y = all_RMSEh_sum_MST(is_validh_sum);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.r');
    legend_str{legend_ctr} = 'Sum half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-r');
    legend_str{legend_ctr} =  sprintf('Sum half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % geo
    x = true_RMSE_MST(is_validh_geo);
    y = all_RMSEh_geo_MST(is_validh_geo);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.c');
    legend_str{legend_ctr} = 'Geo half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-c');
    legend_str{legend_ctr} =  sprintf('Geo half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % mult
    %{
    x = true_RMSE_MST(is_validh_mult);
    y = all_RMSEh_mult_MST(is_validh_mult);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'+g');
    legend_str{legend_ctr} = 'Mult half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, ':g');
    legend_str{legend_ctr} =  sprintf('Mult half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %}
    % max
    x = true_RMSE_MST(is_validh_max);
    y = all_RMSEh_max_MST(is_validh_max);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.b');
    legend_str{legend_ctr} = 'Max half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-b');
    legend_str{legend_ctr} =  sprintf('Max half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % path
    x = true_RMSE_MST(is_validh_path);
    y = all_RMSEh_path_MST(is_validh_path);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.m');
    legend_str{legend_ctr} = 'Path half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-m');
    legend_str{legend_ctr} =  sprintf('Path half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % 
    legend(legend_str);
    set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
            export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
end
%% Only Path RMSE out vs true RMSE
%{
if options_showfig;
    filepath_RMSEtruevsRMSEoutpath = sprintf('%sRMSEdifferentmetrics.png',path_fig);
    figure;
    hold on
    title('Different RMSE metrics?')
    scatter(true_RMSE_MST(:),all_RMSE_path_MST(:),'.m')
    xmin = 0;
    xmax = 1;
    plot([xmin xmax],[xmin xmax],'-k');
    legend_str = {'Path', '1:1','Fit'};
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
    axis([xmin xmax xmin xmax]);
    % plot fits
    % path
    isvalidx = ~isnan(true_RMSE_MST);
    isvalidy = ~isinf(all_RMSE_path_MST);
    isvalid = isvalidx & isvalidy;
    x = true_RMSE_MST(isvalid);
    y = all_RMSE_path_MST(isvalid);
    p = polyfit(x,y,1);
    %p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-m');
    legend_str{11} = sprintf('Path fit (R^2 = %3.2f)',rsq);
    legend(legend_str);
    %export_fig(filepath_RMSEtruevsRMSEout);
    end
%}
%% RMSE out vs RMSE true
if false;%~options_showfig;
    fig_str = 'graph-anal-RMSEoutvsRMSEtrue';
    filepath_png = sprintf('%s%s.png',path_fig, fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz, fig_str);
    figure;
    %title('Is reported RMSE accurate predictor of true RMSE?')
    t_cut = 10;
    is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(all_RMSE_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_MST<t_cut&all_RMSE_MST>10^-3;
    is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(all_RMSE_WMF)&...
        true_RMSE_WMF<t_cut&all_RMSE_WMF<t_cut&all_RMSE_WMF>10^-3;    
    is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(all_RMSE_SVD)&...
        true_RMSE_SVD<t_cut&all_RMSE_SVD<t_cut&all_RMSE_SVD>10^-3;
    xmin = 0;
    xmax = 5;
    loglog([10^-4 xmax],[10^-4 xmax],'-k');
    hold on
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
    axis([0 1 0 1]);
    clear legend_str
    legend_ctr = 1;
    % Linear plot fits
    %{
    plot(true_RMSE_MST(is_valid_MST),all_RMSE_MST(is_valid_MST),'.r')
    plot(true_RMSE_WMF(is_valid_WMF),all_RMSE_WMF(is_valid_WMF),'.g')
    plot(true_RMSE_SVD(is_valid_SVD),all_RMSE_SVD(is_valid_SVD),'.b')
    % MST
    x = true_RMSE_MST(is_valid_MST);
    y = all_RMSE_MST(is_valid_MST);
    % p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsqMST = 1 - SSresid/SStotal;
    legend_str{5} = sprintf('MST fit (R^2 = %3.2f)',rsqMST);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r');
    % WMF
    x = true_RMSE_WMF(is_valid_WMF);
    y = all_RMSE_WMF(is_valid_WMF);
    %p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsqWMF = 1 - SSresid/SStotal;
    legend_str{6} = sprintf('WMF fit (R^2 = %3.2f)',rsqWMF);
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % SVD
    x = true_RMSE_SVD(is_valid_SVD);
    y = all_RMSE_SVD(is_valid_SVD);
    %p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsqSVD = 1 - SSresid/SStotal;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend_str{7} = sprintf('SVD fit (R^2 = %3.2f)',rsqSVD);
    legend(legend_str);
    %}
    % LogLog plot fits 
    x = true_RMSE_MST(is_valid_MST);
    y = all_RMSE_MST(is_valid_MST);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.r');
    legend_str{legend_ctr} = 'MST';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-r');
    legend_str{legend_ctr} =  sprintf('MST fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % WMF
    x = true_RMSE_WMF(is_valid_WMF);
    y = all_RMSE_WMF(is_valid_WMF);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.g');
    legend_str{legend_ctr} = 'WMF';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-g');
    legend_str{legend_ctr} =  sprintf('WMF fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % SVD 
    x = true_RMSE_SVD(is_valid_SVD);
    y = all_RMSE_SVD(is_valid_SVD);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.b');
    legend_str{legend_ctr} = 'SVD';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-b');
    legend_str{legend_ctr} =  sprintf('SVD fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
        set(gcf, 'color','w')
        legend(legend_str);
    if ~exist(filepath_png, 'file');
            export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
end
%% RMSE out vs RMSE in
%{
if options_showfig;
filepath_RMSEoutvsRMSEin = sprintf('%sRMSEoutvsRMSEin.png',path_fig);
figure;
hold on
scatter(all_P_RMSEin(:), all_RMSE_MST(:), '.r');
scatter(all_P_RMSEin(:), all_RMSE_WMF(:), '.g');
scatter(all_P_RMSEin(:), all_RMSE_SVD(:), '.b');
xlabel('Input RMSE');
ylabel('Output RMSE');
legend_str{1} = 'MST';
legend_str{2} = 'WMF';
legend_str{3} = 'SVD';
legend(legend_str);
set(gcf, 'color' ,'w');
%export_fig(filepath_RMSEoutvsRMSEin);
end
%}
%% RMSE true vs RMSE in
if options_showfig;
    fig_str = 'graph-anal-RMSEtruevsRMSEin';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    figure;
    hold on
    legend_ctr = 1;
    clear legend_str;
    %ctitle('Which method is best?')
    t_cut = 1;
    is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(all_P_RMSEin)&...
        true_RMSE_MST<t_cut&all_P_RMSEin<t_cut;
    is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(all_P_RMSEin)&...
        true_RMSE_WMF<t_cut&all_P_RMSEin<t_cut;
    is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(all_P_RMSEin)&...
        true_RMSE_SVD<t_cut&all_P_RMSEin<t_cut;
    plot(all_P_RMSEin(is_valid_MST),true_RMSE_MST(is_valid_MST),'.r')
    legend_str{legend_ctr} = 'MST';
    legend_ctr = legend_ctr + 1;
    plot(all_P_RMSEin(is_valid_WMF),true_RMSE_WMF(is_valid_WMF),'.g')
        legend_str{legend_ctr} = 'WMF';
    legend_ctr = legend_ctr + 1;
    plot(all_P_RMSEin(is_valid_SVD),true_RMSE_SVD(is_valid_SVD),'.b')
        legend_str{legend_ctr} = 'SVD';
    legend_ctr = legend_ctr + 1;
    xmin = 0;
    xmax = 0.1;
    ymin = 0;
    ymax = 1;
    plot([ymin ymax],[ymin ymax],'-k');
    legend_str{legend_ctr} = '1:1';
    legend_ctr = legend_ctr + 1; 
    xlabel('Input RMSE');
    ylabel('True RMSE');
    axis([xmin xmax ymin ymax]);
    set(gcf, 'color' ,'w');
    % MST
    x = all_P_RMSEin(is_valid_MST);
    y = true_RMSE_MST(is_valid_MST);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{legend_ctr} = 'MST fit';
    legend_ctr = legend_ctr + 1;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r');
    % WMF
    x = all_P_RMSEin(is_valid_WMF);
    y = true_RMSE_WMF(is_valid_WMF);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{legend_ctr} = 'WMF fit';
    legend_ctr = legend_ctr + 1;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % SVD
    x = all_P_RMSEin(is_valid_SVD);
    y = true_RMSE_SVD(is_valid_SVD);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{legend_ctr} = 'SVD fit';
    legend_ctr = legend_ctr + 1;
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend(legend_str);
        set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
            export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
    %export_fig(filepath_RMSEtruevsRMSEin);
end
%% RMSE graph vs RMSE pair  (true)
if options_showfig;
    fig_str = 'graph-anal-RMSEgraphvsRMSEpair';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    figure;
    hold on
    %title('Are graphical estimates improved?')
    t_cut = 1;
    is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(true_RMSE_pair)&...
        true_RMSE_MST<t_cut&true_RMSE_pair<t_cut;
    is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(true_RMSE_pair)&...
        true_RMSE_WMF<t_cut&true_RMSE_pair<t_cut;
    is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(true_RMSE_pair)&...
        true_RMSE_SVD<t_cut&true_RMSE_pair<t_cut;
    plot(true_RMSE_pair(is_valid_MST),true_RMSE_MST(is_valid_MST),'.r')
    plot(true_RMSE_pair(is_valid_WMF),true_RMSE_WMF(is_valid_WMF),'.g')
    plot(true_RMSE_pair(is_valid_SVD),true_RMSE_SVD(is_valid_SVD),'.b')
    xmin = 0;
    xmax = 0.5;
    plot([0 xmax],[0 xmax],'-k');
    axis([xmin xmax xmin xmax]);
    legend_str = {'MST','WMF', 'SVD','1:1'};
    xlabel('True Pairwise RMSE [m]');
    ylabel('True Graph RMSE [m]');
    legend(legend_str);
    isnnanpair = ~isnan(true_RMSE_pair);
    isvalid = (true_RMSE_pair<1);
    %
    % MST
    x = true_RMSE_pair(is_valid_MST);
    y = true_RMSE_MST(is_valid_MST);
    %p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{5} = sprintf('MST fit');
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r');
    % WMF
    x = true_RMSE_pair(is_valid_WMF);
    y = true_RMSE_WMF(is_valid_WMF);
    % p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{6} = sprintf('WMF fit');
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % SVD
    x = true_RMSE_pair(is_valid_SVD);
    y = true_RMSE_SVD(is_valid_SVD);
    %p = polyfit(x,y,1);
    p = polyfitZero(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{7} = sprintf('SVD fit');
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend(legend_str);
    axis([0 .5 0 .5]);
    set(gcf, 'color' ,'w');
        set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
            export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
    %}
    %export_fig(filepath_RMSEtruevsRMSEin);
end
%% Percent connected graph vs pair
%{
foo = 1;
true_num_pair = zeros(n_step,n_trial);
true_num_MST = zeros(n_step,n_trial);
true_num_WMF = zeros(n_step,n_trial);
true_num_SVD = zeros(n_step,n_trial);
for step = 1:n_step;
    for trial = 1:n_trial;
        temp_pair = squeeze(true_RMSE_pair(step,trial,:));
        temp_MST = squeeze(true_RMSE_MST(step,trial,:));
        temp_WMF = squeeze(true_RMSE_WMF(step,trial,:));
        temp_SVD = squeeze(true_RMSE_SVD(step,trial,:));
        true_num_pair(step,trial) =  sum((~isnan(temp_pair) & temp_pair<1));
        true_num_MST(step,trial) =  sum((~isnan(temp_MST) & temp_MST<1));
        true_num_WMF(step,trial) =  sum((~isnan(temp_WMF) & temp_WMF<1));
        true_num_SVD(step,trial) =  sum((~isnan(temp_SVD) & temp_SVD<1));
    end
end
true_perc_pair = 100*true_num_pair./n_S;
true_perc_MST = 100*true_num_MST./n_S;
true_perc_WMF = 100*true_num_WMF./n_S;
true_perc_SVD = 100*true_num_SVD./n_S;
figure;
hold on
plot(true_perc_pair(:), true_perc_MST(:), '.r')
plot(true_perc_pair(:), true_perc_WMF(:), '.g')
plot(true_perc_pair(:), true_perc_SVD(:), '.b')
plot([0 100],[0 100],'-k');
axis([0 100 0 100]);
xlabel('Percent detected Pairwise');
ylabel('Percent detected Graphical');
%}
%% RMSE graph vs RMSE pair  (reported)
%{
if options_showfig;
    filepath_RMSEgraphvsRMSEpair = sprintf('%sRMSEgraphvsRMSEpair.png',path_fig);
    figure;
    hold on
    title('Are graphical estimates improved?')
    scatter(all_RMSE_pair(:),all_RMSE_MST(:),'.r')
    scatter(all_RMSE_pair(:),all_RMSE_WMF(:),'.g')
    scatter(all_RMSE_pair(:),all_RMSE_SVD(:),'.b')
    plot([0 .5],[0 .5],'-k');
    legend_str = {'MST','WMF', 'SVD','1:1'};
    xlabel('Reported Pairwise RMSE');
    ylabel('Reported Graph RMSE');
    legend(legend_str);
    axis([0 0.5 0 0.5]);
    set(gcf, 'color' ,'w');
    %export_fig(filepath_RMSEtruevsRMSEin);
end
%}
%% Error vs range
if options_showfig;
    fig_str = 'graph-anal-RMSEvsrange';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    figure;
    hold on
        t_cut = 60;
    is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(all_S_range)&...
        true_RMSE_MST<t_cut&all_S_range<t_cut&all_S_range>2;
    is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(all_S_range)&...
        true_RMSE_WMF<t_cut&all_S_range<t_cut&all_S_range>2;
    is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(all_S_range)&...
        true_RMSE_SVD<t_cut&all_S_range<t_cut&all_S_range>2;
    plot(all_S_range(isvalid),true_RMSE_MST(isvalid),'.r')
    plot(all_S_range(isvalid),true_RMSE_WMF(isvalid),'.g')
    plot(all_S_range(isvalid),true_RMSE_SVD(isvalid),'.b')
    clear legend_str
    legend_str{1} = 'MST';
    legend_str{2} = 'WMF';
    legend_str{3} = 'SVD';
    ylabel('True RMSE [m]');
    xlabel('Sensor Range [m]');
    legend(legend_str);
    xmin = 0;
    xmax = 60;
    axis([xmin xmax 0 1]);
    set(gcf, 'color' ,'w');
    % MST
    x = all_S_range(is_valid_MST);
    y = true_RMSE_MST(is_valid_MST);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{4} = sprintf('MST fit');
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r');
    % WMF
    x = all_S_range(is_valid_WMF);
    y = true_RMSE_WMF(is_valid_WMF);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{5} = sprintf('WMF fit');
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g');
    % SVD
    x = all_S_range(is_valid_SVD);
    y = true_RMSE_SVD(is_valid_SVD);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    legend_str{6} = sprintf('SVD fit');
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b');
    legend(legend_str);
        set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
            export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
end
%%





