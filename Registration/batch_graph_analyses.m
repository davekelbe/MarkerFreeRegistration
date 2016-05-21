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
aux.info_exp = 'GraphAnalRMSE13';
sigma_mean = .1;
n_trial = 30;
aux.info_exp = 'GraphAnalnewtruth';
sigma_mean = 0:.01:.1;
n_trial = 1;
aux.info_exp = 'GraphAnalnewtruthfix';
sigma_mean = 0:.03:.1;
n_trial = 1;
aux.info_exp = 'GraphAnalRMSE14';
sigma_mean = 0:.005:.1;
n_trial = 1;
aux.info_exp = 'GraphAnalRMSE15';
sigma_mean = logspace(-5,-.6);
sigma_mean = sigma_mean(1:5:end);
n_trial = 1;
%}
%{
aux.info_exp = 'GraphAnalRMSE16';
sigma_mean = logspace(-5,-.3,100);
n_trial = 1;
%sigma_mean  = sigma_mean(22:end);
n_step = numel(sigma_mean);
%}
aux.info_exp = 'GraphAnalRMSE17';
sigma_mean = logspace(-5,-.3,100);
n_trial = 1;
n_step = numel(sigma_mean);
%}
%{
aux.info_exp = 'GraphAnalRMSE18';
sigma_mean = logspace(-5,-.3,100);
n_trial = 1;
n_step = numel(sigma_mean);
%}

% SENSITIVITY ANALYSIS STUDY 
aux.sensitivity_r = true;

for step = 85:85;%n_step;
    fprintf('\n Step %g of %g\n',step, n_step);
    aux.info_suffix = sprintf('step-%02.0f',step);
    aux.S_n = sigma_mean(step);
    for trial = 1:n_trial;
        aux.trial = sprintf('%02.0f',trial);
        graph_analysis_rmsefun_generic( aux )
    end
end

path_fig = sprintf('%s%s%s','Z:\Desktop\KelbeGraph\',aux.info_exp,'\');
if ~exist(path_fig, 'dir');
    mkdir(path_fig);
end
path_tikz = sprintf('%s%s%s','Z:\Documents\Research\Documents\KelbeThesis\thesis\');
%% Notx analysis

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
        %match12_RMSE(logical(eye(n_S))) = nan;
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
        all_R_SVD{step,trial}{w} = [];
        all_t_SVD{step,trial}{w} = [];
        filepath_G_t_SVD = sprintf('%s%s',path_mat, 'G_t_SVD.mat');
        load(filepath_G_t_SVD);
        all_t_SVD{step,trial} = G_t_SVD;
        filepath_G_RMSE_SVD = sprintf('%s%s',path_mat, 'G_RMSE_SVD.mat');
        load(filepath_G_RMSE_SVD);
        all_RMSE_SVD(step,trial,:)= G_RMSE_SVD;
    end
end
all_RMSE_MST = all_RMSE_sum_MST;
all_n_T = cellfun(@numel,all_T_x);

if aux.sensitivity_r;
    trial_temp = 1;
    step_temp = 85;
    aux.info_suffix_temp = sprintf('step-%02.0f',step_temp);
    aux.trial_temp = sprintf('%02.0f',trial_temp);
    path_mat_temp = sprintf('%s%s%s%s%s%s%s','D:\Users\djk2312\Documents\Harvard\reg\',aux.info_exp,'\', aux.info_suffix_temp,'\',aux.trial_temp,'\');
    filepath_SENS_ANAL_G_rx_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_rx_WMF.mat');
    load(filepath_SENS_ANAL_G_rx_WMF);
    filepath_SENS_ANAL_G_ry_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_ry_WMF.mat');
    load(filepath_SENS_ANAL_G_ry_WMF);
    filepath_SENS_ANAL_G_rz_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_rz_WMF.mat');
    load(filepath_SENS_ANAL_G_rz_WMF);
    filepath_SENS_ANAL_G_tx_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_tx_WMF.mat');
    load(filepath_SENS_ANAL_G_tx_WMF);
    filepath_SENS_ANAL_G_ty_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_ty_WMF.mat');
    load(filepath_SENS_ANAL_G_ty_WMF);
    filepath_SENS_ANAL_G_tz_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_tz_WMF.mat');
    load(filepath_SENS_ANAL_G_tz_WMF);
    filepath_SENS_ANAL_G_RMSE_WMF = sprintf('%s%s',path_mat_temp, 'SENS_ANAL_G_RMSE_WMF.mat');
    load(filepath_SENS_ANAL_G_RMSE_WMF);
    SENS_RES_G_tx_WMF = nanmean(std(SENS_ANAL_G_tx_WMF,0,2));
    SENS_RES_G_ty_WMF = nanmean(std(SENS_ANAL_G_ty_WMF,0,2));
    SENS_RES_G_tz_WMF = nanmean(std(SENS_ANAL_G_tz_WMF,0,2));
    SENS_RES_G_rx_WMF = nanmean(std(SENS_ANAL_G_rx_WMF,0,2));
    SENS_RES_G_ry_WMF = nanmean(std(SENS_ANAL_G_ry_WMF,0,2));
    SENS_RES_G_rz_WMF = nanmean(std(SENS_ANAL_G_rz_WMF,0,2));
end
%% PP array
all_PP_LCSn = cell(n_step, n_trial, n_S);
all_PP_LCS = cell(n_step, n_trial, n_S);
all_PP_WCSn = cell(n_step, n_trial, n_S);

for step = 1:n_step;
    for trial = 1:n_trial;
        for i = 1:n_S;
            all_PP_WCSn{step,trial,i} = nan(3,all_n_T(step,trial));
            all_PP_WCSn{step,trial,i}(:,all_P_ix{step,trial}{i}) = all_P_WCSn{step,trial}{i};
            all_PP_LCSn{step,trial,i} = nan(3,all_n_T(step,trial));
            all_PP_LCSn{step,trial,i}(:,all_P_ix{step,trial}{i}) = all_P_LCSn{step,trial}{i};
            all_PP_LCS{step,trial,i} = nan(3,all_n_T(step,trial));
            all_PP_LCS{step,trial,i}(:,all_P_ix{step,trial}{i}) = all_P_LCS{step,trial}{i};
        end
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

all_PP_LCSnWi_pairs = cell(n_step, n_trial, n_S,n_S);
all_P_LCSnw_pair = cell(n_step, n_trial, n_S);
all_P_LCSni_pairs = cell(n_step, n_trial, n_S,n_S);
all_P_LCSi_pairs = cell(n_step, n_trial, n_S,n_S);
all_PP_LCSni_pairs = cell(n_step, n_trial, n_S,n_S);
all_PP_LCSi_pairs = cell(n_step, n_trial, n_S,n_S);
all_P_LCSnj12 = cell(n_step, n_trial, n_S);
all_P_LCSnW12 = cell(n_step, n_trial, n_S);
all_P_LCSnw_MST = cell(n_step, n_trial, n_S);
all_PP_LCSw_MST = cell(n_step, n_trial, n_S);
all_PP_LCSw_WMF = cell(n_step, n_trial, n_S);
all_PP_LCSw_SVD = cell(n_step, n_trial, n_S);
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
                all_PP_LCSw_MST{step,trial,j} = all_R_MST{step, trial}{j}*...
                    all_PP_LCS{step,trial,j} + ...
                    repmat(...
                    all_t_MST{step, trial}{j},[1 all_n_T(step,trial)]);
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
                all_PP_LCSw_WMF{step,trial,j} = all_R_WMF{step, trial}{j}*...
                    all_PP_LCS{step,trial,j} + ...
                    repmat(...
                    all_t_WMF{step, trial}{j},[1 all_n_T(step,trial)]);
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
                all_PP_LCSw_SVD{step,trial,j} = all_R_SVD{step, trial}{j}*...
                    all_PP_LCS{step,trial,j} + ...
                    repmat(...
                    all_t_SVD{step, trial}{j},[1 all_n_T(step,trial)]);
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
                    %        all_PP_LCSnWi_pairs{step,trial,i,j} = all_S_R{step,trial}{i}'*...
                    %            all_PP_WCSn{step,trial,j} - all_S_R{step,trial}{i}'*repmat(...
                    %            all_S_t{step,trial}{i}, [1, size(all_PP_WCSn{step,trial,j},2)]);
                    all_P_LCSni_pairs{step,trial,i,j} = all_R_pairs{step, trial}{i,j}*...
                        all_P_LCSn{step,trial}{j} + ...
                        repmat(...
                        all_t_pairs{step, trial}{i,j},[1 P_ntree(j)]);
                    all_PP_LCSni_pairs{step,trial,i,j} = all_R_pairs{step, trial}{i,j}*...
                        all_PP_LCSn{step,trial,j} + ...
                        repmat(...
                        all_t_pairs{step, trial}{i,j},[1 all_n_T(step,trial)]);
                    all_PP_LCSi_pairs{step,trial,i,j} = all_R_pairs{step, trial}{i,j}*...
                        all_PP_LCS{step,trial,j} + ...
                        repmat(...
                        all_t_pairs{step, trial}{i,j},[1 all_n_T(step,trial)]);
                    all_P_LCSi_pairs{step,trial,i,j} = all_R_pairs{step, trial}{i,j}*...
                        all_P_LCS{step,trial}{j} + ...
                        repmat(...
                        all_t_pairs{step, trial}{i,j},[1 P_ntree(j)]);
                end
            end
        end
    end
end
%% Sanity check
%{
step = 3;
trial = 1;
i = 1;
j = 2;
figure;
scatter3(all_PP_LCSnWi_pairs{step,trial,i,j}(1,:), ...
    all_PP_LCSnWi_pairs{step,trial,i,j}(2,:),...
    all_PP_LCSnWi_pairs{step,trial,i,j}(3,:));
hold on
j = 3;
plot3(all_PP_LCSnWi_pairs{step,trial,i,j}(1,:), ...
    all_PP_LCSnWi_pairs{step,trial,i,j}(2,:),...
    all_PP_LCSnWi_pairs{step,trial,i,j}(3,:), '+');
j = 4;
plot3(all_PP_LCSnWi_pairs{step,trial,i,j}(1,:), ...
    all_PP_LCSnWi_pairs{step,trial,i,j}(2,:),...
    all_PP_LCSnWi_pairs{step,trial,i,j}(3,:), 'x');
%}
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
true_RMSE_MST_NEW = nan(n_step, n_trial, n_S);
true_RMSE_WMF = nan(n_step, n_trial, n_S);
true_RMSE_SVD = nan(n_step, n_trial, n_S);
true_RMSE_pair = nan(n_step, n_trial, n_S);
true_RMSE_pair2 = nan(n_step, n_trial, n_S);
true_RMSEnnS_pairs = nan(n_step, n_trial, n_S,n_S);
true_RMSEnn_pairs = nan(n_step, n_trial, n_S,n_S);
true_RMSEn1_pairs = nan(n_step, n_trial, n_S,n_S);
true_RMSEn2_pairs = nan(n_step, n_trial, n_S,n_S);
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
                %{
                true_RMSE_MST(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_MST{step,trial,j}).^2)));
                %}
                true_RMSE_MST(step,trial,j) = ...
                    sqrt(nanmean(sum((all_PP_LCS{step,trial,w} -...
                    all_PP_LCSw_MST{step,trial,j}).^2)));
                
            end
            % WMF
            if isvalid_WMF(step,trial,j);
                %{
                true_RMSE_WMF(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_WMF{step,trial,j}).^2)));
                %}
                true_RMSE_WMF(step,trial,j) = ...
                    sqrt(nanmean(sum((all_PP_LCS{step,trial,w} -...
                    all_PP_LCSw_WMF{step,trial,j}).^2)));
            end
            % SVD
            if isvalid_SVD(step,trial,j);
                %{
                true_RMSE_SVD(step,trial,j) = ...
                    sqrt(nanmean(sum((all_P_WCSn{step,trial}{j} -...
                    all_P_LCSnW_SVD{step,trial,j}).^2)));
                %}
                true_RMSE_SVD(step,trial,j) = ...
                    sqrt(nanmean(sum((all_PP_LCS{step,trial,w} -...
                    all_PP_LCSw_SVD{step,trial,j}).^2)));
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
                    %      true_RMSEnnS_pairs(step,trial,i,j) = ...
                    %          sqrt(nanmean(sum((all_PP_LCSn{step,trial,i} -...
                    %          all_PP_LCSnWi_pairs{step,trial,i,j}).^2)));
                    % Exploratory noise noise
                    true_RMSEnn_pairs(step,trial,i,j) = ...
                        sqrt(nanmean(sum((all_PP_LCSn{step,trial,i} -...
                        all_PP_LCSni_pairs{step,trial,i,j}).^2)));
                    % Exploratory noise no noise
                    true_RMSEn1_pairs(step,trial,i,j) = ...
                        sqrt(nanmean(sum((all_PP_LCS{step,trial,i} -...
                        all_PP_LCSni_pairs{step,trial,i,j}).^2)));
                    % Exploratory no noise noise
                    true_RMSEn2_pairs(step,trial,i,j) = ...
                        sqrt(nanmean(sum((all_PP_LCSn{step,trial,i} -...
                        all_PP_LCSi_pairs{step,trial,i,j}).^2)));
                    % Used in pub fig
                    true_RMSE_pairs(step,trial,i,j) = ...
                        sqrt(nanmean(sum((all_PP_LCS{step,trial,i} -...
                        all_PP_LCSi_pairs{step,trial,i,j}).^2)));
                end
            end
        end
    end
end
%% Sanity check
%{
step = 10;
trial = 1;
j = 2;
i = 1;
figure;
hold on
axis equal
view(0,90);
for t = 1:all_n_T(step);
    if isnan(all_PP_LCS{step,trial,i}(1,t)); continue; end
    h = filledCircle([all_PP_LCS{step,trial,i}(1,t) all_PP_LCS{step,trial,i}(2,t)],...
        all_T_r{step,trial}(t),1000, 'r');
    set(h, 'FaceAlpha', 0.2);
end
for t = 1:all_n_T(step);
    if isnan(all_PP_LCSi_pairs{step,trial,i,j}(1,t)); continue; end
    h = filledCircle([all_PP_LCSi_pairs{step,trial,i,j}(1,t) all_PP_LCSi_pairs{step,trial,i,j}(2,t)],...
        all_T_r{step,trial}(t),1000, 'b');
    set(h, 'FaceAlpha', 0.2);
end
%}

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
%% SENSITIVITY STUDY ON exponential parameter 



%%
%{
figure;
isvalidnnS = ~isnan(true_RMSEnnS_pairs) & ~isnan(all_RMSE_pairs);
plot(true_RMSEnnS_pairs(isvalidnnS), all_RMSE_pairs(isvalidnnS), '.k')
%}
%% Which truth data is best? Linear
%{
figure;
hold on
clear legend_str
legend_ctr = 1;
isvalidnn = ~isnan(true_RMSEnn_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1;
xnn = true_RMSEnn_pairs(isvalidnn);
ynn = all_RMSE_pairs(isvalidnn);
scatter(xnn,ynn,5, 'r', 'filled');
legend_str{legend_ctr} = 'NN';
legend_ctr = legend_ctr + 1;
isvalidn1 = ~isnan(true_RMSEn1_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1;
xn1 = true_RMSEn1_pairs(isvalidn1);
yn1 = all_RMSE_pairs(isvalidn1);
scatter(xn1,yn1,5, 'g', 'filled');
legend_str{legend_ctr} = 'N-';
legend_ctr = legend_ctr + 1;
isvalid = ~isnan(true_RMSE_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1;
x = true_RMSE_pairs(isvalid);
y = all_RMSE_pairs(isvalid);
scatter(x,y,5, 'b', 'filled');
legend_str{legend_ctr} = '--';
legend_ctr = legend_ctr + 1;
plot([0 .4],[ 0 .4], '-k');
legend_str{legend_ctr} = '1:1';
legend_ctr = legend_ctr + 1;
xstep = [0 .4];
[p,~,~,~,stats] = regress(ynn,[xnn ones(size(xnn))]);
yfit = polyval(p,xstep);
plot(xstep, yfit, '-r');
legend_str{legend_ctr} =  sprintf('NN Fit (R^2 = %3.2f)',stats(1));
legend_ctr = legend_ctr + 1;
[p,~,~,~,stats] = regress(yn1,[xn1 ones(size(xn1))]);
yfit = polyval(p,xstep);
plot(xstep, yfit, '-g');
legend_str{legend_ctr} =  sprintf('N1 Fit (R^2 = %3.2f)',stats(1));
legend_ctr = legend_ctr + 1;
[p,~,~,~,stats] = regress(y,[x ones(size(x))]);
yfit = polyval(p,xstep);
plot(xstep, yfit, '-b');
legend_str{legend_ctr} =  sprintf('-- Fit (R^2 = %3.2f)',stats(1));
legend_ctr = legend_ctr + 1;
legend(legend_str);
xlabel('True Pairwise RMSE [m]');
ylabel('Reported Pairwise RMSE [m]');
%}
%% Which truth data is best Log Log ?
%{
figure;
xstep = logspace(-5,-.6); %10^-4:.05:xmax+.05;
clear legend_str
legend_ctr = 1;
isvalidnn = ~isnan(true_RMSEnn_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1&all_RMSE_pairs>10^-5;
xnn = true_RMSEnn_pairs(isvalidnn);
ynn = all_RMSE_pairs(isvalidnn);
loglog(xnn,ynn,'.r');
hold on
legend_str{legend_ctr} = 'NN';
legend_ctr = legend_ctr + 1;
isvalidn1 = ~isnan(true_RMSEn1_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1&all_RMSE_pairs>10^-5;
xn1 = true_RMSEn1_pairs(isvalidn1);
yn1 = all_RMSE_pairs(isvalidn1);
loglog(xn1,yn1,'.g');
legend_str{legend_ctr} = 'N-';
legend_ctr = legend_ctr + 1;
isvalid = ~isnan(true_RMSE_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1&all_RMSE_pairs>10^-5;
x = true_RMSE_pairs(isvalid);
y = all_RMSE_pairs(isvalid);
loglog(x,y,'.b');
legend_str{legend_ctr} = '--';
legend_ctr = legend_ctr + 1;
plot(xstep,xstep, '-k');
legend_str{legend_ctr} = '1:1';
legend_ctr = legend_ctr + 1;
%
[p,~,~,~,stats] = regress(log(ynn),[log(xnn) ones(size(xnn))]);
yfit = exp(polyval(p,log(xstep)));
loglog(xstep, yfit, '-r');
legend_str{legend_ctr} =  sprintf('NN Fit (R^2 = %3.2f)',stats(1));
legend_ctr = legend_ctr + 1;
%
[p,~,~,~,stats] = regress(log(yn1),[log(xn1) ones(size(xn1))]);
yfit = exp(polyval(p,log(xstep)));
loglog(xstep, yfit, '-g');
legend_str{legend_ctr} =  sprintf('N- Fit (R^2 = %3.2f)',stats(1));
legend_ctr = legend_ctr + 1;
%
[p,~,~,~,stats] = regress(log(ynn),[log(x) ones(size(x))]);
yfit = exp(polyval(p,log(xstep)));
loglog(xstep, yfit, '-b');
legend_str{legend_ctr} =  sprintf('-- Fit (R^2 = %3.2f)',stats(1));
legend_ctr = legend_ctr + 1;
%
xlim([xstep(1) xstep(end)])
loglog(xstep, xstep, '-k');
legend_str{legend_ctr} =  sprintf('1:1');
legend_ctr = legend_ctr + 1;
legend(legend_str, 'location', 'northwest');
xlabel('True Pairwise RMSE [m]');
ylabel('Reported Pairwie RMSE [m]');
%}
%% Reported Pairwise vs True Pairwise - Log Log - Pub
ms = 8;
lw = 2;
if options_showfig
    figure;
    fig_str = 'graph-anal-pairoutvspairtrue-log';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    xstep = logspace(-5,-.6);
    clear legend_str
    legend_ctr = 1;
    isvalid = ~isnan(true_RMSE_pairs) & ~isnan(all_RMSE_pairs) & all_RMSE_pairs<1&all_RMSE_pairs>10^-5;
    x = true_RMSE_pairs(isvalid);
    y = all_RMSE_pairs(isvalid);%./3.1930;
    issub = randperm(numel(x),min(3000, numel(x)));
    loglog(x(issub),y(issub),'.','Color', [.5 .5 .5], 'markersize', ms);
    hold on
    legend_str{legend_ctr} = 'Data';
    legend_ctr = legend_ctr + 1;
    % Regression
    %[p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    %mdl =fitlm(log(x), log(y));
    % t_crit
    btest = [1 exp(1)];
    [p, stats, t0slope, t0int, tcrit ] = regress_with_test( log(x),log(y),btest );
    if abs(t0int) > tcrit;
        fprintf('\n Reject null hypothesis intercept\n');
    else
        fprintf('\n Fail to reject null hypothesis intercept');
    end
    if abs(t0slope) > tcrit;
        fprintf('\n Reject null hypothesis slope\n');
    else
        fprintf('\n Fail to reject null hypothesis slope');
    end
    %yresid = sqrt((log(y) - polyval(p,log(x))).^2);
    %figure; scatter(log(x),log(y), 10, yresid, 'filled'); colorbar
    %[plin,~,~,~,stats] = regress(y,[x ones(size(x))]);
    %yresidlin = sqrt((y - polyval(plin,x)).^2);
    %yfitlin = polyval(plin,xstep);
    %figure; scatter(x,y,10,yresidlin, 'filled'); hold on
    %plot(xstep,yfitlin, '-r');
    yfit = exp(polyval(p,log(xstep)));
    loglog(xstep, yfit, '-k', 'linewidth', lw);
    legend_str{legend_ctr} =  sprintf('Fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % 1:1
    loglog(xstep,xstep, '--', 'Color', [.3 .3 .3], 'linewidth', lw);
    legend_str{legend_ctr} = '1:1';
    legend_ctr = legend_ctr + 1;
    xlim([xstep(1) xstep(end)])
    legend(legend_str, 'location', 'northwest');
    xlabel('True Pairwise RMSE [m]');
    ylabel('Reported Pairwie RMSE [m]');
    set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
        export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        matlab2tikz(filepath_tikz);
    end
end
%% Different RMSE metrics Linear
%{
if options_showfig;
    clear legend_str
    legend_ctr = 1;
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
    %is_valid_mult  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_mult_MST)&...
    %    true_RMSE_MST<t_cut&all_RMSE_mult_MST<t_cut&all_RMSE_mult_MST>10^-10;
    is_valid_max  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_max_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_max_MST<t_cut&all_RMSE_max_MST>10^-10;
    is_valid_path  = ~isnan(true_RMSE_MST)&~isinf(all_RMSE_path_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_path_MST<t_cut&all_RMSE_path_MST>10^-10;
    xstep = logspace(-5,-.6);
    xmin = xstep(1);
    xmax = xstep(end);
    hold on
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
    axis([xmin xmax xmin xmax]);
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
    % Plot fits: Half
    %
    % sum
    x = true_RMSE_MST(is_validh_sum);
    y = all_RMSEh_sum_MST(is_validh_sum);
    plot(x,y,'+r')
    legend_str{legend_ctr} =  sprintf('Sum half');
    legend_ctr = legend_ctr + 1;
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(xstep, yfit, '--r');
    legend_str{legend_ctr} =  sprintf('Sum half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % geo
     x = true_RMSE_MST(is_validh_geo);
    y = all_RMSEh_geo_MST(is_validh_geo);
    plot(x,y,'+c')
    legend_str{legend_ctr} =  sprintf('Geo half');
    legend_ctr = legend_ctr + 1;
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(xstep, yfit, '--c');
    legend_str{legend_ctr} =  sprintf('Geo half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % max
    x = true_RMSE_MST(is_validh_max);
    y = all_RMSEh_max_MST(is_validh_max);
    plot(x,y,'+b')
    legend_str{legend_ctr} =  sprintf('Max half');
    legend_ctr = legend_ctr + 1;
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(xstep, yfit, '--b');
    legend_str{legend_ctr} =  sprintf('Max half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % Path
    x = true_RMSE_MST(is_validh_path);
    y = all_RMSEh_path_MST(is_validh_path);
    plot(x,y,'+m')
    legend_str{legend_ctr} =  sprintf('Path half');
    legend_ctr = legend_ctr + 1;
    [p,~,~,~,stats] = regress(y,[x ones(size(x))]);
    yfit = polyval(p,xstep);
    plot(xstep, yfit, '--m');
    legend_str{legend_ctr} =  sprintf('Path half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %
    plot([xmin xmax], [xmin xmax], '-k');
    legend_str{legend_ctr} =  sprintf('1:1');
    legend_ctr = legend_ctr + 1;
    legend(legend_str, 'location', 'southeast');
    set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
        %  export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        %   matlab2tikz(filepath_tikz);
    end
end
%}
%% Sanity check
%{
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
%}
%% Different RMSE metrics LogLog
%{
if options_showfig;
    clear legend_str
    fig_str = 'graph-anal-differentmetrics';
    filepath_png = sprintf('%s%s.png',path_fig,fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
    figure;
    xstep = logspace(-5,-.6);
    xmin = xstep(1);
    xmax = xstep(end);
    %title('Different RMSE metrics?')
    %scatter(true_RMSE_MST(:),all_RMSE_MST(:),'+k')
    t_cut = 20;
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
    clear legend_str;
    loglog([xmin xmax],[xmin xmax],'-k', 'linewidth', lw);
    legend_ctr = 1;
    legend_str{legend_ctr} = '1:1';
    legend_ctr = legend_ctr + 1;
    hold on
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
    axis([xmin xmax xmin xmax]);
    % Plot half fits
    % sum
    x = true_RMSE_MST(is_validh_sum);
    y = all_RMSEh_sum_MST(is_validh_sum);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.r', 'markersize', ms);
    legend_str{legend_ctr} = 'Sum half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-r', 'linewidth', lw);
    legend_str{legend_ctr} =  sprintf('Sum half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % geo
    x = true_RMSE_MST(is_validh_geo);
    y = all_RMSEh_geo_MST(is_validh_geo);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.c', 'markersize', ms);
    legend_str{legend_ctr} = 'Geo half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-c', 'linewidth', lw);
    legend_str{legend_ctr} =  sprintf('Geo half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % max
    x = true_RMSE_MST(is_validh_max);
    y = all_RMSEh_max_MST(is_validh_max);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.b', 'markersize', ms);
    legend_str{legend_ctr} = 'Max half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-b', 'linewidth', lw);
    legend_str{legend_ctr} =  sprintf('Max half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    % path
    x = true_RMSE_MST(is_validh_path);
    y = all_RMSEh_path_MST(is_validh_path);
    [p,~,~,~,stats] = regress(log(y),[log(x) ones(size(x))]);
    yfit = exp(polyval(p,log(xstep)));
    loglog(x,y,'.m', 'markersize', ms);
    legend_str{legend_ctr} = 'Path half';
    legend_ctr = legend_ctr + 1;
    hold on
    loglog(xstep, yfit, '-m', 'linewidth', lw);
    legend_str{legend_ctr} =  sprintf('Path half fit (R^2 = %3.2f)',stats(1));
    legend_ctr = legend_ctr + 1;
    %
    legend(legend_str, 'location', 'northwest');
    set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
        %export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        % matlab2tikz(filepath_tikz);
    end
end
%}
%% RMSE out vs RMSE true Pub
if options_showfig;
    clear legend_str
    fig_str = 'graph-anal-RMSEoutvsRMSEtrue';
    filepath_png = sprintf('%s%s.png',path_fig, fig_str);
    filepath_tikz = sprintf('%s%s.tex',path_tikz, fig_str);
    figure;
    %title('Is reported RMSE accurate predictor of true RMSE?')
    t_cut = inf;
    is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(all_RMSE_MST)&...
        true_RMSE_MST<t_cut&all_RMSE_MST<t_cut&all_RMSE_MST>10^-5;
    is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(all_RMSE_WMF)&...
        true_RMSE_WMF<t_cut&all_RMSE_WMF<t_cut&all_RMSE_WMF>10^-5;
    is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(all_RMSE_SVD)&...
        true_RMSE_SVD<t_cut&all_RMSE_SVD<t_cut&all_RMSE_SVD>10^-5;
    xstep = logspace(-5,-.6);
    xmin = xstep(1);
    xmax = xstep(end);
    xlim([xmin xmax])
    clear legend_str
    legend_ctr = 1;
    % LogLog plot fits
    x_MST = true_RMSE_MST(is_valid_MST);
    y_MST = all_RMSE_MST(is_valid_MST);
    btest = [1 0]; 
    [p, stats_MST, t0slope, t0int, tcrit ] = regress_with_test( log(x_MST),log(y_MST),btest );
    if abs(t0int) > tcrit;
        fprintf('\n MST Reject null hypothesis intercept\n');
    else
        fprintf('\n MST Fail to reject null hypothesis intercept');
    end
    if abs(t0slope) > tcrit;
        fprintf('\n MST Reject null hypothesis slope\n');
    else
        fprintf('\n MST Fail to reject null hypothesis slope');
    end
    yfit_MST = exp(polyval(p,log(xstep)));
    loglog(x_MST,y_MST,'o', 'markersize', ms/2, 'color', [.75 .75 .75], 'markerfacecolor', [.85 .85 .85]);
    hold on
    legend_str{legend_ctr} = 'MST';
    legend_ctr = legend_ctr + 1;
    % WMF
    x_WMF = true_RMSE_WMF(is_valid_WMF);
    y_WMF = all_RMSE_WMF(is_valid_WMF);
    btest = [1 0]; 
    [p, stats_WMF, t0slope, t0int, tcrit ] = regress_with_test( log(x_WMF),log(y_WMF),btest );
    if abs(t0int) > tcrit;
        fprintf('\n WMF Reject null hypothesis intercept\n');
    else
        fprintf('\n WMF Fail to reject null hypothesis intercept');
    end
    if abs(t0slope) > tcrit;
        fprintf('\n WMF Reject null hypothesis slope\n');
    else
        fprintf('\n WMF Fail to reject null hypothesis slope');
    end
    yfit_WMF = exp(polyval(p,log(xstep)));
    loglog(x_WMF,y_WMF,'s', 'markersize', ms/2, 'color', [.6 .6 .6], 'markerfacecolor', [.7 .7 .7]);
    legend_str{legend_ctr} = 'WMF';
    legend_ctr = legend_ctr + 1;
    % SVD
    x_SVD = true_RMSE_SVD(is_valid_SVD);
    y_SVD = all_RMSE_SVD(is_valid_SVD);
    btest = [1 0]; 
    [p, stats_SVD, t0slope, t0int, tcrit ] = regress_with_test( log(x_SVD),log(y_SVD),btest );
    if abs(t0int) > tcrit;
        fprintf('\n SVD Reject null hypothesis intercept\n');
    else
        fprintf('\n SVD Fail to reject null hypothesis intercept');
    end
    if abs(t0slope) > tcrit;
        fprintf('\n SVD Reject null hypothesis slope\n');
    else
        fprintf('\n SVD Fail to reject null hypothesis slope');
    end
    yfit_SVD = exp(polyval(p,log(xstep)));
    loglog(x_SVD,y_SVD,'^', 'markersize', ms/2, 'color', [.35 .35 .35], 'markerfacecolor', [.45 .45 .45]);
    legend_str{legend_ctr} = 'SVD';
    legend_ctr = legend_ctr + 1;
    % MST
    loglog(xstep, yfit_MST, '-o', 'MarkerEdgeColor','k',...
        'MarkerFaceColor',[.85 .85 .85], 'linewidth', lw/1.5, 'color', 'k');
    legend_str{legend_ctr} =  sprintf('MST fit (R^2 = %3.2f)',stats_MST(1));
    legend_ctr = legend_ctr + 1;
    % WMF
    loglog(xstep, yfit_WMF, '-s', 'MarkerEdgeColor','k',...
        'MarkerFaceColor',[.7 .7 .7], 'linewidth', lw/1.5, 'color', 'k');
    legend_str{legend_ctr} =  sprintf('WMF fit (R^2 = %3.2f)',stats_WMF(1));
    legend_ctr = legend_ctr + 1;
    % SVD
    loglog(xstep, yfit_SVD, '-^',  'MarkerEdgeColor','k',...
        'MarkerFaceColor',[.45 .45 .45], 'linewidth', lw/1.5, 'color', 'k');
    legend_str{legend_ctr} =  sprintf('SVD fit (R^2 = %3.2f)',stats_SVD(1));
    legend_ctr = legend_ctr + 1;
    % 1:1
    loglog([xmin xmax],[xmin xmax],'--', 'color', [.3 .3 .3], 'linewidth', lw);
    legend_str{legend_ctr} = '1:1';
    legend_ctr = legend_ctr + 1;
    % Appearance
    set(gcf, 'color','w')
    xlim([xmin xmax]);
    legend(legend_str, 'location', 'northwest');
    xlabel('True RMSE [m]');
    ylabel('Reported RMSE [m]');
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
%% RMSE true vs RMSE in linear
%{
if options_showfig;
    clear legend_str
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
    plot(all_P_RMSEin(is_valid_MST),true_RMSE_MST(is_valid_MST),'.r', 'markersize',ms)
    legend_str{legend_ctr} = 'MST';
    legend_ctr = legend_ctr + 1;
    plot(all_P_RMSEin(is_valid_WMF),true_RMSE_WMF(is_valid_WMF),'.g', 'markersize',ms)
    legend_str{legend_ctr} = 'WMF';
    legend_ctr = legend_ctr + 1;
    plot(all_P_RMSEin(is_valid_SVD),true_RMSE_SVD(is_valid_SVD),'.b', 'markersize',ms)
    legend_str{legend_ctr} = 'SVD';
    legend_ctr = legend_ctr + 1;
    xmin = 0;
    xmax = 0.1;
    ymin = 0;
    ymax = 1;
    plot([ymin ymax],[ymin ymax],'--','color', [.3 .3 .3], 'linewidth', lw);
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
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-r', 'linewidth', lw);
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
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-g', 'linewidth', lw);
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
    plot([xmin xmax],  polyval(p,[xmin xmax]), '-b', 'linewidth', lw);
    legend(legend_str, 'location', 'northwest');
    set(gcf, 'color','w')
    if ~exist(filepath_png, 'file');
        %export_fig(filepath_png);
    end
    if ~exist(filepath_tikz, 'file');
        %matlab2tikz(filepath_tikz);
    end
    %export_fig(filepath_RMSEtruevsRMSEin);
end
%}
%% RMSE true vs RMSE in Pub log log
clear legend_str
fig_str = 'graph-anal-RMSEtruevsRMSEinloglog';
filepath_png = sprintf('%s%s.png',path_fig, fig_str);
filepath_tikz = sprintf('%s%s.tex',path_tikz, fig_str);
figure;
%title('Is reported RMSE accurate predictor of true RMSE?')
t_cut = inf;
is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(all_P_RMSEin)&...
    true_RMSE_MST<t_cut&all_P_RMSEin<t_cut&all_RMSE_MST>10^-5;
is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(all_P_RMSEin)&...
    true_RMSE_WMF<t_cut&all_P_RMSEin<t_cut&all_RMSE_MST>10^-5;
is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(all_P_RMSEin)&...
    true_RMSE_SVD<t_cut&all_P_RMSEin<t_cut&all_RMSE_MST>10^-5;
xstep = logspace(-5,-.6);
xmin = xstep(1);
xmax = xstep(end);
xlim([xmin xmax])
clear legend_str
legend_ctr = 1;
% LogLog ploft fits
x_MST = all_P_RMSEin(is_valid_MST);
y_MST = true_RMSE_MST(is_valid_MST);
btest = [1 0]; 
[p, stats_MST, t0slope, t0int, tcrit ] = regress_with_test( log(x_MST),log(y_MST),btest );
if abs(t0int) > tcrit;
    fprintf('\n MST Reject null hypothesis intercept\n');
else
    fprintf('\n MST Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n MST Reject null hypothesis slope\n');
else
    fprintf('\n MST Fail to reject null hypothesis slope\n');
end
yfit_MST = exp(polyval(p,log(xstep)));
loglog(x_MST,y_MST,'o', 'markersize', ms/2, 'color', [.75 .75 .75], 'markerfacecolor', [.85 .85 .85]);
hold on
legend_str{legend_ctr} = 'MST';
legend_ctr = legend_ctr + 1;
% WMF
x_WMF = all_P_RMSEin(is_valid_WMF);
y_WMF = true_RMSE_WMF(is_valid_WMF);
btest = [1 0]; 
[p, stats_WMF, t0slope, t0int, tcrit ] = regress_with_test( log(x_WMF),log(y_WMF),btest );
if abs(t0int) > tcrit;
    fprintf('\n WMF Reject null hypothesis intercept\n');
else
    fprintf('\n WMF Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n WMF Reject null hypothesis slope\n');
else
    fprintf('\n WMF Fail to reject null hypothesis slope\n');
end
yfit_WMF = exp(polyval(p,log(xstep)));
loglog(x_WMF,y_WMF,'s', 'markersize', ms/2, 'color', [.6 .6 .6], 'markerfacecolor', [.7 .7 .7]);
legend_str{legend_ctr} = 'WMF';
legend_ctr = legend_ctr + 1;
% SVD
x_SVD = all_P_RMSEin(is_valid_SVD);
y_SVD = true_RMSE_SVD(is_valid_SVD);
btest = [1 0]; 
[p, stats_SVD, t0slope, t0int, tcrit ] = regress_with_test( log(x_SVD),log(y_SVD),btest );
if abs(t0int) > tcrit;
    fprintf('\n SVD Reject null hypothesis intercept\n');
else
    fprintf('\n SVD Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n SVD Reject null hypothesis slope\n');
else
    fprintf('\n SVD Fail to reject null hypothesis slope\n');
end
yfit_SVD = exp(polyval(p,log(xstep)));
loglog(x_SVD,y_SVD,'^', 'markersize', ms/2, 'color', [.35 .35 .35], 'markerfacecolor', [.45 .45 .45]);
legend_str{legend_ctr} = 'SVD';
legend_ctr = legend_ctr + 1;
% MST
loglog(xstep, yfit_MST, '-o', 'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.85 .85 .85], 'linewidth', lw/1.5, 'color', 'k');
legend_str{legend_ctr} =  sprintf('MST fit (R^2 = %3.2f)',stats_MST(1));
legend_ctr = legend_ctr + 1;
% WMF
loglog(xstep, yfit_WMF, '-s', 'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.7 .7 .7], 'linewidth', lw/1.5, 'color', 'k');
legend_str{legend_ctr} =  sprintf('WMF fit (R^2 = %3.2f)',stats_WMF(1));
legend_ctr = legend_ctr + 1;
% SVD
loglog(xstep, yfit_SVD, '-^',  'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.45 .45 .45], 'linewidth', lw/1.5, 'color', 'k');
legend_str{legend_ctr} =  sprintf('SVD fit (R^2 = %3.2f)',stats_SVD(1));
legend_ctr = legend_ctr + 1;
% 1:1
loglog([xmin xmax],[xmin xmax],'--', 'color', [.3 .3 .3], 'linewidth', lw);
legend_str{legend_ctr} = '1:1';
legend_ctr = legend_ctr + 1;
% Appearance
set(gcf, 'color','w')
xlim([xmin xmax]);
legend(legend_str, 'location', 'northwest');
xlabel('RMSE in [m]');
ylabel('True RMSE [m]');
if ~exist(filepath_png, 'file');
    export_fig(filepath_png);
end
if ~exist(filepath_tikz, 'file');
    matlab2tikz(filepath_tikz);
end
%% RMSE graph vs RMSE pair  (true)
%{
clear legend_str
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
%% RMSE true vs RMSE pair true in Pub log log
clear legend_str
fig_str = 'graph-anal-RMSEtruevsRMSEpairloglog';
filepath_png = sprintf('%s%s.png',path_fig, fig_str);
filepath_tikz = sprintf('%s%s.tex',path_tikz, fig_str);
figure;
%title('Is reported RMSE accurate predictor of true RMSE?')
t_cut = 1;
is_valid_MST  = ~isnan(true_RMSE_MST)&~isnan(true_RMSE_pair)&...
    true_RMSE_MST<t_cut&true_RMSE_pair<t_cut&true_RMSE_MST>10^-5;
is_valid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(true_RMSE_pair)&...
    true_RMSE_WMF<t_cut&true_RMSE_pair<t_cut&true_RMSE_WMF>10^-5;
is_valid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(true_RMSE_pair)&...
    true_RMSE_SVD<t_cut&true_RMSE_pair<t_cut&true_RMSE_SVD>10^-5;
xstep = logspace(-5,-.6);
xmin = xstep(1);
xmax = xstep(end);
xlim([xmin xmax])
clear legend_str
legend_ctr = 1;
% LogLog plot fits
x_MST = true_RMSE_pair(is_valid_MST);
y_MST = true_RMSE_MST(is_valid_MST);
btest = [1 0]; 
[p, stats_MST, t0slope, t0int, tcrit ] = regress_with_test( log(x_MST),log(y_MST),btest );
if abs(t0int) > tcrit;
    fprintf('\n MST Reject null hypothesis intercept\n');
else
    fprintf('\n MST Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n MST Reject null hypothesis slope\n');
else
    fprintf('\n MST Fail to reject null hypothesis slope\n');
end
yfit_MST = exp(polyval(p,log(xstep)));
loglog(x_MST,y_MST,'o', 'markersize', ms/2, 'color', [.75 .75 .75], 'markerfacecolor', [.85 .85 .85]);
hold on
legend_str{legend_ctr} = 'MST';
legend_ctr = legend_ctr + 1;
% WMF
x_WMF = true_RMSE_pair(is_valid_WMF);
y_WMF = true_RMSE_WMF(is_valid_WMF);
btest = [1 0]; 
[p, stats_WMF, t0slope, t0int, tcrit ] = regress_with_test( log(x_WMF),log(y_WMF),btest );
if abs(t0int) > tcrit;
    fprintf('\n WMF Reject null hypothesis intercept\n');
else
    fprintf('\n WMF Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n WMF Reject null hypothesis slope\n');
else
    fprintf('\n WMF Fail to reject null hypothesis slope\n');
end
yfit_WMF = exp(polyval(p,log(xstep)));
loglog(x_WMF,y_WMF,'s', 'markersize', ms/2, 'color', [.6 .6 .6], 'markerfacecolor', [.7 .7 .7]);
legend_str{legend_ctr} = 'WMF';
legend_ctr = legend_ctr + 1;
% SVD
x_SVD = true_RMSE_pair(is_valid_SVD);
y_SVD = true_RMSE_SVD(is_valid_SVD);
btest = [1 0]; 
[p, stats_SVD, t0slope, t0int, tcrit ] = regress_with_test( log(x_SVD),log(y_SVD),btest );
if abs(t0int) > tcrit;
    fprintf('\n SVD Reject null hypothesis intercept\n');
else
    fprintf('\n SVD Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n SVD Reject null hypothesis slope\n');
else
    fprintf('\n SVD Fail to reject null hypothesis slope\n');
end
yfit_SVD = exp(polyval(p,log(xstep)));
loglog(x_SVD,y_SVD,'^', 'markersize', ms/2, 'color', [.35 .35 .35], 'markerfacecolor', [.45 .45 .45]);
legend_str{legend_ctr} = 'SVD';
legend_ctr = legend_ctr + 1;
% MST
loglog(xstep, yfit_MST, '-o', 'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.85 .85 .85], 'linewidth', lw/1.5, 'color', 'k');
legend_str{legend_ctr} =  sprintf('MST fit (R^2 = %3.2f)',stats_MST(1));
legend_ctr = legend_ctr + 1;
% WMF
loglog(xstep, yfit_WMF, '-s', 'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.7 .7 .7], 'linewidth', lw/1.5, 'color', 'k');
legend_str{legend_ctr} =  sprintf('WMF fit (R^2 = %3.2f)',stats_WMF(1));
legend_ctr = legend_ctr + 1;
% SVD
loglog(xstep, yfit_SVD, '-^',  'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.45 .45 .45], 'linewidth', lw/1.5, 'color', 'k');
legend_str{legend_ctr} =  sprintf('SVD fit (R^2 = %3.2f)',stats_SVD(1));
legend_ctr = legend_ctr + 1;
% 1:1
loglog([xmin xmax],[xmin xmax],'--', 'color', [.1 .1 .1], 'linewidth', lw);
legend_str{legend_ctr} = '1:1';
legend_ctr = legend_ctr + 1;
% Appearance
set(gcf, 'color','w')
xlim([xmin xmax]);
legend(legend_str, 'location', 'northwest');
xlabel('True Pairwise RMSE [m]');
ylabel('True Graph RMSE [m]');
if ~exist(filepath_png, 'file');
    export_fig(filepath_png);
end
if ~exist(filepath_tikz, 'file');
    matlab2tikz(filepath_tikz);
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
fig_str = 'graph-anal-RMSEvsrange';
filepath_png = sprintf('%s%s.png',path_fig,fig_str);
filepath_tikz = sprintf('%s%s.tex',path_tikz,fig_str);
figure;
hold on
t_cut = 60;
isvalid_MST  = ~isnan(true_RMSE_MST)&~isnan(all_S_range)&...
    true_RMSE_MST<t_cut&all_S_range<t_cut&all_S_range>-inf;
isvalid_WMF  = ~isnan(true_RMSE_WMF)&~isnan(all_S_range)&...
    true_RMSE_WMF<t_cut&all_S_range<t_cut&all_S_range>-inf;
isvalid_SVD  = ~isnan(true_RMSE_SVD)&~isnan(all_S_range)&...
    true_RMSE_SVD<t_cut&all_S_range<t_cut&all_S_range>-inf;
plot(all_S_range(isvalid_MST),true_RMSE_MST(isvalid_MST),'.r')
plot(all_S_range(isvalid_WMF),true_RMSE_WMF(isvalid_WMF),'.g')
plot(all_S_range(isvalid_SVD),true_RMSE_SVD(isvalid_SVD),'.b')
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
%%





