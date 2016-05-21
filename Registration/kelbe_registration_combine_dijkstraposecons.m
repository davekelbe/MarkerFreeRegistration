function [  ] = kelbe_registration_combine_dijkstraposecons( aux )
% Only load valid plots
%info_valid_plot = {'13','15'};
%      info_valid_plot = cell(1,25);
%      for s = 1:5;
%          info_valid_plot{s} = sprintf('%02.0f', s);
%      end
%      info_site = 31;
info_valid_plot = aux.info_valid_plot;
info_site = aux.info_site;

tree = 1;
%% KELBE REGISTRATION
% Clear without breakpoints
% Default settings
set(0,'defaultfigureposition', [895   169   760   651]')
options_verbose = true;
options_imagepoints = false;
options_initialmatch = false;
options_unique = false;
options_loadmatch = false;
path_save = 'Z:\Desktop\Registration\Figures\';
path_tikz = 'Z:\Desktop\Registration\tikz\';
path_matvar = 'D:\Users\djk2312\Documents\matvar\';
options_loadvar =false;
%if ~options_loadvar;
% Initialize program-level options and paths
% Outputs
% options               options
% paths                 paths
clear D ctr d filepath_tree plot t tree
%% Parameters
t_rad = 0.2; %0.06^2; % Trees can be matched if their radius is within threshold
t_coll = 0.05; % Collinnearity threshold: don't consider triangles that are collinear
t_eig_error = 1e1;
t_RANSAC_xyz = 0.4^2;
t_RANSAC_nsearch = 100;
t_flagr = 20; % deg
t_flagt = 1; % m


% Initialize program-level options and paths
% Outputs
% t_*               thresholds
%% Load points
% Input to function are stem maps derived from lidar.m
fprintf('\nLoad points\n');

% Filename Lookup for Harvard Data
info_exp = 'Harvard';
info_suffix = 'reg';
info_slash = '\';
path_site = sprintf('%s%s%s%s%s%03.0f%s','D:\Users\djk2312\Documents\',...
    info_exp, info_slash, info_suffix,info_slash,info_site,info_slash);
D = dir(path_site);
ctr = 1;
% Only load valid plots
n_S = numel(info_valid_plot);
n_tree = zeros(n_S,1);
P_LCS = cell(n_S,1);
P_rad = cell(n_S,1);
P_plot = zeros(n_S,1);
path_ply = cell(n_S,1);
filepath_ply = cell(n_S,1);

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
    if numel(tree)==1;
        if isempty(tree.r)
            continue
        end
    end
    if numel(tree)<=6;
        continue
    end
    n_tree(ctr) = numel(tree);
    plot = str2double(info_plot);
    P_LCS{ctr} = nan(3,n_tree(ctr));
    P_rad{ctr} = nan(n_tree(ctr),1);
    P_plot(ctr) = plot;
    for t = 1:n_tree(ctr);
        P_LCS{ctr}(:,t) = tree(t).loc(:,1);
        P_rad{ctr}(t) = tree(t).r(1);
    end
    path_ply{ctr} = sprintf('%s%s%sply%s',path_site,info_plot,info_slash,info_slash);
    filepath_ply{ctr} = sprintf('%spoints_full_%03.0f-%02.0f.ply', path_ply{ctr}, info_site, plot);
    ctr = ctr + 1;
end
isvalid = (P_plot~=0);
P_LCS = P_LCS(isvalid);
P_rad = P_rad(isvalid);
P_plot = P_plot(isvalid);
n_S = sum(isvalid);
%[P_LCS, P_rad, P_plot, n_tree, n_S ] = generate_test_data;

i_xmin = min(cellfun(@(x) min(x(1,:)),P_LCS));
i_xmax = max(cellfun(@(x) max(x(1,:)),P_LCS));
i_ymin = min(cellfun(@(x) min(x(2,:)),P_LCS));
i_ymax = max(cellfun(@(x) max(x(2,:)),P_LCS));
i_zmin = min(cellfun(@(x) min(x(3,:)),P_LCS));
i_zmax = max(cellfun(@(x) max(x(3,:)),P_LCS));
options_axesval = [i_xmin i_xmax i_ymin i_ymax i_zmin i_zmax];

% Colormap for sensors
P_color = jet(n_S);

% Individual camera views
if false;%options_verbose && options_imagepoints;
    for s = 1:n_S;
        clear legend_str
        figure
        hold on
        plot3(0,0,0,'^k','markersize',10,...
            'markerfacecolor',P_color(s,:));
        hdummy = plot3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),'ok','markersize',5,...
            'markerfacecolor',P_color(s,:));
        set(hdummy, 'visible', 'off');
        for t = 1:numel(P_rad{s});
            h = filledCircle([P_LCS{s}(1,t); P_LCS{s}(2,t)]',P_rad{s}(t),1000,P_color(s,:));
        end
        %scatter3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),30,...
        %    color_P_index(truth_P_index{s},:),'filled');
        %axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
        axis auto
        axisval = axis;
        %set(gca, 'xtick',
        xlabel('x Position relative to plot center [m]');
        ylabel('y Position relative to plot center [m]');
        zlabel('z Position relative to plot center [m]');
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
        % filepath_save = sprintf('%sLCS_%02.0f.eps',path_save, P_plot(s));
        % saveas(gcf,filepath_save,'psc2')
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
clear D ctr d filepath_tree plot t tree i_xmax i_xmin i_ymin i_ymax i_zmin i_zmax
%if false;
%% Find radius similarity
%{
    match_israd = cell(n_S);
    match_rad_diff = cell(n_S);
    for i = 1:n_S-1;
        for j = i+1:n_S; % Symmetric
            rep_radI = repmat(P_rad{i},[1,n_tree(j)]);
            rep_radJ = repmat(P_rad{j}',[n_tree(i),1]);
            rad_diff = abs((rep_radI - rep_radJ)./((rep_radI + rep_radJ)/2));
            match_rad_diff{i,j} = rad_diff;
            match_israd{i,j} = (rad_diff < t_rad);
            match_israd{j,i} =  match_israd{i,j}';
        end
    end
    
    clear i j rad_diff rep_radI rep_radJ
    % Outputs
    % match_israd                cell array with {i,j} logical array of matching radii
%}
%% Find connection network
%{
    % All possible triangles combinations for the full sets
    
    T_comb = cell(n_S,1);
    
    if max(n_tree) < 255;
        FDT_T_comb = @(x) uint8(x);
    elseif max(n_tree) < 65535;
        FDT_T_comb = @(x) uint16(x);
    end
    
    for s = 1:n_S;
        T_comb{s} = FDT_T_comb(combnk(1:n_tree(s),3));
    end
    
    [n_comb,~] = cellfun(@size,T_comb);
    
%{
s = 1;
hold on
for c = 1:200;% n_comb(s);
    plot([P_LCS{s}(1,T_comb{s}(c,1)) P_LCS{s}(1,T_comb{s}(c,2))],...
         [P_LCS{s}(2,T_comb{s}(c,1)) P_LCS{s}(2,T_comb{s}(c,2))],'-k')
    plot([P_LCS{s}(1,T_comb{s}(c,1)) P_LCS{s}(1,T_comb{s}(c,3))],...
         [P_LCS{s}(2,T_comb{s}(c,1)) P_LCS{s}(2,T_comb{s}(c,3))],'-k')
    plot([P_LCS{s}(1,T_comb{s}(c,2)) P_LCS{s}(1,T_comb{s}(c,3))],...
         [P_LCS{s}(2,T_comb{s}(c,2)) P_LCS{s}(2,T_comb{s}(c,3))],'-k')
end
%}
    
    clear c s
    % Outputs:
    %       T_comb{i}               - Combinations of tree's; n_P_comb{i} x 3
    %       n_comb(i)               - number of combinations of trees
    %       FDT_T_comb              - FDT = Datatype function
%}
%% Sort triangle set combinations by decreasing radius
%{
    T_LCS = cell(n_S,1);
    T_R = cell(n_S,1);
    
    for s = 1:n_S;
        T_R_temp =  reshape(P_rad{s}(T_comb{s}),[n_comb(s),3]);
        [T_R{s}, sortix] = sort(T_R_temp,2,'descend');
        rowix = repmat((1:n_comb(s))',[1,3]);
        linear = sub2ind([n_comb(s),3],rowix(:), sortix(:));
        T_comb{s} = reshape(T_comb{s}(linear), n_comb(s),3); % Redefine T_comb
        T_LCS{s} = zeros(n_comb(s),3,3);
        T_LCS{s}(:,:,1) = reshape(P_LCS{s}(1,T_comb{s}),[n_comb(s),3]); % x values
        T_LCS{s}(:,:,2) = reshape(P_LCS{s}(2,T_comb{s}),[n_comb(s),3]); % y values
        T_LCS{s}(:,:,3) = reshape(P_LCS{s}(3,T_comb{s}),[n_comb(s),3]); % z values
    end
    
    clear s linear rowix T_R_temp sortix
    % Outputs:
    %       T_LCS(i)               - Triangle points sorted
    %       T_R(i)                 - Triangle radii sorted
    %       T_comb{s}              - [UPDATED] sorted combinations based on radius
%}
%% Find eigenvalues of each triangle
%{
    T_eig = cell(n_S,1);
    for s = 1:n_S;
        T_eig{s} = zeros(n_comb(s),2);
        for t = 1:n_comb(s);
            temp = eig(cov(squeeze(T_LCS{s}(t,:,:))));
            T_eig{s}(t,:) = temp(2:3);
        end
    end
    
    T_normeig = cell(n_S,1);
    T_isncoll = cell(n_S,1);
    for s = 1:n_S;
        T_normeig{s} = T_eig{s}./repmat(sum(T_eig{s},2),[1,2]);
        T_isncoll{s} = T_normeig{s}(:,1)> t_coll;
    end
    
%{
%% Histogram of normalized (smaller) eigenvalues
figure;
edges = linspace(0,1,20);
histogram(T_normeig{1}(:,1),edges)
xlabel('percent variance');
ylabel('count');
%% Simulate data with different covariance matrices
white_data = [cosd(30) -sind(30); -cosd(30) -sind(30); 0 1];
percent_var = 0.05:0.05:0.75;
n_var = numel(percent_var);
var_color = jet(n_var);
figure;
hold on
legend_str = cell(n_var,1);
for i = 1:n_var
    C = [percent_var(i) 0; 0 1-percent_var(i)];
    legend_str{i} = sprintf('%2.2f', percent_var(i));
    % C = [.25 0; 0 1.1607];
    % C = C./trace(C);
    L = chol(C);
    cov_data = white_data*L;
    scatter(cov_data(:,1), cov_data(:,2),20,var_color(i,:),'filled');
end
for i = 1:n_var;
     C = [percent_var(i) 0; 0 1-percent_var(i)];
    L = chol(C);
    cov_data = white_data*L;
	plot([cov_data(:,1); cov_data(1,1)],[cov_data(:,2); cov_data(1,2)],'color',var_color(i,:),'linewidth',2);
end
legend(legend_str);
axis equal
%}
    
    clear s t temp edges legend_str cov_data white_data percent_var n_var var_color
    clear C L i
    % Outputs:
    %       T_eig(i)               - Triangle eigenvalues
    %       T_normeig(i)           - Normalized triangle eigenvalues
    %       T_isncoll(i)           - List of triangles which are not collinear
%}
%% Create 1D array of valid triangles after filtering by collinearity
%{
    m_combix = cell(n_S,1); % index to valid combinations
    m_eig = cell(n_S,1); % 1d array of eigenvalues for valid combinations
    n_m = zeros(n_S,1); % number of valid combinations
    for i = 1:n_S;
        m_combix{i} = 1:n_comb(i);
        m_combix{i} = m_combix{i}(T_isncoll{i}); % Filter by collinearity
        m_eig{i} = T_eig{i}(m_combix{i},:); % Lowercase m denotes filtering by collinearity
        n_m(i) = numel(m_combix{i});
    end
    
    clear i
    % Outputs:
    %       m_combix                - [cleared next block] index of combinations which are not collinear
    %       m_eig                   - [cleared next block] Eigenvalues
    %       n_m                     - [cleared next block] number of valid combinations
%}
%% Find likely RANSAC pairs by filtering eigenvalues, radius
%{
    % Could also utilize 4-point combinations, graph theory, etc.
    % Also add radius information
    
%{
More intuitive way: to check consistency
i = 1; j = 2;
rep_raI = repmat(T_R{i}(:,1), [1,n_comb(j)]);
rep_raJ = repmat(T_R{j}(:,1)', [n_comb(i),1]);
error_ra = abs((rep_raI-rep_raJ)./((rep_raI + rep_raJ)/2));
is_ra = (error_ra < t_rad);
rep_rbI = repmat(T_R{i}(:,2), [1,n_comb(j)]);
rep_rbJ = repmat(T_R{j}(:,2)', [n_comb(i),1]);
error_rb = abs((rep_rbI-rep_rbJ)./((rep_rbI + rep_rbJ)/2));
is_rb = (error_rb < t_rad);
rep_rcI = repmat(T_R{i}(:,3), [1,n_comb(j)]);
rep_rcJ = repmat(T_R{j}(:,3)', [n_comb(i),1]);
error_rc = abs((rep_rcI-rep_rcJ)./((rep_rcI + rep_rcJ)/2));
is_rc = (error_rc < t_rad);
is_rad = is_ra & is_rb & is_rc;
%}
    
    M_bestixI = cell(n_S,n_S);
    M_bestixJ = cell(n_S,n_S);
    M_combI = cell(n_S,n_S);
    M_combJ = cell(n_S,n_S);
    M_eigerror = cell(n_S,n_S);
    n_M = zeros(n_S,n_S);
    for i = 1:n_S;
        for j = i+1:n_S; % Symmetric
            % Replicate eigenvalues to determine error
            M_eigaI = repmat(m_eig{i}(:,1), [1,n_m(j)]);
            M_eigbI = repmat(m_eig{i}(:,2), [1,n_m(j)]);
            M_eigaJ = repmat(m_eig{j}(:,1)', [n_m(i),1]);
            M_eigbJ = repmat(m_eig{j}(:,2)', [n_m(i),1]);
            M_error = (M_eigaI - M_eigaJ).^2 + (M_eigbI - M_eigbJ).^2;
            % Brute force removal
            is_brute = (M_error < t_eig_error);
            %xx denotes 1d array of valid
            xxM_error = M_error(is_brute);
            % Work in 1d instead of 2d matrix
            [row, col] = find(is_brute);
            xx_combixI = m_combix{i}(row)'; %1d array of all valid comb indices
            xx_combixJ = m_combix{j}(col)';
            xx_combI = T_comb{i}(xx_combixI,:);
            xx_combJ = T_comb{j}(xx_combixJ,:);
            n_brute = numel(row);
            % test xx_combI
%{
        temp_eig = zeros(n_brute,3);
        for t = 1:n_brute;
            temp_eig(t,:) = eig(cov(P_LCS{i}(:,xx_combI(t,:))'))';
        end
        temp_coll = temp_eig(:,2)./sum(temp_eig,2);
%}
            % Now remove all with insufficient radius similarity
            xx_israd = true(n_brute,1);
            for r = 1:3;
                %linear2 = sub2ind([n_tree(i), n_tree(j)], xxcombI(:,1), xxcombJ(:,1));
                linear = double(xx_combI(:,r)) + double((xx_combJ(:,r)-1))*n_tree(i);
                xx_israd = match_israd{i,j}(linear) & xx_israd;
            end
%{
        figure;
        subplot(1,2,1);imagesc(is_rad); title('Manual')
        subplot(1,2,2);imagesc(reshape(xx_israd,n_comb(1), n_comb(2))); title('Index');
%}
            xx_combixI = xx_combixI(xx_israd);
            xx_combixJ = xx_combixJ(xx_israd);
            %
            xx_combI = xx_combI(xx_israd,:);
            xx_combJ = xx_combJ(xx_israd,:);
            %
            xxM_error = xxM_error(xx_israd);
            % Now sort remaining
            % Sort eigenvalue error
            [xxM_error_sort,sortix] = sort(xxM_error);
            % Update indices of best triangles
            M_eigerror{i,j} = xxM_error_sort;
            M_bestixI{i,j} = xx_combixI(sortix);
            M_bestixJ{i,j} = xx_combixJ(sortix);
            M_combI{i,j} = xx_combI(sortix,:);
            M_combJ{i,j} = xx_combJ(sortix,:);
            n_M(i,j) = numel(sortix);
        end
    end
    
    % make n_M non directed
    n_M = n_M + n_M';
    
    clear match_rad_diff
    clear m_eig m_combix n_m
    clear M_eigaI M_eigaJ M_eigbI M_eigbJ M_error is_brute
    clear row col xx_combixI xx_combixJ xx_combI xx_combJ n_brute
    clear xx_israd linear xxM_error xxM_error_sort sortix r i j
    % M_bestixI{i,j}          - index of best i triangles for j-> i
    % M_bestixJ{i,j}          - index of best j triangles for j-> i
    % M_combI{i,j}            - triangle combinations
    % M_combJ{i,j}            - triangle combinations
    % M_eigerror{i,j}         - Eigenvalue error for each triangle match
%}
%% Construct ordered tri radius, locations for clarity
%{
    M_triI = cell(n_S,n_S);
    M_triJ = cell(n_S,n_S);
    M_triradI = cell(n_S,n_S);
    M_triradJ = cell(n_S,n_S);
    for i = 1:n_S;
        for j = i+1:n_S;
            minsize = min(n_M(i,j),t_RANSAC_nsearch);
            M_triI{i,j} = zeros(minsize,3,3);
            M_triI{i,j}(:,1,:) = P_LCS{i}(:,M_combI{i,j}(1:minsize,1))';
            M_triI{i,j}(:,2,:) = P_LCS{i}(:,M_combI{i,j}(1:minsize,2))';
            M_triI{i,j}(:,3,:) = P_LCS{i}(:,M_combI{i,j}(1:minsize,3))';
            M_triradI{i,j} = P_rad{i}(M_combI{i,j}(1:minsize,:));
            M_triJ{i,j} = zeros(minsize,3,3);
            M_triJ{i,j}(:,1,:) = P_LCS{j}(:,M_combJ{i,j}(1:minsize,1))';
            M_triJ{i,j}(:,2,:) = P_LCS{j}(:,M_combJ{i,j}(1:minsize,2))';
            M_triJ{i,j}(:,3,:) = P_LCS{j}(:,M_combJ{i,j}(1:minsize,3))';
            M_triradJ{i,j} = P_rad{j}(M_combJ{i,j}(1:minsize,:));
            M_eigerror{i,j} = M_eigerror{i,j}(1:minsize);
            M_combI{i,j}  = M_combI{i,j}(1:minsize,:);
            M_combJ{i,j}  = M_combJ{i,j}(1:minsize,:);
            M_bestixI{i,j}  = M_bestixI{i,j}(1:minsize);
            M_bestixJ{i,j}  = M_bestixJ{i,j}(1:minsize);
        end
    end
%{
i = 1;
j = 2;
M_raderror = abs((M_triradI{i,j}-M_triradJ{i,j})./((M_triradI{i,j}+M_triradJ{i,j})/2));
M_eigerror_check = zeros(n_M(i,j),1);
for t = 1:n_M(i,j);
    eI = eig(cov(squeeze(M_triI{i,j}(t,:,:))));
    eJ = eig(cov(squeeze(M_triJ{i,j}(t,:,:))));
    M_eigerror_check(t) = sum((eI-eJ).^2);
%{
    triI = squeeze(M_triI{i,j}(t,:,:));
    triJ = squeeze(M_triJ{i,j}(t,:,:));
    figure;
    plot_triangle(triI);
    plot_triangle(triJ);
    hold on; axis equal;
%}
end
%}
    
    clear i j
    % Outputs:
    % M_triI{i,j}             - triangle of points for set I
    % M_triJ{i,j}             - triangle of points for set J
    % M_triradI{i,j}          - radius of triangle of points for set I
    % M_triradJ{i,j}          - radius of triangle of points for set J
    %         save(sprintf('%sline637.mat', path_matvar));
    %     else
    %         load(sprintf('%sline637.mat',path_matvar));
    %     end
%}
%% Add additional permutations if similar radii
%{
    possible_is = [0 0 0; 1 0 0; 1 1 0; 1 1 1 ; 0 1 0; 0 1 1; 0 0 1];
    n_possible = size(possible_is,1);
    for i = 1:n_S;
        for j = i+1:n_S;
            rad_diff1 = abs((M_triradI{i,j} - circshift(M_triradJ{i,j}, [0 -1]))./...
                ((M_triradI{i,j} + circshift(M_triradJ{i,j}, [0 -1]))/2));
            rad_diff2 = abs((M_triradJ{i,j} - circshift(M_triradI{i,j}, [0 -1]))./...
                ((M_triradJ{i,j} + circshift(M_triradI{i,j}, [0 -1]))/2));
            is_rad1 = (rad_diff1 < t_rad);
            is_rad2 = (rad_diff2 < t_rad);
            I1J2 = is_rad1(:,1);
            I2J3 = is_rad1(:,2);
            I2J1 = is_rad2(:,1);
            I3J2 = is_rad2(:,2);
            I1J2andI2J1 = I1J2&I2J1;
            I2J3andI3J2 = I2J3&I3J2;
            %I1J2andI2J1andI2J3andI3J2 = I1J2andI2J1&I2J3andI3J2 ;
            M_triI{i,j} = cat(1,M_triI{i,j}, M_triI{i,j}(I1J2andI2J1,[1 2 3],:),...
                M_triI{i,j}(I2J3andI3J2,[1 2 3],:));
            M_triJ{i,j} = cat(1,M_triJ{i,j}, M_triJ{i,j}(I1J2andI2J1,[2 1 3],:),...
                M_triJ{i,j}(I2J3andI3J2,[1 3 2],:));
            M_triradI{i,j} = cat(1,M_triradI{i,j}, M_triradI{i,j}(I1J2andI2J1,[1 2 3]),...
                M_triradI{i,j}(I2J3andI3J2,[1 2 3]));
            M_triradJ{i,j} = cat(1,M_triradJ{i,j}, M_triradJ{i,j}(I1J2andI2J1,[2 1 3]),...
                M_triradJ{i,j}(I2J3andI3J2,[1 3 2]));
            M_eigerror{i,j} = cat(1,M_eigerror{i,j}, M_eigerror{i,j}(I1J2andI2J1),...
                M_eigerror{i,j}(I2J3andI3J2));
            M_combI{i,j} = cat(1,M_combI{i,j}, M_combI{i,j}(I1J2andI2J1,[ 1 2 3]),...
                M_combI{i,j}(I2J3andI3J2,[1 2 3]));
            M_combJ{i,j} = cat(1,M_combJ{i,j}, M_combJ{i,j}(I1J2andI2J1, [2 1 3]),...
                M_combJ{i,j}(I2J3andI3J2, [1 3 2]));
            M_bestixI{i,j} = cat(1,M_bestixI{i,j}, M_bestixI{i,j}(I1J2andI2J1),...
                M_bestixI{i,j}(I2J3andI3J2));
            M_bestixJ{i,j} = cat(1,M_bestixJ{i,j}, M_bestixJ{i,j}(I1J2andI2J1),...
                M_bestixJ{i,j}(I2J3andI3J2));
            [M_eigerror{i,j}, sortix] = sort( M_eigerror{i,j});
            M_triI{i,j} = M_triI{i,j}(sortix,:,:);
            M_triJ{i,j} = M_triJ{i,j}(sortix,:,:);
            M_triradI{i,j} = M_triradI{i,j}(sortix,:);
            M_triradJ{i,j} = M_triradJ{i,j}(sortix,:);
            M_combI{i,j} = M_combI{i,j}(sortix,:);
            M_combJ{i,j} = M_combJ{i,j}(sortix,:);
            M_bestixI{i,j} = M_bestixI{i,j}(sortix);
            M_bestixJ{i,j} = M_bestixJ{i,j}(sortix);
        end
    end
    n_M = cellfun(@numel, M_eigerror);
%}
%% OLD CODE Unnecessary figures
%{
% Both are in local coordinate system
%{
t = 1; i = 1; j = 2;
figure; hold on; axis equal;
for p= 1:n_tree(i);
    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
end
for p= 1:n_tree(j);
    filledCircle([P_LCS{j}(1,p); P_LCS{j}(2,p)]',P_rad{j}(p),1000,P_color(j,:));
end
plot_triangle(squeeze(M_triI{i,j}(t,:,:)));
plot_triangle(squeeze(M_triJ{i,j}(t,:,:)));
plot_triangle_connection(squeeze(M_triI{i,j}(t,:,:)),squeeze(M_triJ{i,j}(t,:,:)));
%}
% Check first few and make sure the data are correct % above plot is better
%{
i = 1;
j = 2;
r = 1;
for r = 1:1:10;
figure('position', [76, 354, 1498, 585]);
subplot(1,3,1);
hold on
hdummy = plot3(P_LCS{i}(1,1),P_LCS{i}(2,1),P_LCS{i}(3,1),'ok','markersize',5,...
    'markerfacecolor',P_color(i,:));
set(hdummy, 'visible', 'off');
triI3 =  squeeze(T_LCS{i}(M_bestixI{i,j}(r),:,:));
triI4 = [triI3; triI3(1,:)];
plot(triI4(:,1),triI4(:,2),'color', P_color(i,:))
for t = 1:numel(P_rad{i});
    h = filledCircle([P_LCS{i}(1,t); P_LCS{i}(2,t)]',P_rad{i}(t),1000,P_color(i,:));
end
axis auto
axisval = axis;
xlabel('x Position relative [m]');
ylabel('y Position relative [m]');
zlabel('z Position relative [m]');
view(0,90);
grid on
%titlestr = sprintf('Scan %g',P_plot(s));
%title(titlestr);
legend_str{1} = sprintf('Stem map I');
legend_str{2} = sprintf('Triangle');%, P_plot(s));
legend(legend_str,'location','northeast');
legend boxoff       % Hides the legend's axes (legend border and background)
subplot(1,3,2);
hold on
hdummy = plot3(P_LCS{j}(1,1),P_LCS{j}(2,1),P_LCS{j}(3,1),'ok','markersize',5,...
    'markerfacecolor',P_color(j,:));
set(hdummy, 'visible', 'off');
triJ3 =  squeeze(T_LCS{j}(M_bestixJ{i,j}(r),:,:));
triJ4 = [triJ3; triJ3(1,:)];
plot(triJ4(:,1),triJ4(:,2),'color', P_color(j,:))
for t = 1:numel(P_rad{j});
    h = filledCircle([P_LCS{j}(1,t); P_LCS{j}(2,t)]',P_rad{j}(t),1000,P_color(j,:));
end
axis auto
axisval = axis;
xlabel('x Position relative [m]');
ylabel('y Position relative [m]');
zlabel('z Position relative [m]');
view(0,90);
grid on
%titlestr = sprintf('Scan %g',P_plot(s));
%title(titlestr);
legend_str{1} = sprintf('Stem map J');
legend_str{2} = sprintf('Triangle');%, P_plot(s));
legend(legend_str,'location','northeast');
legend boxoff       % Hides the legend's axes (legend border and background)
subplot(1,3,3);
hold on
tI = M_bestixI{i,j}(r);
tJ = M_bestixJ{i,j}(r);
%eig(cov(triI3));
%eig(cov(triJ3));
[Rhat,that] = rigid_transform_3D(triJ3,triI3);
dataJ = P_LCS{j}';
dataJt = ((Rhat*dataJ')+ repmat(that,1,size(dataJ,1))); %Hard code 3
triJt = ((Rhat*triJ4')+ repmat(that,1,size(triJ4,1)));
triJt = triJt';
hdummy = plot3(P_LCS{i}(1,P_ix1{i}),P_LCS{i}(2,P_ix1{i}),P_LCS{i}(3,P_ix1{i}),'ok','markersize',5,...
    'markerfacecolor',P_color(i,:));
set(hdummy, 'visible', 'off');
hdummy = plot3(P_LCS{j}(1,P_ix1{j}),P_LCS{j}(2,P_ix1{j}),P_LCS{j}(3,P_ix1{j}),'ok','markersize',5,...
    'markerfacecolor',P_color(j,:));
set(hdummy, 'visible', 'off');
plot(triI4(:,1),triI4(:,2),'color',P_color(i,:))
plot(triJt(:,1),triJt(:,2),'color',P_color(j,:))
for t = 1:numel(P_rad{i});
    h = filledCircle([P_LCS{i}(1,t); P_LCS{i}(2,t)]',P_rad{i}(t),1000,P_color(i,:));
end
for t = 1:numel(P_rad{j});
    h = filledCircle([dataJt(1,t); dataJt(2,t)]',P_rad{j}(t),1000,P_color(j,:));
end
axis auto
axisval = axis;
xlabel('x Position relative [m]');
ylabel('y Position relative [m]');
zlabel('z Position relative [m]');
view(0,90);
grid on
%titlestr = sprintf('Scan %g',P_plot(s));
%title(titlestr);
legend_str{1} = sprintf('Stem map I');
legend_str{1} = sprintf('Stem map J');
legend_str{2} = sprintf('TriangleI');%, P_plot(s));
legend_str{2} = sprintf('TriangleJ');%, P_plot(s));
legend(legend_str,'location','northeast');
legend boxoff       % Hides the legend's axes (legend border and background)
end
%}
%}
%% OLD CODE Permutations??
%{
perm_of_3 = perms(1:3);
n_perm_of_3 = size(perm_of_3,1);
perm_of_3a = perm_of_3;
perm_of_3b = repmat(perm_of_3(1,:),[n_perm_of_3,1]);
n_p = size(perm_of_3a,1);
%}
%% Find pairwise matches between full sets
%{
    fprintf('\nFind pairwise matches \n');
    
    tree_ix = cell(n_S,1);
    for s = 1:n_S;
        tree_ix{s} = 1:n_tree(s);
    end
    
    % Initialize graph and corresponding functional relationships
    match_i = cell(n_S); %base
    match_j = cell(n_S); %mobile
    match_R = cell(n_S);
    match_t = cell(n_S);
    match_nit = nan(n_S);
    match1_R = cell(n_S);
    match1_t = cell(n_S);
    P_ixI1 = cell(n_S);
    P_ixJ1 = cell(n_S);
    P_ixI2 = cell(n_S);
    P_ixJ2 = cell(n_S);
    triI_best = [];
    % Find rotation and tranlation from j to i along with point matching pairs
    for i = 1:n_S-1;
        for j = i+1:n_S; % Symmetric
            fprintf('\n\tMatching %g to %g\n',j,i);
            if i==j;
                continue
            end
            % RANSAC setup
            n_match_best = 0;
            m1_best = [];
            m2_best = [];
            n_search = n_M(i,j);% min(t_RANSAC_nsearch, n_M(i,j));
            
            for r =1:n_search;
                
                triI = squeeze(M_triI{i,j}(r,:,:));
                triJ = squeeze(M_triJ{i,j}(r,:,:));
                
                % Check input Plotting
%{
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    h = filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,'r');
                    set(h, 'FaceAlpha', 0.5);
                end
                for p= 1:n_tree(j);
                    h = filledCircle([P_LCS{j}(1,p); P_LCS{j}(2,p)]',P_rad{j}(p),1000,'b');
                    set(h, 'FaceAlpha', 0.5)
                end
                if ~isempty(triI_best);
                    plot(triI_best(:,1), triI_best(:,2), 'xk', 'markersize', 10, 'linewidth', 2);
                    plot(triJ_best(:,1), triJ_best(:,2), 'xk', 'markersize', 10, 'linewidth', 2);
                end
                plot_triangle(triI);
                plot_triangle(triJ);
                plot_triangle_connection(triI,triJ);
                foo = 1;
%}
                
                [Rhat,that] = rigid_transform_3D(triJ,triI);
                dataJt = ((Rhat*P_LCS{j})+ repmat(that,1,n_tree(j)));
                % Add rotation information
                % Use n_match, round-trip or other method to determine fit
                P_LCSIrep = repmat(P_LCS{i}',[1,1,n_tree(j)]);
                P_LCSJtrep = zeros(n_tree(i),3,n_tree(j));
                for t = 1:n_tree(i);
                    P_LCSJtrep(t,:,:) = dataJt;
                end
                error_xyz = squeeze(sum((P_LCSIrep-P_LCSJtrep).^2,2));
                
                [m1,m2] = find(error_xyz<t_RANSAC_xyz& match_israd{i,j});
                n_match = numel(m1);
                
                % Check error calculation
%{
                triJt = ((Rhat*triJ')+ repmat(that,1,3))';
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
                end
                for p= 1:n_tree(j);
                    filledCircle([dataJt(1,p); dataJt(2,p)]',P_rad{j}(p),1000,P_color(j,:));
                end
                dist_color = plyintensity2color( error_xyz(1,:), 'jet' )./255;
                for c = 1:n_tree(j);
                    plot([P_LCS{i}(1,1), dataJt(1,c)],...
                        [P_LCS{i}(2,1), dataJt(2,c)],...
                        'color', dist_color(c,:), 'linewidth', 2);
                end
%}
                
                if n_match > n_match_best;
                    n_match_best = n_match;
                    R_best = Rhat;
                    t_best = that;
                    triI_best = triI;
                    triJ_best = triJ;
                    m1_best = m1;
                    m2_best = m2;
                    nit_best = r;
                end
                
                % Check output Plotting
%{
                dataJt = ((Rhat*P_LCS{j})+ repmat(that,1,n_tree(j)));
                triJt = ((Rhat*triJ')+ repmat(that,1,3))';
                figure; hold on; axis equal;
                for m = 1:numel(m1_best);
                    m1 = m1_best(m);
                    m2 = m2_best(m);
                    
                    h = filledCircle([(P_LCS{i}(1,m1) +dataJt(1,m2))/2;...
                         (P_LCS{i}(2,m1)+ dataJt(2,m2))/2]',3*P_rad{i}(m1),1000,'y');
                    set(h, 'FaceAlpha', 0.4);
                end
                for p= 1:n_tree(i);
                    h = filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,'r');
                    set(h, 'FaceAlpha', 0.4);
                end
                for p= 1:n_tree(j);
                    h = filledCircle([dataJt(1,p); dataJt(2,p)]',P_rad{j}(p),1000,'b');
                    set(h, 'FaceAlpha', 0.4);
                end
                plot_triangle(triI);
                plot_triangle(triJt);
                plot_triangle_connection(triI,triJt);
%}
            end
            
            if n_match_best > 3;
                % Update disjoint sets
                r = nit_best;
                P_ixI1{i,j} = M_combI{i,j}(r,:,:);
                P_ixJ1{i,j} = M_combJ{i,j}(r,:,:);
                P_ixI2{i,j} = setdiff(tree_ix{i}, P_ixI1{i,j});
                P_ixJ2{i,j} = setdiff(tree_ix{j}, P_ixJ1{i,j});
                P_ixI1{j,i} = M_combJ{i,j}(r,:,:);
                P_ixJ1{j,i} = M_combI{i,j}(r,:,:);
                P_ixI2{j,i} = setdiff(tree_ix{j}, P_ixI1{j,i});
                P_ixJ2{j,i} = setdiff(tree_ix{i}, P_ixJ1{j,i});
                match1_R{i,j} = R_best;
                match1_t{i,j} = t_best;
                % Re-estimate model with inliers
                [R_reest,t_reest] = rigid_transform_3D(P_LCS{j}(:,m2_best)',P_LCS{i}(:,m1_best)');
                match_R{i,j} = R_reest;
                match_t{i,j} = t_reest;
                match_i{i,j} = m1_best;
                match_j{i,j} = m2_best;
                match_nit(i,j) = nit_best;
                % Check output Plotting
%{
            P_LCSI1 = P_LCS{i}(:, P_ixI1{i,j});
            P_LCSJ1 = P_LCS{j}(:, P_ixJ1{i,j});
            P_LCSI2 = P_LCS{i}(:, P_ixI2{i,j});
            P_LCSJ2 = P_LCS{j}(:, P_ixJ2{i,j});
            P_radI1 = P_rad{i}( P_ixI1{i,j});
            P_radJ1 = P_rad{j}( P_ixJ1{i,j});
            P_radI2 = P_rad{i}( P_ixI2{i,j});
            P_radJ2 = P_rad{j}( P_ixJ2{i,j});
            r = match_nit(i,j);
            triI = squeeze(M_triI{i,j}(r,:,:));
            triJ = squeeze(M_triJ{i,j}(r,:,:));
            figure; hold on; axis equal;
            for p= 1:3;
                h = filledCircle([P_LCSI1(1,p); P_LCSI1(2,p)]',P_radI1(p),1000,'r');
                set(h, 'FaceAlpha', 0.5);
            end
            for p= 1:3;
                h = filledCircle([P_LCSJ1(1,p); P_LCSJ1(2,p)]',P_radJ1(p),1000,'b');
                set(h, 'FaceAlpha', 0.5)
            end
            for p= 1:n_tree(i)-3;
                h = filledCircle([P_LCSI2(1,p); P_LCSI2(2,p)]',P_radI2(p),1000,'r');
                set(h, 'FaceAlpha', 0.3);
            end
            for p= 1:n_tree(j)-3;
                h = filledCircle([P_LCSJ2(1,p); P_LCSJ2(2,p)]',P_radJ2(p),1000,'b');
                set(h, 'FaceAlpha', 0.3)
            end
            plot(P_LCSI1(1,:), P_LCSI1(2,:), 'xk', 'markersize', 10, 'linewidth', 2);
            plot(P_LCSJ1(1,:), P_LCSJ1(2,:), 'xk', 'markersize', 10, 'linewidth', 2);

%{
                for p= 1:n_tree(i);
                    h = filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,'r');
                    set(h, 'FaceAlpha', 0.5);
                end
                for p= 1:n_tree(j);
                    h = filledCircle([P_LCS{j}(1,p); P_LCS{j}(2,p)]',P_rad{j}(p),1000,'b');
                    set(h, 'FaceAlpha', 0.5)
                end
%}
            plot_triangle(triI);
            plot_triangle(triJ);
            plot_triangle_connection(triI,triJ);
            foo = 1;
%}
            end
        end
    end
    n_ix1 = cellfun(@numel, P_ixI1);
    n_ix2 = cellfun(@numel, P_ixI2);
    
    % Plot first 10 matches
%{
    colortri = jet(10);
    figure; hold on; axis equal;
    i = 1; j = 2;
    for p= 1:n_tree(i);
        h = filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,'r');
        set(h, 'FaceAlpha', 0.5);
    end
    for p= 1:n_tree(j);
        h = filledCircle([P_LCS{j}(1,p); P_LCS{j}(2,p)]',P_rad{j}(p),1000,'b');
        set(h, 'FaceAlpha', 0.5)
    end
    for r =1:10;
        triI = squeeze(M_triI{i,j}(r,:,:));
        triJ = squeeze(M_triJ{i,j}(r,:,:));
        h = plot_triangle(triI); set(h, 'color', colortri(r,:));
        h = plot_triangle(triJ); set(h, 'color', colortri(r,:));
        plot_triangle_connection(triI,triJ);
    end
    % Cross out best points
    r = match_nit(i,j);
    triI = squeeze(M_triI{i,j}(r,:,:));
    triJ = squeeze(M_triJ{i,j}(r,:,:));
    plot(triI(:,1), triI(:,2), 'xk', 'markersize', 10, 'linewidth', 2);
    plot(triJ(:,1), triJ(:,2), 'xk', 'markersize', 10, 'linewidth', 2);
%}
    
    % Declare match from i-i identity
    for i = 1:n_S;
        match_R{i,i} = eye(3);
        match_t{i,i} = zeros(3,1);
        match_i{i,i} = 1:numel(P_LCS{i});
        match_j{i,i} = 1:numel(P_LCS{i});
    end
    
    % Make R,t non-directed
    %
    for i = 1:n_S;
        for j = i:n_S;
            if ~isempty(match_R{i,j});
                match_R{j,i} = match_R{i,j}';
                match_t{j,i} = -(match_R{i,j}')*match_t{i,j};
            end
        end
    end
    %
    
    % Decompose for easy visualiztion
    match_rx = 180*cellfun(@decompose_rotation_rx,match_R)/pi;
    match_ry = 180*cellfun(@decompose_rotation_ry,match_R)/pi;
    match_rz = 180*cellfun(@decompose_rotation_rz,match_R)/pi;
    match_tx = cellfun(@decompose_translation_tx,match_t);
    match_ty = cellfun(@decompose_translation_ty,match_t);
    match_tz = cellfun(@decompose_translation_tz,match_t);
    
    % Declare match from i-i identity
    for i = 1:n_S;
        match1_R{i,i} = eye(3);
        match1_t{i,i} = zeros(3,1);
        match1_i{i,i} = 1:numel(P_LCS{i});
        match1_j{i,i} = 1:numel(P_LCS{i});
    end
    
    % Make R,t non-directed
    %
    for i = 1:n_S;
        for j = i:n_S;
            if ~isempty(match1_R{i,j});
                match1_R{j,i} = match1_R{i,j}';
                match1_t{j,i} = -(match1_R{i,j}')*match1_t{i,j};
            end
        end
    end
    %
    
    match1_rx = 180*cellfun(@decompose_rotation_rx,match1_R)/pi;
    match1_ry = 180*cellfun(@decompose_rotation_ry,match1_R)/pi;
    match1_rz = 180*cellfun(@decompose_rotation_rz,match1_R)/pi;
    match1_tx = cellfun(@decompose_translation_tx,match1_t);
    match1_ty = cellfun(@decompose_translation_ty,match1_t);
    match1_tz = cellfun(@decompose_translation_tz,match1_t);
    
    filepath_nit = sprintf('%s%s',path_mat, 'match_nit.mat');
    save(filepath_nit, 'match_nit');
    clear n_match_best m1_best m2_best n_search i j r Rhat that dataJt
    clear P_LCSIrep P_LCSJtrep t error_xyz m1 m2 n_match
    % Outputs
    % match_R               {i,j} is pairwise rotation from j into i
    % match_t               {i,j} is pairwise translation from j into i
    % match_i               {i,j} is pairwise matches from i
    % match_j               {i,j} is pairwise matches from j
    % match_rx              {i,j} is pairwise rotation x from j into i
    % match_ry              {i,j} is pairwise rotation y from j into i
    % match_rz              {i,j} is pairwise rotation z from j into i
%}
%% JUST PLOTTING Plot pose for each i,j
%{
    x_axis = [0 0 0; 1 0 0]';
    y_axis = [0 0 0; 0 1 0]';
    z_axis = [0 0 0; 0 0 1]';
    for i= 1:n_S;
        figure;
        title(sprintf('Pairwise Pose Estimates for Site %3.0f-%g', info_site,i));
        hold on;
        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis equal
        for j = 1:n_S;
            x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
            y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
            z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
        end
    end
%}
%% Update LCS and rad arrays for disjoint sets
%{
    P_LCSI1 = cell(n_S);
    P_LCSI2 = cell(n_S);
    P_LCSJ1 = cell(n_S);
    P_LCSJ2 = cell(n_S);
    P_radI1 = cell(n_S);
    P_radI2 = cell(n_S);
    P_radJ1 = cell(n_S);
    P_radJ2 = cell(n_S);
    for i = 1:n_S;
        for j = 1:n_S;
            P_LCSI1{i,j} = P_LCS{i}(:, P_ixI1{i,j});
            P_LCSJ1{i,j} = P_LCS{j}(:, P_ixJ1{i,j});
            P_LCSI2{i,j} = P_LCS{i}(:, P_ixI2{i,j});
            P_LCSJ2{i,j} = P_LCS{j}(:, P_ixJ2{i,j});
            P_radI1{i,j} = P_rad{i}( P_ixI1{i,j});
            P_radJ1{i,j} = P_rad{j}( P_ixJ1{i,j});
            P_radI2{i,j} = P_rad{i}( P_ixI2{i,j});
            P_radJ2{i,j} = P_rad{j}( P_ixJ2{i,j});
        end
    end
    n2_treeI = cellfun(@numel, P_radI2);
    n2_treeJ = cellfun(@numel, P_radJ2);
    
%{
i = 1; j = 2;
figure; hold on; axis equal;
for p= 1:3;
    h = filledCircle([P_LCSI1{i,j}(1,p); P_LCSI1{i,j}(2,p)]',P_radI1{i,j}(p),1000,'r');
    set(h, 'FaceAlpha', 0.5);
end
for p= 1:3;
    h = filledCircle([P_LCSJ1{i,j}(1,p); P_LCSJ1{i,j}(2,p)]',P_radJ1{i,j}(p),1000,'b');
    set(h, 'FaceAlpha', 0.5)
end
for p= 1:n_tree(i)-3;
    h = filledCircle([P_LCSI2{i,j}(1,p); P_LCSI2{i,j}(2,p)]',P_radI2{i,j}(p),1000,'r');
    set(h, 'FaceAlpha', 0.3);
end
for p= 1:n_tree(j)-3;
    h = filledCircle([P_LCSJ2{i,j}(1,p); P_LCSJ2{i,j}(2,p)]',P_radJ2{i,j}(p),1000,'b');
    set(h, 'FaceAlpha', 0.3)
end
plot(P_LCSI1{i,j}(1,:), P_LCSI1{i,j}(2,:), 'xk', 'markersize', 10, 'linewidth', 2);
plot(P_LCSJ1{i,j}(1,:), P_LCSJ1{i,j}(2,:), 'xk', 'markersize', 10, 'linewidth', 2);
foo = 1;
%}
%}
%% Start from scratch
%{
M2_triI = cell(n_S,n_S);
M2_triJ = cell(n_S,n_S);
M2_triradI = cell(n_S,n_S);
M2_triradJ = cell(n_S,n_S);
for i = 1:n_S;
    for j = i+1:n_S;
        is_validI = ~any(ismember(M_combI{i,j}, P_ixI1{i,j}),2);
        is_validJ = ~any(ismember(M_combJ{i,j}, P_ixJ1{i,j}),2);
        is_validIJ = is_validI&is_validJ;
        M2_triI{i,j} = M_triI{i,j}(is_validIJ,:,:);
        M2_triJ{i,j} = M_triJ{i,j}(is_validIJ,:,:);
        M2_triradI{i,j} = M_triradI{i,j}(is_validIJ,:);
        M2_triradJ{i,j} = M_triradJ{i,j}(is_validIJ,:);
    end
end
n_MI2 = cellfun(@(x) size(x,1), M2_triradI);
n_MJ2 = cellfun(@(x) size(x,1), M2_triradJ);
% Make n_MJ2 nondirected
n_MI2 = n_MI2 + n_MI2';
n_MJ2 = n_MJ2 + n_MJ2';

% Outputs:
% M2_triI{i,j}             - triangle of points for set 1I
% M2_triJ{i,j}             - triangle of points for set 1J
% M2_triradI{i,j}          - radius of triangle of points for set 1I
% M2_triradJ{i,j}          - radius of triangle of points for set 1J
%}
%% Find pairwise matches between disjoint sets for set 2
%{
fprintf('\nFind pairwise matches \n');

% Initialize graph and corresponding functional relationships
match2_i = cell(n_S); %base
match2_j = cell(n_S); %mobile
match2_R = cell(n_S);
match2_t = cell(n_S);

% Find rotation and tranlation from j to i along with point matching pairs
for i = 1:n_S; %why -1 before?
    for j = i+1:n_S; % Symmetric
        fprintf('\n\tMatching %g to %g\n',j,i);
        
        if i==j;
            continue
        end
        % RANSAC setup
        n_match_best = 0;
        m1_best = [];
        m2_best = [];
        
        n_search = min([n_MJ2(i,j),t_RANSAC_nsearch, n_MI2(i,j)]);
        
        
        for r = 1:n_search;
            
            triI = squeeze(M2_triI{i,j}(r,:,:));
            triJ = squeeze(M2_triJ{i,j}(r,:,:));
            
            % Check input Plotting
%{
            figure; hold on; axis equal;
            for p= 1:n2_treeI(i,j);
                h = filledCircle([P_LCSI2{i,j}(1,p); P_LCSI2{i,j}(2,p)]',P_radI2{i,j}(p),1000,'r');
                set(h, 'FaceAlpha', 0.4)
            end
            for p= 1:n2_treeJ(i,j);
                h = filledCircle([P_LCSJ2{i,j}(1,p); P_LCSJ2{i,j}(2,p)]',P_radJ2{i,j}(p),1000,'b');
                set(h, 'FaceAlpha', 0.4)
            end
            plot(P_LCSI1{i,j}(1,:), P_LCSI1{i,j}(2,:), 'xk', 'markersize', 10, 'linewidth', 2);
            plot(P_LCSJ1{i,j}(1,:), P_LCSJ1{i,j}(2,:), 'xk', 'markersize', 10, 'linewidth', 2);
            h = plot_triangle(P_LCSI1{i,j}'); set(h, 'color', 'k');
            h = plot_triangle(P_LCSJ1{i,j}'); set(h, 'color', 'k');
             plot_triangle_connection(P_LCSI1{i,j}',P_LCSJ1{i,j}');
            h = plot_triangle(triI); set(h, 'color', 'r');
            h = plot_triangle(triJ); set(h, 'color', 'b');
            plot_triangle_connection(triI,triJ);
            foo = 1;
%}
            
            [Rhat,that] = rigid_transform_3D(triJ,triI);
            dataJt = ((Rhat*P_LCSJ2{i,j})+ repmat(that,1,n2_treeJ(i,j)));
            % Add rotation information
            % Use n_match, round-trip or other method to determine fit
            P_LCSIrep = repmat(P_LCSI2{i,j}',[1,1,n2_treeJ(i,j)]);
            P_LCSJtrep = zeros(n2_treeI(i,j),3,n2_treeJ(i,j));
            for t = 1:n2_treeI(i,j);
                P_LCSJtrep(t,:,:) = dataJt;
            end
            error_xyz = squeeze(sum((P_LCSIrep-P_LCSJtrep).^2,2));
            [m1,m2] = find(error_xyz<t_RANSAC_xyz);
            n_match = numel(m1);
            
            % Check error calculation
%{
                triJt = ((Rhat*triJ')+ repmat(that,1,3))';
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
                end
                for p= 1:n_tree(j);
                    filledCircle([dataJt(1,p); dataJt(2,p)]',P_rad{j}(p),1000,P_color(j,:));
                end
                dist_color = plyintensity2color( error_xyz(1,:), 'jet' )./255;
                for c = 1:n_tree(j);
                    plot([P_LCS{i}(1,1), dataJt(1,c)],...
                        [P_LCS{i}(2,1), dataJt(2,c)],...
                        'color', dist_color(c,:), 'linewidth', 2);
                end
%}
            
            if n_match > n_match_best;
                n_match_best = n_match;
                R_best = Rhat;
                t_best = that;
                m1_best = m1;
                m2_best = m2;
            end
            
            % Check output Plotting
%{
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
                end
                for p= 1:n_tree(j);
                    filledCircle([dataJt(1,p); dataJt(2,p)]',P_rad{j}(p),1000,P_color(j,:));
                endRb
                plot_triangle(triI);
                plot_triangle(triJt);
                plot_triangle_connection(triI,triJt);
%}
        end
        
        if n_match_best >= 3;
            % Re-estimate model with inliers
            [R_reest,t_reest] = rigid_transform_3D(P_LCSJ2{i,j}(:,m2_best)',P_LCSI2{i,j}(:,m1_best)');
            match2_R{i,j} = R_reest;
            match2_t{i,j} = t_reest;
            match2_i{i,j} = m1_best;
            match2_j{i,j} = m2_best;
            rnew =  match2_R{i,j}';
            tnew = -(match2_R{i,j}')*match2_t{i,j};
            match12_R = rnew*match1_R{i,j};
            match12_t = (rnew*match1_t{i,j}) + tnew;
            P_LCSt =  (match12_R*P_LCS{i})+ repmat(match12_t,1,n_tree(i));
            P_dist = sqrt(sum((P_LCSt - P_LCS{i}).^2));
            match12_RMSE = nanmean(P_dist);
            fprintf('\n %d -> %d: RMSE = %4.3f\n', j, i,match12_RMSE);
            foo =1;
        else
            foo = 1;
        end
    end
end

% Declare match from i-i identity
for i = 1:n_S;
    match2_R{i,i} = eye(3);
    match2_t{i,i} = zeros(3,1);
    match2_i{i,i} = 1:numel(P_LCS{i});
    match2_j{i,i} = 1:numel(P_LCS{i});
end

% Make R,t non-directed
%
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(match2_R{i,j});
            match2_R{j,i} = match2_R{i,j}';
            match2_t{j,i} = -(match2_R{i,j}')*match2_t{i,j};
        end
    end
end
%

% Flip R,t to go the opposite direction (J->I)
match2_Rflip = cell(n_S,n_S);
match2_tflip = cell(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        match2_Rflip{i,j} = match2_R{i,j}';
        match2_tflip{i,j} = -(match2_R{i,j}')*match2_t{i,j};
    end
end
match2_R = match2_Rflip;
match2_t = match2_tflip;

match2_rx = 180*cellfun(@decompose_rotation_rx,match2_R)/pi;
match2_ry = 180*cellfun(@decompose_rotation_ry,match2_R)/pi;
match2_rz = 180*cellfun(@decompose_rotation_rz,match2_R)/pi;
match2_tx = cellfun(@decompose_translation_tx,match2_t);
match2_ty = cellfun(@decompose_translation_ty,match2_t);
match2_tz = cellfun(@decompose_translation_tz,match2_t);

clear match2_tflip match2_Rflip
clear filepath_match_R filepath_match_t filepath_match_i filepath_match_j
clear n_match_best Rhat that m1_best m2_best n_search
clear r i j h n_match_best triI triJ Rhat that triJt dataJt
clear  P_LCSIrep t P_LCSJtrep  error_xyz m1 m2 n_match
clear R_best t_best
% Outputs
% match2_R               (i,j) is pairwise rotation from j into i
% match2_t               (i,j) is pairwise translation from j into i
% match2_i               (i,j) is pairwise matches from i
% match2_j               (i,j) is pairwise matches from j
%}
%% Determine RMSE error of stem centers
%{
isempty1 = cellfun(@isempty, match1_R);
isempty2 = cellfun(@isempty, match2_R);
match12_R = cell(n_S,n_S);
match12_t = cell(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        if ~isempty1(i,j) && ~isempty2(i,j);
            match12_R{i,j} = match2_R{i,j}*match1_R{i,j};
            match12_t{i,j} = (match2_R{i,j}*match1_t{i,j}) + match2_t{i,j};
        end
    end
end

%match12_R = cellfun(@(x,y) x*y,match2_R, match1_R, 'uniformOutput', false);
%match12_t = cellfun(@(x,y,z) (x*y) + z,match2_R, match1_t, match2_t, 'uniformOutput', false);
match12_rx = 180*cellfun(@decompose_rotation_rx,match12_R)/pi;
match12_ry = 180*cellfun(@decompose_rotation_ry,match12_R)/pi;
match12_rz = 180*cellfun(@decompose_rotation_rz,match12_R)/pi;
match12_er = abs(match12_rx) + abs(match12_ry) + abs(match12_rz);
%match12_et = cellfun(@sum,cellfun(@abs,match12_t,'uniformOutput',false));
%match12_etxy = cellfun(@sum,cellfun(@abs,cellfun(@(x) x(1:2),match12_t, 'uniformoutput', false),'uniformOutput',false));

% OR, could generate synthetic points about volume and R1t1R2t2 them
match12_RMSE = nan(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        if ~isnan(match12_er(i,j));
            P_LCSt =  (match12_R{i,j}*P_LCS{i})+ repmat(match12_t{i,j},1,n_tree(i));
            match12_RMSE(i,j) = sqrt(nanmean(sum((P_LCSt - P_LCS{i}).^2)));
        end
    end
end
%}
%% Flag one way transformations which are not similar to full set
%{
match_rxdiff = abs(match_rx - match1_rx);
match_rydiff = abs(match_ry - match1_ry);
match_rzdiff = abs(match_rz - match1_rz);
match_rdiff = (match_rxdiff + match_rydiff + match_rzdiff)/3;
match_rdiffmax = max(cat(3,match_rxdiff , match_rydiff , match_rzdiff),[],3);

isempty0 = cellfun(@isempty, match_R);
match_tdiff = zeros(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        if ~isempty1(i,j) && ~isempty0(i,j);
            match_tdiff(i,j) = sqrt(sum(match_t{i,j} - match1_t{i,j}).^2);
        end
    end
end
%  match_tdiff = cellfun(@(x,y) sqrt(sum((x-y).^2)), match_t, match1_t);
match_flag = (match_rxdiff>t_flagr) | (match_rydiff>t_flagr) | ...
    (match_rzdiff>t_flagr) | (match_tdiff > t_flagt);

%abs(match_t - match1_t);
%match_tydiff = abs(match_ty - match1_ty);
%match_tzdiff = abs(match_tz - match1_tz);
%}
%% Plot all pairwise sensors poses relative to i
%{
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
[x_sph, y_sph, z_sph] = sphere;

for i= 14;%:n_S;
    err_alpha = match12_RMSE(i,:);
    err_color = (double(vec2cmap(err_alpha, 'jet', 0,2 )))./255;
    figure;
    title(sprintf('Pairwise Pose Estimates for Site %3.0f-%g', info_site,i));
    hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    for j = 1:n_S;
        if ~isempty0(i,j);
            x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
            y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
            z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
                z_sph+match_t{i,j}(3));
            alpha(0.2)
            set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
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
    set(gcf, 'color', 'white')
    filepath_fig = sprintf('%spose_c2_%03.0f-%02.0f.png',...
        path_save, info_site, i);
    %  export_fig(filepath_fig, '-m1');
end
%}
%% Figure to check false alarms/ misses
%{
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
i = 13;
j = 16;
P_LCSt =  (match12_R{i,j}*P_LCS{i})+ repmat(match12_t{i,j},1,n_tree(i));
P_dist = sqrt(sum((P_LCSt - P_LCS{i}).^2));
P_color = (double(vec2cmap(P_dist, 'jet' )))./255;
figure;
title(sprintf('RMSE = %3.2f', match12_RMSE(i,j)));
hold on;
axis equal;
x_axist = match12_R{i,j}*x_axis + repmat(match12_t{i,j},[1,2]);
y_axist = match12_R{i,j}*y_axis + repmat(match12_t{i,j},[1,2]);
z_axist = match12_R{i,j}*z_axis + repmat(match12_t{i,j},[1,2]);
plot3(x_axis(1,:),x_axis(2,:),x_axis(3,:),'-r', 'linewidth',2)
plot3(y_axis(1,:),y_axis(2,:),y_axis(3,:),'-g', 'linewidth',2)
plot3(z_axis(1,:),z_axis(2,:),z_axis(3,:),'-b', 'linewidth',2)
plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
for p= 1:n_tree(i);
    h =  filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(p,:));
    set(h, 'facealpha', 0.5);
end
for p= 1:n_tree(j);
    h = filledCircle([P_LCSt(1,p); P_LCSt(2,p)]',P_rad{i}(p),1000,P_color(p,:));
    set(h, 'facealpha', 0.5);
end
% Show in WCS
P_LCSt =  (match_R{i,j}*P_LCS{j})+ repmat(match_t{i,j},1,n_tree(j));
% P_dist = sqrt(sum((P_LCSt - P_LCS{i}).^2));
% P_color = (double(vec2cmap(P_dist, 'jet' )))./255;
figure;
title(sprintf('RMSE = %3.2f', match12_RMSE(i,j)));
hold on;
axis equal;
x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
x_axis1t = match1_R{i,j}*x_axis + repmat(match1_t{i,j},[1,2]);
y_axis1t = match1_R{i,j}*y_axis + repmat(match1_t{i,j},[1,2]);
z_axis1t = match1_R{i,j}*z_axis + repmat(match1_t{i,j},[1,2]);
plot3(x_axis(1,:),x_axis(2,:),x_axis(3,:),'-r', 'linewidth',2)
plot3(y_axis(1,:),y_axis(2,:),y_axis(3,:),'-g', 'linewidth',2)
plot3(z_axis(1,:),z_axis(2,:),z_axis(3,:),'-b', 'linewidth',2)
plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
plot3(x_axis1t(1,:),x_axis1t(2,:),x_axis1t(3,:),'-r', 'linewidth',2)
plot3(y_axis1t(1,:),y_axis1t(2,:),y_axis1t(3,:),'-g', 'linewidth',2)
plot3(z_axis1t(1,:),z_axis1t(2,:),z_axis1t(3,:),'-b', 'linewidth',2)
for p= 1:n_tree(i);
    h =  filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,'r');
    set(h, 'facealpha', 0.5);
end
for p= 1:n_tree(j);
    h = filledCircle([P_LCSt(1,p); P_LCSt(2,p)]',P_rad{i}(p),1000,'b');
    set(h, 'facealpha', 0.5);
end
% Text on axes
textloc = (x_axis(:,2) + y_axis(:,2))/2;
textstr = sprintf('%g', i);
text(textloc(1), textloc(2), textloc(3), textstr);
textloc = (x_axist(:,2) + y_axist(:,2))/2;
textstr = sprintf('%g', j);
text(textloc(1), textloc(2), textloc(3), textstr);
textloc = (x_axis1t(:,2) + y_axis1t(:,2))/2;
textstr = sprintf('%g OW', j);
text(textloc(1), textloc(2), textloc(3), textstr);
% Plot scanner
h =  filledCircle([0,0],10,1000,'r');
set(h, 'facealpha', 0);
h =  filledCircle(match_t{i,j}(1:2),10,1000,'b');
set(h, 'facealpha', 0);
axis equal


clear i j P_dist P_LCSt
% Outputs
% match12_R               {i,j} Effective rotation
% match12_t               {i,j} Effective translation
%   save('line1309.mat');
%else
%    load('line1309.mat');
%end
%}
%% Plot pose for each i,j
%{
%
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
[x_sph, y_sph, z_sph] = sphere;

for i= 13;%:n_S;
    err_alpha = match12_RMSE(i,:);
    err_color = (double(vec2cmap(err_alpha, 'jet', 0,2 )))./255;
    figure;
    title(sprintf('Pairwise Pose Estimates for Site %3.0f-%g', info_site,i));
    hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    for j = 1:n_S;
        x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
        y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
        z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
        
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
        h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
            z_sph+match_t{i,j}(3));
        alpha(0.2)
        set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
    end
    set(gcf, 'color', 'white')
    filepath_fig = sprintf('%spose_c2_%03.0f-%02.0f.png',...
        path_save, info_site, i);
    %  export_fig(filepath_fig, '-m1');
end
%}
%% Determine Pose Consistency in WCS
% Load variables
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
match_R = match1_R;
match_t = match1_t;
foo = 1; 
   %% Effective (WCS) R and t using Dijkstra's
    %
% Shortest path to WCS
% Weighted undirected adjacency matrix
G_path = cell(1,n_S);
for i = 1:n_S
for j = 1:n_S;
    [~,G_path{i,j},~] = graphshortestpath(sparse(match12_RMSE),j,i);
end
end

% Effective R and t back to WCS
G_R = cell(1,n_S);
G_t = cell(1,n_S);
for i = 1:n_S;
for j = 1:n_S;
    path = G_path{i,j};
    Rtemp = eye(3);
    ttemp = zeros(3,1);
    for k = 1:numel(path)-1;
        if any(size(match_R{i,k})~=size(Rtemp));
            break
        end
        Rtemp = match_R{path(k+1),path(k)}*Rtemp;
        ttemp = match_R{path(k+1),path(k)}*(ttemp)+ match_t{path(k+1),path(k)};
    end
    G_R{i,j} = Rtemp;
    G_t{i,j} = ttemp;
end
end

%{
    % Declare match from i-i identity
    for i = 1:n_S;
        match_R{i,i} = eye(3);
        match_t{i,i} = zeros(3,1);
        match_i{i,i} = 1:numel(P_LCS{i});
        match_j{i,i} = 1:numel(P_LCS{i});
    end
%}

 for i = 1:n_S;
        filepath_fig = sprintf('%sgraph_%03.0f-%02.0f.png',...
            path_save, info_site, i);
    if exist(filepath_fig, 'file');
        continue
    end
        x_axis = [0 0 0; 1 0 0]';
        y_axis = [0 0 0; 0 1 0]';
        z_axis = [0 0 0; 0 0 1]';
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
       % err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        figure;
        title(sprintf('Graph Locations for Site %3.0f-%g', info_site,i));
        hold on;
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        axis equal
        for j= 1:n_S;
            hold on
            x_axist = G_R{i,j}*x_axis + repmat(G_t{i,j},[1,2]);
            y_axist = G_R{i,j}*y_axis + repmat(G_t{i,j},[1,2]);
            z_axist = G_R{i,j}*z_axis + repmat(G_t{i,j},[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            %h = surf(x_sph+G_t{i,j}(1), y_sph+G_t{i,j}(2), ...
            %    z_sph+G_t{i,j}(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
        set(gcf, 'color', 'white')
        filepath_fig = sprintf('%sgraph_%03.0f-%02.0f.png',...
            path_save, info_site, i);
        export_fig(filepath_fig, '-m1');
        close(gcf);
 end
    
        G_rx =  180*cellfun(@decompose_rotation_rx,G_R)/pi;
        G_ry =  180*cellfun(@decompose_rotation_ry,G_R)/pi;
        G_rz =  180*cellfun(@decompose_rotation_rz,G_R)/pi;
        G_tx =  cellfun(@decompose_translation_tx,G_t);
        G_ty =  cellfun(@decompose_translation_ty,G_t);
        G_tz =  cellfun(@decompose_translation_tz,G_t); 
        
        filepath_G_R = sprintf('%s%s',path_mat, 'G_R.mat');
        if ~exist(filepath_G_R, 'file');
        save(filepath_G_R, 'G_R');
        filepath_G_t = sprintf('%s%s',path_mat, 'G_t.mat');
        save(filepath_G_t, 'G_t');
        filepath_G_tx = sprintf('%s%s',path_mat, 'G_tx.mat');
        save(filepath_G_tx, 'G_tx');        
        filepath_G_ty = sprintf('%s%s',path_mat, 'G_ty.mat');
        save(filepath_G_ty, 'G_ty');        
        filepath_G_tz = sprintf('%s%s',path_mat, 'G_tz.mat');
        save(filepath_G_tz, 'G_tz');
        filepath_G_rx = sprintf('%s%s',path_mat, 'G_rx.mat');
        save(filepath_G_rx, 'G_rx');        
        filepath_G_ry = sprintf('%s%s',path_mat, 'G_ry.mat');
        save(filepath_G_ry, 'G_ry');        
        filepath_G_rz = sprintf('%s%s',path_mat, 'G_rz.mat');
        save(filepath_G_rz, 'G_rz');  
        end
        clear i j path Rtemp ttemp k Gsparse
% Outputs
% G_path               (j) best path to i
% G_R                  (j) effective Rotation to i
% G_t                  (j) effective translation to i
% G_rx 
% G_tx 
    %}
%% Create explicit S_T_LCS cell array
%
% Define S_P_LCS array
%             save(sprintf('%sline1238.mat', path_matvar));
%       else
%           load(sprintf('%sline1238.mat',path_matvar));
%       end

% Redefine to use identical code 
match_R = G_R;
match_t = G_t; 

match_rx = 180*cellfun(@decompose_rotation_rx,match_R)/pi;
match_ry = 180*cellfun(@decompose_rotation_ry,match_R)/pi;
match_rz = 180*cellfun(@decompose_rotation_rz,match_R)/pi;

S_P_LCS = cell(n_S,1);
S_r_LCS = cell(n_S,1);
for i = 1:n_S;
    S_P_LCS{i} = zeros(3,n_S);
    S_r_LCS{i} = zeros(3,n_S);
    for j = 1:n_S;
        if isempty(match_t{i,j});
            continue
        end
        S_P_LCS{i}(:,j) = match_t{i,j}';
        S_r_LCS{i}(:,j) = [match_rx(i,j); match_ry(i,j); match_rz(i,j)]; % not used
    end
end
S_R_LCS = match_R;

% Plot S_P_LC
%{
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
[x_sph, y_sph, z_sph] = sphere;
for i= 1;%:n_S;
    clear alpha
    err_alpha = match12_RMSE(:,i);
    err_color = (double(vec2cmap(err_alpha, 'jet', 0,2 )))./255;
    figure;
    title(sprintf('Pairwise Pose Estimates for Site %3.0f-%g', info_site,i));
    hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    for j = 1:n_S;
        x_axist = match_R{i,j}*x_axis + repmat(S_P_LCS{i}(:,j),[1,2]);
        y_axist = match_R{i,j}*y_axis + repmat(S_P_LCS{i}(:,j),[1,2]);
        z_axist = match_R{i,j}*z_axis + repmat(S_P_LCS{i}(:,j),[1,2]);
        
        plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
        h = surf(x_sph+S_P_LCS{i}(1,j), y_sph+S_P_LCS{i}(2,j), ...
            z_sph+S_P_LCS{i}(3,j));
        alpha(0.2)
        set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
    end
    set(gcf, 'color', 'white')
    filepath_fig = sprintf('%spose_c2_%03.0f-%02.0f.png',...
        path_save, info_site, i);
    % export_fig(filepath_fig, '-m1');
end
%}
%
% Outputs:
%       S_P_LCS{i}                - locations of sensors in LCSi
%       S_R_LCS                   - rotations of sensors in LCSi
%}
%% Find connection network for scan poses
%
% All possible triangles combinations for scan locations

% * Note that S_T_comb{s} is identical
if n_S < 255;
    FDT_S_T_comb = @(x) uint8(x);
elseif n_S < 65535;
    FDT_S_T_comb = @(x) uint16(x);
end
if n_S==0;
    return
end
S_T_comb = uint8(combnk(1:n_S,3));
S_n_comb = size(S_T_comb,1);

%{
s = 1;
hold on
for c = 1:200;% S_n_comb;
    plot([S_P_LCS{s}(1,S_T_comb{s}(c,1)) S_P_LCS{s}(1,T_comb{s}(c,2))],...
         [S_P_LCS{s}(2,S_T_comb{s}(c,1)) S_P_LCS{s}(2,T_comb{s}(c,2))],'-k')
    plot([S_P_LCS{s}(1,S_T_comb{s}(c,1)) S_P_LCS{s}(1,T_comb{s}(c,3))],...
         [S_P_LCS{s}(2,S_T_comb{s}(c,1)) S_P_LCS{s}(2,T_comb{s}(c,3))],'-k')
    plot([S_P_LCS{s}(1,S_T_comb{s}(c,2)) S_P_LCS{s}(1,T_comb{s}(c,3))],...
         [S_P_LCS{s}(2,S_T_comb{s}(c,2)) S_P_LCS{s}(2,T_comb{s}(c,3))],'-k')
end
%}

% Outputs:
%       S_T_comb                  - Combinations of scans; n_P_comb{i} x 3
%       S_n_comb                  - number of combinations of scans
%       FDT_S_T_comb              - FDT = Datatype function
foo = 1;
%}
%% Make explicit arrays of triangles for scans
%
S_T_LCS = cell(n_S,1);
% Note these are not sorted as before... now need permutations
for s = 1:n_S;
    S_T_LCS{s} = zeros(S_n_comb,3,3);
    S_T_LCS{s}(:,:,1) = reshape(S_P_LCS{s}(1,S_T_comb),[S_n_comb,3]); % x values
    S_T_LCS{s}(:,:,2) = reshape(S_P_LCS{s}(2,S_T_comb),[S_n_comb,3]); % y values
    S_T_LCS{s}(:,:,3) = reshape(S_P_LCS{s}(3,S_T_comb),[S_n_comb,3]); % z values
end

clear s
% Outputs:
%       S_T_LCS                   - Triplets of scan pts for each comb
foo = 1;
%}
%% Find eigenvalues of each triangle
%
S_T_eig = cell(n_S,1);
for s = 1:n_S;
    S_T_eig{s} = zeros(S_n_comb,2);
    for t = 1:S_n_comb;
        temp = eig(cov(squeeze(S_T_LCS{s}(t,:,:))));
        S_T_eig{s}(t,:) = temp(2:3);
    end
end

S_T_normeig = cell(n_S,1);
S_T_isncoll = cell(n_S,1);
for s = 1:n_S;
    S_T_normeig{s} = S_T_eig{s}./repmat(sum(S_T_eig{s},2),[1,2]);
    S_T_isncoll{s} = S_T_normeig{s}(:,1)> t_coll;
end

% Figure
%{
%}

clear s t temp edges legend_str cov_data white_data percent_var n_var var_color
clear C L i
% Outputs:
%       S_T_eig(i)               - Triangle eigenvalues
%       S_T_normeig(i)           - Normalized triangle eigenvalues
%       S_T_isncoll(i)           - List of triangles which are not collinear
%}
%% Create 1D array of valid triangles after filtering by collinearity
%
S_m_combix = cell(n_S,1); % index to valid combinations
S_m_eig = cell(n_S,1); % 1d array of eigenvalues for valid combinations
S_n_m = zeros(n_S,1); % number of valid combinations
for i = 1:n_S;
    S_m_combix{i} = 1:S_n_comb;
    S_m_combix{i} = S_m_combix{i}(S_T_isncoll{i}); % Filter by collinearity
    S_m_eig{i} = S_T_eig{i}(S_m_combix{i},:); % Lowercase m denotes filtering by collinearity
    S_n_m(i) = numel(S_m_combix{i});
end

clear i
% Outputs:
%       m_combix                - [cleared next block] index of combinations which are not collinear
%       m_eig                   - [cleared next block] Eigenvalues
%       n_m                     - [cleared next block] number of valid combinations
%}
%% Find likely RANSAC pairs by filtering eigenvalues
%
S_M_bestixI = cell(n_S,n_S);
S_M_bestixJ = cell(n_S,n_S);
S_M_combI = cell(n_S,n_S);
S_M_combJ = cell(n_S,n_S);
S_M_eigerror = cell(n_S,n_S);
S_n_M = zeros(n_S,n_S);
for i = 1:n_S;
    for j = i+1:n_S; % Symmetric
        % Replicate eigenvalues to determine error
        M_eigaI = repmat(S_m_eig{i}(:,1), [1,S_n_m(j)]);
        M_eigbI = repmat(S_m_eig{i}(:,2), [1,S_n_m(j)]);
        M_eigaJ = repmat(S_m_eig{j}(:,1)', [S_n_m(i),1]);
        M_eigbJ = repmat(S_m_eig{j}(:,2)', [S_n_m(i),1]);
        M_error = (M_eigaI - M_eigaJ).^2 + (M_eigbI - M_eigbJ).^2;
        % Brute force removal
        is_brute = (M_error < t_eig_error);
        %xx denotes 1d array of valid
        xxM_error = M_error(is_brute);
        % Work in 1d instead of 2d matrix
        [row, col] = find(is_brute);
        xx_combixI = S_m_combix{i}(row)'; %1d array of all valid comb indices
        xx_combixJ = S_m_combix{j}(col)';
        xx_combI = S_T_comb(xx_combixI,:);
        xx_combJ = S_T_comb(xx_combixJ,:);
        % Now sort remaining
        % Sort eigenvalue error
        [xxM_error_sort,sortix] = sort(xxM_error);
        % Update indices of best triangles
        S_M_eigerror{i,j} = xxM_error_sort;
        S_M_bestixI{i,j} = xx_combixI(sortix);
        S_M_bestixJ{i,j} = xx_combixJ(sortix);
        S_M_combI{i,j} = xx_combI(sortix,:);
        S_M_combJ{i,j} = xx_combJ(sortix,:);
        S_n_M(i,j) = numel(sortix);
    end
end
clear match_israd match_rad_diff
clear m_eig m_combix n_m
clear M_eigaI M_eigaJ M_eigbI M_eigbJ M_error is_brute
clear row col xx_combixI xx_combixJ xx_combI xx_combJ n_brute
clear xx_israd linear xxM_error xxM_error_sort sortix r i j
% S_M_bestixI{i,j}          - index of best i triangles for j-> i
% S_M_bestixJ{i,j}          - index of best j triangles for j-> i
% S_M_combI{i,j}            - triangle combinations
% S_M_combJ{i,j}            - triangle combinations
% S_M_eigerror{i,j}         - Eigenvalue error for each triangle match
%}
%% Construct ordered tri radius, locations for clarity for scans
%
S_M_triI = cell(n_S,n_S);
S_M_triJ = cell(n_S,n_S);
S_M_triradI = cell(n_S,n_S);
S_M_triradJ = cell(n_S,n_S);
for i = 1:n_S;
    for j = i + 1:n_S;
        S_M_triI{i,j} = zeros(S_n_M(i,j),3,3);
        S_M_triI{i,j}(:,1,:) = S_P_LCS{i}(:,S_M_combI{i,j}(:,1))';
        S_M_triI{i,j}(:,2,:) = S_P_LCS{i}(:,S_M_combI{i,j}(:,2))';
        S_M_triI{i,j}(:,3,:) = S_P_LCS{i}(:,S_M_combI{i,j}(:,3))';
        S_M_triJ{i,j} = zeros(S_n_M(i,j),3,3);
        S_M_triJ{i,j}(:,1,:) = S_P_LCS{j}(:,S_M_combJ{i,j}(:,1))';
        S_M_triJ{i,j}(:,2,:) = S_P_LCS{j}(:,S_M_combJ{i,j}(:,2))';
        S_M_triJ{i,j}(:,3,:) = S_P_LCS{j}(:,S_M_combJ{i,j}(:,3))';
    end
end
%{
i = 1;
j = 2;
M_raderror = abs((M_triradI{i,j}-M_triradJ{i,j})./((M_triradI{i,j}+M_triradJ{i,j})/2));
M_eigerror_check = zeros(n_M(i,j),1);
for t = 1:n_M(i,j);
    eI = eig(cov(squeeze(M_triI{i,j}(t,:,:))));
    eJ = eig(cov(squeeze(M_triJ{i,j}(t,:,:))));
    M_eigerror_check(t) = sum((eI-eJ).^2);
%{
    triI = squeeze(M_triI{i,j}(t,:,:));
    triJ = squeeze(M_triJ{i,j}(t,:,:));
    figure;
    plot_triangle(triI);
    plot_triangle(triJ);
    hold on; axis equal;
%}
end
%}

clear i j
% Outputs:
% S_M_triI{i,j}             - triangle of points for set I
% S_M_triJ{i,j}             - triangle of points for set J
% S_M_triradI{i,j}          - radius of triangle of points for set I
% S_M_triradJ{i,j}          - radius of triangle of points for set J
%}
%% Find pairwise matches between full scans
%
fprintf('\nFind pairwise matches \n');

% Initialize graph and corresponding functional relationships
S_match_i = cell(n_S); %base
S_match_j = cell(n_S); %mobile
S_match_R = cell(n_S);
S_match_t = cell(n_S);
S_match_nit = nan(n_S);
% Find rotation and tranlation from j to i along with point matching pairs
for i = 1:n_S-1;
    for j = i+1:n_S; % Symmetric
        fprintf('\n\tMatching Scan %g to %g\n',j,i);
        
        % RANSAC setup
        n_match_best = 0;
        m1_best = [];
        m2_best = [];
        n_search = min(t_RANSAC_nsearch, S_n_M(i,j));
        
        for r = 1:n_search;
            
            triI = squeeze(S_M_triI{i,j}(r,:,:));
            triJ = squeeze(S_M_triJ{i,j}(r,:,:));
            
            % Check input Plotting
            %{
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
                end
                for p= 1:n_tree(j);
                    filledCircle([P_LCS{j}(1,p); P_LCS{j}(2,p)]',P_rad{j}(p),1000,P_color(j,:));
                end
                plot_triangle(triI);
                plot_triangle(triJ);
                plot_triangle_connection(triI,triJ);
            %}
            
            [Rhat,that] = rigid_transform_3D(triJ,triI);
            dataJt = ((Rhat*S_P_LCS{j})+ repmat(that,1,n_S)); %Note change
            % Add rotation information
            % Use n_match, round-trip or other method to determine fit
            P_LCSIrep = repmat(S_P_LCS{i}',[1,1,n_S]);
            P_LCSJtrep = zeros(n_S,3,n_S);
            for t = 1:n_S;
                P_LCSJtrep(t,:,:) = dataJt;
            end
            error_xyz = squeeze(sum((P_LCSIrep-P_LCSJtrep).^2,2));
            [m1,m2] = find(error_xyz<t_RANSAC_xyz);
            n_match = numel(m1);
            
            % Check error calculation
            %{
                triJt = ((Rhat*triJ')+ repmat(that,1,3))';
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
                end
                for p= 1:n_tree(j);
                    filledCircle([dataJt(1,p); dataJt(2,p)]',P_rad{j}(p),1000,P_color(j,:));
                end
                dist_color = plyintensity2color( error_xyz(1,:), 'jet' )./255;
                for c = 1:n_tree(j);
                    plot([P_LCS{i}(1,1), dataJt(1,c)],...
                        [P_LCS{i}(2,1), dataJt(2,c)],...
                        'color', dist_color(c,:), 'linewidth', 2);
                end
            %}
            
            if n_match > n_match_best;
                n_match_best = n_match;
                R_best = Rhat;
                t_best = that;
                m1_best = m1;
                m2_best = m2;
                nit_best = r;
            end
            
            % Check output Plotting
            %{
                figure; hold on; axis equal;
                for p= 1:n_tree(i);
                    filledCircle([P_LCS{i}(1,p); P_LCS{i}(2,p)]',P_rad{i}(p),1000,P_color(i,:));
                end
                for p= 1:n_tree(j);
                    filledCircle([dataJt(1,p); dataJt(2,p)]',P_rad{j}(p),1000,P_color(j,:));
                end
                plot_triangle(triI);
                plot_triangle(triJt);
                plot_triangle_connection(triI,triJt);
            %}
        end
        
        if n_match_best > 0;
            S_match_R{i,j} = R_best;
            S_match_t{i,j} = t_best;
            S_match_i{i,j} = m1_best;
            S_match_j{i,j} = m2_best;
            S_match_nit(i,j) = nit_best;
        end
    end
end

% Declare match from i-i identity
for i = 1:n_S;
    S_match_R{i,i} = eye(3);
    S_match_t{i,i} = zeros(3,1);
    S_match_i{i,i} = 1:numel(P_LCS{i}); % not valid
    S_match_j{i,i} = 1:numel(P_LCS{i});
end

% Make R,t non-directed
for i = 1:n_S;
    for j = i:n_S;
        if ~isempty(match_R{i,j});
            S_match_R{j,i} = match_R{i,j}';
            S_match_t{j,i} = -(match_R{i,j}')*match_t{i,j};
        end
    end
end

clear n_match_best m1_best m2_best n_search i j r Rhat that dataJt
clear P_LCSIrep P_LCSJtrep t error_xyz m1 m2 n_match
% Outputs
% match_R               {i,j} is pairwise rotation from j into i
% match_t               {i,j} is pairwise translation from j into i
% match_i               {i,j} is pairwise matches from i
% match_j               {i,j} is pairwise matches from j
%            save(sprintf('%sline1238.mat', path_matvar));
%       else
%           load(sprintf('%sline1238.mat',path_matvar));
%       end
%}
%%
%
S_t = cell(n_S,n_S);
S_R = cell(n_S,n_S);
for i = 1:n_S; % Choose reference CS
    for j = 1:n_S;
        if isempty(S_match_R{i,j});
            continue
        end
        S_t{i,j} = S_match_R{i,j}*S_P_LCS{j} + repmat(S_match_t{i,j}, 1,n_S);
        S_R{i,j} = cell(n_S,1);
        for s = 1:n_S;
            if isempty(S_R_LCS{j,s});
                continue
            end
            S_R{i,j}{s} = S_match_R{i,j}*S_R_LCS{j,s}; % inverse for j<s??
        end
    end
end

% Plot S_P_LC
%
x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
[x_sph, y_sph, z_sph] = sphere;
for i = 3:n_S;
    filepath_fig = sprintf('%salldijk_%03.0f-%02.0f.png',...
        path_save, info_site, i);
    if exist(filepath_fig, 'file');
   %     continue
    end
    clear alpha
    figure;
    title(sprintf('Pairwise Pose Estimates for Site %3.0f-%g', info_site,i));
    hold on;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    axis equal
    for j= 1:n_S;
        %  err_alpha = match12_RMSE(i,:); % should change to use global cmap
        %  err_color = (double(vec2cmap(err_alpha, 'jet', 0,2 )))./255;
        for s = 1:n_S;
            if isempty(S_R{i,j});
                continue
            end
            if isempty(S_R{i,j}{s});
                continue
            end
            x_axist = S_R{i,j}{s}*x_axis + repmat(S_t{i,j}(:,s),[1,2]);
            y_axist = S_R{i,j}{s}*y_axis + repmat(S_t{i,j}(:,s),[1,2]);
            z_axist = S_R{i,j}{s}*z_axis + repmat(S_t{i,j}(:,s),[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', s);
            text(textloc(1), textloc(2), textloc(3), textstr);
            foo = 1;
            %    h = surf(x_sph+S_t{i,j}(1,s), y_sph+S_t{i,j}(2,s), ...
            %        z_sph+S_t{i,j}(3,s));
            %     alpha(0.2)
            %     set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
    end
    set(gcf, 'color', 'white')
    export_fig(filepath_fig, '-m1');
    delete(gcf);
end
%          save(sprintf('%sline1949.mat', path_matvar));
%     else
%         load(sprintf('%sline1949.mat',path_matvar));
%     end%% Find true location based on
%%
filepath_S_match_R = sprintf('%s%s',path_mat, 'S_match_R.mat');
save(filepath_S_match_R, 'S_match_R');
filepath_S_match_t = sprintf('%s%s',path_mat, 'S_match_t.mat');
save(filepath_S_match_t, 'S_match_t');
filepath_S_R = sprintf('%s%s',path_mat, 'S_R.mat');
save(filepath_S_R, 'S_R');
filepath_S_t = sprintf('%s%s',path_mat, 'S_t.mat');
save(filepath_S_t, 'S_t');
% else
%     filepath_S_match_R = sprintf('%s%s',path_mat, 'S_match_R.mat');
%     load(filepath_S_match_R);
%     filepath_S_match_t = sprintf('%s%s',path_mat, 'S_match_t.mat');
%     load(filepath_S_match_t);
%     filepath_S_R = sprintf('%s%s',path_mat, 'S_R.mat');
%     load(filepath_S_R);
%     filepath_S_t = sprintf('%s%s',path_mat, 'S_t.mat');
%     load(filepath_S_t);
% end


S_rx = cell(n_S,n_S);
S_ry = cell(n_S,n_S);
S_rz = cell(n_S,n_S);
S_tx = cell(n_S,n_S);
S_ty = cell(n_S,n_S);
S_tz = cell(n_S,n_S);
for i = 1:n_S;
    for j  = 1:n_S;
        if isempty(S_R{i,j});
            S_rx{i,j} = nan;
            S_ry{i,j} = nan;
            S_rz{i,j} = nan;
            S_tx{i,j} = nan;
            S_ty{i,j} = nan;
            S_tz{i,j} = nan;
            continue
        end
        S_rx{i,j} =  180*cellfun(@decompose_rotation_rx,S_R{i,j})/pi;
        S_ry{i,j} =  180*cellfun(@decompose_rotation_ry,S_R{i,j})/pi;
        S_rz{i,j} =  180*cellfun(@decompose_rotation_rz,S_R{i,j})/pi;
        S_tx{i,j} =  S_t{i,j}(1,:)';
        S_ty{i,j} =  S_t{i,j}(2,:)';
        S_tz{i,j} =  S_t{i,j}(3,:)';
    end
end

Sall_rx = zeros(n_S,n_S,n_S);
Sall_ry = zeros(n_S,n_S,n_S);
Sall_rz = zeros(n_S,n_S,n_S);
Sall_tx = zeros(n_S,n_S,n_S);
Sall_ty = zeros(n_S,n_S,n_S);
Sall_tz = zeros(n_S,n_S,n_S);

for i = 1:n_S;
    for j = 1:n_S;
        Sall_rx(i,j,:) = S_rx{i,j};
        Sall_ry(i,j,:) = S_ry{i,j};
        Sall_rz(i,j,:) = S_rz{i,j};
        Sall_tx(i,j,:) = S_tx{i,j};
        Sall_ty(i,j,:) = S_ty{i,j};
        Sall_tz(i,j,:) = S_tz{i,j};
    end
end

% Median - doesn't work
%{
i = 13;
s = 1;
clear alpha
figure;
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
scatter3(Sall_tx(i,:,s), Sall_ty(i,:,s), Sall_tz(i,:,s));
for j = 1:n_S;
        textstr = sprintf('%g', s);
        text(Sall_tx(i,j,s), Sall_ty(i,j,s), Sall_tz(i,j,s), textstr);
end
set(gcf, 'color', 'white')
plot3(nanmedian(Sall_tx(i,:,s)), nanmedian(Sall_ty(i,:,s)), nanmedian(Sall_tz(i,:,s)), 'xr', 'markersize',20,'linewidth',2);
    
med_rx = squeeze(nanmedian(Sall_rx,2));
med_ry = squeeze(nanmedian(Sall_ry,2));
med_rz = squeeze(nanmedian(Sall_rz,2));
med_tx = squeeze(nanmedian(Sall_tx,2));
med_ty = squeeze(nanmedian(Sall_ty,2));
med_tz = squeeze(nanmedian(Sall_tz,2));
med_is = cell(n_S,n_S);
med_conf = zeros(n_S,n_S);
t_medrad = 5;
t_meddist = .5;
for i =1:n_S;
    for j = 1:n_S;
        med_is{i,j} = ...
            (Sall_rx(i,:,j)-med_rx(i,j) < t_medrad) &...
            (Sall_ry(i,:,j)-med_ry(i,j) < t_medrad) &...
            (Sall_rz(i,:,j)-med_rz(i,j) < t_medrad) &...
            (Sall_tx(i,:,j)-med_tx(i,j) < t_meddist) &...
            (Sall_ty(i,:,j)-med_ty(i,j) < t_meddist) &...
            (Sall_tz(i,:,j)-med_tz(i,j) < t_meddist);
        med_conf(i,j) = sum(med_is{i,j})./n_S;
    end
end

med_rx_rad = deg2rad(med_rx);
med_ry_rad = deg2rad(med_ry);
med_rz_rad = deg2rad(med_rz);
med_R = cell(n_S,n_S);
med_t = cell(n_S,n_S);
for i = 1:n_S;
    for j = 1:n_S;
        med_R{i,j} = compose_rotation(med_rx_rad(i,j), med_ry_rad(i,j), med_rz_rad(i,j));
        med_t{i,j} = [med_tx(i,j) med_ty(i,j) med_tz(i,j)]';
    end
end
%}

ran_rx = zeros(n_S,n_S);
ran_ry = zeros(n_S,n_S);
ran_rz = zeros(n_S,n_S);
ran_tx = zeros(n_S,n_S);
ran_ty = zeros(n_S,n_S);
ran_tz = zeros(n_S,n_S);
ran_tdist = 1;
ran_trad = 5;
ran_conf = zeros(n_S,n_S);
% RANSAC
for i = 1:n_S;
    for j = 1:n_S;
        rxdiff = abs(repmat(Sall_rx(i,:,j),[n_S,1]) - repmat(Sall_rx(i,:,j)',[1,n_S]));
        rydiff = abs(repmat(Sall_ry(i,:,j),[n_S,1]) - repmat(Sall_ry(i,:,j)',[1,n_S]));
        rzdiff = abs(repmat(Sall_rz(i,:,j),[n_S,1]) - repmat(Sall_rz(i,:,j)',[1,n_S]));
        txdiff = abs(repmat(Sall_tx(i,:,j),[n_S,1]) - repmat(Sall_tx(i,:,j)',[1,n_S]));
        tydiff = abs(repmat(Sall_ty(i,:,j),[n_S,1]) - repmat(Sall_ty(i,:,j)',[1,n_S]));
        tzdiff = abs(repmat(Sall_tz(i,:,j),[n_S,1]) - repmat(Sall_tz(i,:,j)',[1,n_S]));
        tdiff = sqrt(txdiff.^2 + tydiff.^2 + tzdiff.^2);
        ran_inlier = (tdiff  < ran_tdist) & (rxdiff < ran_trad) & ...
            (rydiff < ran_trad) & (rzdiff < ran_trad);
        is_best = (sum(ran_inlier,2) == max(sum(ran_inlier,2)));
        ran_rx(i,j) = mean(Sall_rx(i,is_best,j));
        ran_ry(i,j) = mean(Sall_ry(i,is_best,j));
        ran_rz(i,j) = mean(Sall_rz(i,is_best,j));
        ran_tx(i,j) = mean(Sall_tx(i,is_best,j));
        ran_ty(i,j) = mean(Sall_ty(i,is_best,j));
        ran_tz(i,j) = mean(Sall_tz(i,is_best,j));
        ran_conf(i,j) = sum(is_best)./n_S;
        %{
        figure;
        hold on
        x_axis = [0 0 0; 1 0 0]';
        y_axis = [0 0 0; 0 1 0]';
        z_axis = [0 0 0; 0 0 1]';
        %[x_sph, y_sph, z_sph] = sphere;
        %err_alpha = ran_conf(i,:); % should change to use global cmap
        %err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        %clear alpha
        for s = 1:n_S;
            Sall_R = compose_rotation(deg2rad(Sall_rx(i,s,j)), deg2rad(Sall_ry(i,s,j)), deg2rad(Sall_rz(i,s,j)));
            Sall_t = [Sall_tx(i,s,j) Sall_ty(i,s,j) Sall_tz(i,s,j)]';
            x_axist = Sall_R*x_axis + repmat(Sall_t,[1,2]);
            y_axist = Sall_R*y_axis + repmat(Sall_t,[1,2]);
            z_axist = Sall_R*z_axis + repmat(Sall_t,[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', s);
            text(textloc(1), textloc(2), textloc(3), textstr);
            %h = surf(x_sph+ran_t(1), y_sph+ran_t(2), ...
            %        z_sph+ran_t(3));
            %alpha(0.2)
            %set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
        %}
    end    
    end
    
    ran_rx_rad = deg2rad(ran_rx);
    ran_ry_rad = deg2rad(ran_ry);
    ran_rz_rad = deg2rad(ran_rz);
    ran_R = cell(n_S,n_S);
    ran_t = cell(n_S,n_S);
    for i = 1:n_S;
        for j = 1:n_S;
            ran_R{i,j} = compose_rotation(ran_rx_rad(i,j), ran_ry_rad(i,j), ran_rz_rad(i,j));
            ran_t{i,j} = [ran_tx(i,j) ran_ty(i,j) ran_tz(i,j)]';
        end
    end
    
    filepath_ran_rx = sprintf('%s%s',path_mat, 'ran_rx.mat');
    save(filepath_ran_rx, 'ran_rx');
    filepath_ran_ry = sprintf('%s%s',path_mat, 'ran_ry.mat');
    save(filepath_ran_ry, 'ran_ry');
    filepath_ran_rz = sprintf('%s%s',path_mat, 'ran_rz.mat');
    save(filepath_ran_rz, 'ran_rz');
    filepath_ran_tx = sprintf('%s%s',path_mat, 'ran_tx.mat');
    save(filepath_ran_tx, 'ran_tx');
    filepath_ran_ty = sprintf('%s%s',path_mat, 'ran_ty.mat');
    save(filepath_ran_ty, 'ran_ty');
    filepath_ran_tz = sprintf('%s%s',path_mat, 'ran_tz.mat');
    save(filepath_ran_tz, 'ran_tz');
    filepath_ran_conf = sprintf('%s%s',path_mat, 'ran_conf.mat');
    save(filepath_ran_conf, 'ran_conf');
    filepath_ran_R = sprintf('%s%s',path_mat, 'ran_R.mat');
    save(filepath_ran_R, 'ran_R');
    filepath_ran_t = sprintf('%s%s',path_mat, 'ran_t.mat');
    save(filepath_ran_t, 'ran_t');
    
    % else
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
    % end
    
    % Plot consensus locations
    %
    for i = 1:n_S;
        x_axis = [0 0 0; 1 0 0]';
        y_axis = [0 0 0; 0 1 0]';
        z_axis = [0 0 0; 0 0 1]';
        [x_sph, y_sph, z_sph] = sphere;
        err_alpha = ran_conf(i,:); % should change to use global cmap
        err_color = (double(vec2cmap(err_alpha, 'jet', 0,1 )))./255;
        clear alpha
        figure;
        title(sprintf('Consensus Locations for Site %3.0f-%g', info_site,i));
        hold on;
        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis equal
        for j= 1:n_S;
            hold on
            x_axist = ran_R{i,j}*x_axis + repmat(ran_t{i,j},[1,2]);
            y_axist = ran_R{i,j}*y_axis + repmat(ran_t{i,j},[1,2]);
            z_axist = ran_R{i,j}*z_axis + repmat(ran_t{i,j},[1,2]);
            
            plot3(x_axist(1,:),x_axist(2,:),x_axist(3,:),'-r', 'linewidth',2)
            plot3(y_axist(1,:),y_axist(2,:),y_axist(3,:),'-g', 'linewidth',2)
            plot3(z_axist(1,:),z_axist(2,:),z_axist(3,:),'-b', 'linewidth',2)
            textloc = (x_axist(:,2) + y_axist(:,2))/2;
            textstr = sprintf('%g', j);
            text(textloc(1), textloc(2), textloc(3), textstr);
            h = surf(x_sph+ran_t{i,j}(1), y_sph+ran_t{i,j}(2), ...
                z_sph+ran_t{i,j}(3));
            alpha(0.2)
            set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
        end
        set(gcf, 'color', 'white')
        filepath_fig = sprintf('%sconsensus_c2_%03.0f-%02.0f.png',...
            path_save, info_site, i);
        export_fig(filepath_fig, '-m1');
        close(gcf);
    end
    foo = 1;
    
    %% Remove bad matches: Instead, let Dijstra decide
    %{
t_n_minmatch = 4;
t_rmin = 20;

match_nmatch = cellfun(@numel,match_i);
match_nmatch = match_nmatch + tril(match_nmatch.',-1);
match_rx = 180*cellfun(@decompose_rotation_rx,match_R)/pi;
match_ry = 180*cellfun(@decompose_rotation_ry,match_R)/pi;
match_rz = 180*cellfun(@decompose_rotation_rz,match_R)/pi;

is_valid = (match_nmatch > t_n_minmatch) & ...
    (match_rx < t_rmin) & (match_ry < t_rmin) & (match_rz < t_rmin);
% is valid should be symmetric
is_valid = (is_valid & is_valid');
match_R(~is_valid) = {[]};
match_t(~is_valid) = {[]};
match_i(~is_valid) = {[]};
match_j(~is_valid) = {[]};

% Plot paths
    %{
P_plot_str = cell(n_S,1);
for s = 1:n_S;
    P_plot_str{s} = sprintf('Plot %g',P_plot(s));
end
view(biograph(triu(is_valid,1),P_plot_str,'ShowArrows','off' ));
    %}

% Outputs
% is_valid              adjacency matrix for scan connections
% t_n_minmatch          minimum number of matches required
% t_rmin                minimum radius for a match
% match_nmatch          number of matches for a given link
clear i j match_rx match_ry match_rz
    %}
    %% Effective (WCS) R and t using Dijkstra's
    %{
% Shortest path to WCS
% Weighted undirected adjacency matrix
G_path = cell(1,n_S);
for j = 2:n_S;
    [~,G_path{1,j},~] = graphshortestpath(sparse(match12_RMSE),j,1);
end
G_path{1,1} = 1;

% Effective R and t back to WCS
G_R = cell(1,n_S);
G_t = cell(1,n_S);
i = 1;
for j = i:n_S;
    path = G_path{1,j};
    Rtemp = eye(3);
    ttemp = zeros(3,1);
    for k = 1:numel(path)-1;
        if any(size(match_R{i,k})~=size(Rtemp));
            break
        end
        Rtemp = match_R{path(k+1),path(k)}*Rtemp;
        ttemp = match_R{path(k+1),path(k)}*(ttemp)+ match_t{path(k+1),path(k)};
    end
    G_R{i,j} = Rtemp;
    G_t{i,j} = ttemp;
end
clear i j path Rtemp ttemp k Gsparse
% Outputs
% G_path               (j) best path to i
% G_R                  (j) effective Rotation to i
% G_t                  (j) effective translation to i
    %}
    %% Transform tree points to WCS
    %{
P_WCS = cell(n_S,1);
for j = 1:n_S;
    if ~isempty(G_R{1,j});
        P_WCS{j} = G_R{1,j}*P_LCS{j} + repmat(G_t{1,j},[1,n_tree(j)]);
    end
end

% Visualization
    %{
figure;
hold on
%i = 1; % Match_Pi is the points which match after effective RT to WCS
for j = 1:n_S;
        for p= 1:n_tree(j);
            filledCircle([P_WCS{j}(1,p); P_WCS{j}(2,p)]',P_rad{j}(p),1000,P_color(j,:));
        end
end
axis equal
title('Points in WCS - Initial Transformation');
grid on
    %}

clear j
% Outputs
% P_WCS                 (j) P_LCS transformed to WCS
    %}
    %% Write transformed PLY
    %{
filepath_plyt = cell(n_S,1);
for j = 1:n_S;
    fprintf('\nWriting ply %g of %g\n',j,n_S);
    filepath_plyt{j} = sprintf('%st.ply',filepath_ply{j}(1:end-4));
    [vertex, ~] = read_ply(filepath_ply{j});
    data_x = vertex(:,1);
    data_y = vertex(:,2);
    data_z = vertex(:,3);
    data2_xyz = (G_R{j}*[data_x data_y data_z]') + ...
        repmat(G_t{j},1,numel(data_x));
    write2ply(filepath_plyt{j},data2_xyz', vertex(:,4:6));
end
%
clear j vertex data_x data_y data_z data2_xyz
% Outputs
% PLY written to disk

%       save('line1392.mat');
%else
%    load('line1392.mat');
%end
    %}
    %% Effective R and t back to WCS
    %{
% Return path could be 1-2-3-4-1 if better than 1-2-3-4-3-2-1
% Way there
G2_R = cell(1,n_S);
G2_t = cell(1,n_S);
i = 13;
for j = 1:n_S;
    path = fliplr(G_path{1,j});
    Rtemp = eye(3);
    ttemp = zeros(3,1);
    for k = 1:numel(path)-1;
        if any(size(match2_R{i,k})~=size(Rtemp));
            break
        end
        Rtemp = match2_R{path(k+1),path(k)}*Rtemp;
        ttemp = match2_R{path(k+1),path(k)}*(ttemp)+ match2_t{path(k+1),path(k)};
    end
    G2_R{i,j} = Rtemp;
    G2_t{i,j} = ttemp;
end
% Way back
G1_R = cell(1,n_S);
G1_t = cell(1,n_S);
for j = 1:n_S;
    path = G_path{1,j};
    Rtemp = eye(3);
    ttemp = zeros(3,1);
    for k = 1:numel(path)-1;
        if any(size(match1_R{i,k})~=size(Rtemp));
            break
        end
        Rtemp = match1_R{path(k),path(k+1)}*Rtemp;
        ttemp = match1_R{path(k),path(k+1)}*(ttemp)+ match1_t{path(k),path(k+1)};
    end
    G1_R{i,j} = Rtemp;
    G1_t{i,j} = ttemp;
end

G12_R = cellfun(@(x,y) x*y, G1_R,G2_R, 'uniformOutput', false);
G12_T = cellfun(@(x,y,z) x*y + z, G1_R,G2_t,G1_t, 'UniformOutput', false);

clear G2_R G2_t path Rtemp ttemp k j i G1_R G1_t
% Outputs
% G12_R                 Circular rotation along disjoint path
% G12_T                 Circular translation along disjoint path
%       save('line1440.mat');
%else
%    load('line1440.mat');
%end
    %}
    %% Confidence Outputs
    %{
% Tranformation Parameters
G12_rx = 180*cellfun(@decompose_rotation_rx,G12_R)/pi;
G12_ry = 180*cellfun(@decompose_rotation_ry,G12_R)/pi;
G12_rz = 180*cellfun(@decompose_rotation_rz,G12_R)/pi;
G12_tx = cellfun(@(x) x(1), G12_T);
G12_ty = cellfun(@(x) x(2), G12_T);
G12_tz = cellfun(@(x) x(3), G12_T);
G12_r = [G12_rx; G12_ry; G12_rz];
G12_t = [G12_tx; G12_ty; G12_tz];
clear G12_rx G12_ry G12_rz G12_tx G12_ty G12_tz
% Outputs
% G12_r                 Circular rotation along disjoint path (arr)
% G12_t                 Circular translation along disjoint path (arr)

% RMSE of tree locations
G12_P_RMSE = zeros(n_S,1);
for j = 1:n_S;
    P_LCSt =  G12_R{j}*P_LCS{j}+ repmat(G12_T{j},1,n_tree(j));
    P_dist = sqrt(sum((P_LCSt - P_LCS{j}).^2));
    G12_P_RMSE(j) = mean(P_dist);
end
clear j P_dist P_LCSt
% Outputs
% G12_P_RMSE            RMSE of tree locations
%

% Sensor Offset
G12_S_dist = sqrt(sum(G12_t.^2,1));
% Outputs
% G12_S_dist            Sensor offset
    %}
    %% Point Cloud Vector Field
    %{
% Initialize volume
interval = 1;
vec_x = transpose(floor(options_axesval(1)):interval:ceil(options_axesval(2)));
vec_y = transpose(floor(options_axesval(3)):interval:ceil(options_axesval(4)));
vec_z = transpose(-5:interval:options_axesval(2));
[vec_X, vec_Y, vec_Z] = meshgrid(vec_x, vec_y, vec_z);
vec_xx = vec_X(:);
vec_yy = vec_Y(:);
vec_zz = vec_Z(:);
vec_xxyyzz = [vec_xx vec_yy vec_zz]';
n_vec = numel(vec_xx);

% Compute error vector
vec_uuvvww = cell(n_S,1);
for j = 2:n_S;
    vec_t =  G12_R{j}*vec_xxyyzz+ repmat(G12_T{j},1,n_vec);
    vec_uuvvww{j} = vec_t - vec_xxyyzz;
end
vec_mag = cellfun(@(x) sqrt(sum(x.^2,1)), vec_uuvvww, 'uniformOutput', false);

% Plot
    %{
j = 2;
max_mag = max(vec_mag{j});
n_quivercolor = 25;
quiverstep = linspace(0,max_mag,n_quivercolor);
quivercolor = jet(n_quivercolor-1);
figure;
hold on
legend_str = cell(n_quivercolor-1,1);
for c = 1:n_quivercolor-1;
    is = vec_mag{j}>=quiverstep(c) & vec_mag{j}<=quiverstep(c+1);
    quiver3(vec_xx(is), vec_yy(is), vec_zz(is),...
        vec_uuvvww{j}(1,is)', vec_uuvvww{j}(2,is)', vec_uuvvww{j}(3,is)',...
        'color', quivercolor(c,:))
    legend_str{c} = sprintf('%2.2f < e < %2.2f', quiverstep(c), quiverstep(c+1));
end
xlabel('x relative [m]');
ylabel('y relative [m]');
zlabel('z relative [m]');
legend(legend_str,'location','bestoutside');
    %}
clear interval vec_x vec_y vec_z vec_X vec_Y vec_Z vec_xx vec_yy vec_zz
clear vec_xxyyzz n_vec vec_uuvvww j vec_t
clear max_mag n_quivercolor quiverstep quivercolor legend_str is
% Outputs
% vector field figure
    %}
    %% Per-point Error
    %{
filepath_plyte = cell(n_S,1);
for j = 1:n_S;
    fprintf('\nWriting ply %g of %g\n',j,n_S);
    filepath_plyte{j} = sprintf('%ste.ply',filepath_ply{j}(1:end-4));
    [vertex, ~] = read_ply(filepath_ply{j});
    data_x = vertex(:,1);
    data_y = vertex(:,2);
    data_z = vertex(:,3);
    n_data = numel(data_x);
    data2_xyz = (G_R{j}*[data_x data_y data_z]') + ...
        repmat(G_t{j},1,n_data);
    data12_xyz = (G12_R{j}*[data_x data_y data_z]') + ...
        repmat(G12_T{j},1,n_data);
    data_e = sqrt(sum((data12_xyz - [data_x data_y data_z]').^2,1));
    color = plyintensity2color(data_e, 'jet');
    if max(color(:))<.001;
        color = repmat([ 0 0 127.5], n_data,1);
    end
    write2ply(filepath_plyte{j},data2_xyz',color);
end
clear j vertex dataJ_x dataJ_y dataJ_z n_data
clear data2_xyz data12_xyz data_e color
% Outputs
% filepath_plyte        filepath for ply transformed with color = error
% PLY with color according to error written to disk
    %}
    %% Per-point Vector Field
    %{
j = 2;
[vertex, ~] = read_ply(filepath_ply{j});
data_x = vertex(:,1);
data_y = vertex(:,2);
data_z = vertex(:,3);
n_data = numel(data_x);
data12_xyz = (G12_R{j}*[data_x data_y data_z]') + ...
    repmat(G12_T{j},1,n_data);
data_uvw = data12_xyz - [data_x data_y data_z]';
data_e = sqrt(sum((data_uvw).^2,1));

% Plot
figure;
hold on
max_mag = max(data_e);
n_quivercolor = 25;
quiverstep = linspace(0,max_mag,n_quivercolor);
quivercolor = jet(n_quivercolor-1);
legend_str = cell(n_quivercolor-1,1);
for c = 1:n_quivercolor-1;
    is = data_e>=quiverstep(c) & data_e<=quiverstep(c+1);
    quiver3(data_x(is), data_y(is), data_z(is),...
        data_uvw(1,is)', data_uvw(2,is)', data_uvw(3,is)',...
        'color', quivercolor(c,:))
    legend_str{c} = sprintf('%2.2f < e < %2.2f', quiverstep(c), quiverstep(c+1));
end
xlabel('x relative [m]');
ylabel('y relative [m]');
zlabel('z relative [m]');
legend(legend_str,'location','bestoutside');
%
clear j vertex data_x data_y data_z n_data data12_xyz
clear data_uv data_e n_quivercolor quiverstep quivercolor legend_str
clear c is
% Outputs
% Outputs
% per point vector field figure
    %}
    
end



