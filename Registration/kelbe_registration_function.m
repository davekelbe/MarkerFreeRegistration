function [ match_R, match_t ] = kelbe_registration_function( aux )
%UNTITLED Summary of this function goes here
%   This is identical to lines 1-1275 of kelbe_registration 
match_R = [];
match_t = [];

P_LCS = aux.P_LCS;
%P_color = aux.P_color; 
%P_plot = aux.P_plot; 
P_rad = aux.P_rad; 
%filepath_ply = aux.filepath_ply;
%info_exp = aux.info_exp; 
%info_plot = aux.info_plot; 
%info_site = aux.info_site; 
%info_slash = aux.info_slash;
%info_suffix = aux.info_suffix;
%info_valid_plot = aux.info_valid_plot;
n_S = aux.n_S;
n_tree = aux.n_tree;
%options_axesval = aux.options_axesval;
%options_imagepoints = aux.options_imagepoints;        
%options_initialmatch = aux.options_initialmatch;           
%options_loadmatch = aux.options_loadmatch;            
%options_loadvar = aux.options_loadvar;         
%options_unique = aux.options_unique;                
%options_verbose = aux.options_verbose;           
%p_std = aux.p_std;                          
%path_mat = aux.path_mat;          
%path_ply = aux.path_ply;             
%path_save = aux.path_save;       
%path_site = aux.path_site;            
%path_tikz = aux.path_tikz;    
t_RANSAC_nsearch = aux.t_RANSAC_nsearch; 
%t_RANSAC_rad = aux.t_RANSAC_rad;       
t_RANSAC_xyz = aux.t_RANSAC_xyz;          
t_coll = aux.t_coll;         
t_eig_error = aux.t_eig_error;    
t_rad = aux.t_rad;

    %% Determine disjoint input sets
    %{
    P_ix1 = cell(n_S,1);
    P_ix2 = cell(n_S,1);
    
    for s = 1:n_S;
        ix_randi = randperm(n_tree(s));
        mid = floor(n_tree(s)/2);
        P_ix1{s} = ix_randi(1:mid);
        P_ix2{s} = ix_randi(mid + 1:end);
    end
    
    n_ix1 = cellfun(@numel,P_ix1);
    n_ix2 = cellfun(@numel,P_ix2);
    
    % Disjoint sets
    if false;%options_verbose && options_imagepoints;
        s=1;
        clear legend_str
        figure
        hold on
        plot3(0,0,0,'^k','markersize',10,...
            'markerfacecolor',color(s,:));
        hdummy = plot3(P_LCS{s}(1,P_ix1{s}),P_LCS{s}(2,P_ix1{s}),P_LCS{s}(3,P_ix1{s}),'ok','markersize',5,...
            'markerfacecolor','r');
        set(hdummy, 'visible', 'off');
        hdummy = plot3(P_LCS{s}(1,P_ix2{s}),P_LCS{s}(2,P_ix2{s}),P_LCS{s}(3,P_ix2{s}),'ok','markersize',5,...
            'markerfacecolor','b');
        set(hdummy, 'visible', 'off');
        hdummy = plot([P_LCS{s}(1,1) P_LCS{s}(1,2)],...
            [P_LCS{s}(2,1) P_LCS{s}(2,2)],'-r')
        set(hdummy, 'visible', 'off');
        hdummy = plot([P_LCS{s}(1,1) P_LCS{s}(1,2)],...
            [P_LCS{s}(2,1) P_LCS{s}(2,2)],'-b')
        set(hdummy, 'visible', 'off');
        for t = 1:numel(P_ix1{s});
            h = filledCircle([P_LCS{s}(1,P_ix1{s}(t)); P_LCS{s}(2,P_ix1{s}(t))]',P_rad{s}(P_ix1{s}(t)),1000,'r');
        end
        for t = 1:numel(P_ix2{s});
            h = filledCircle([P_LCS{s}(1,P_ix2{s}(t)); P_LCS{s}(2,P_ix2{s}(t))]',P_rad{s}(P_ix2{s}(t)),1000,'b');
        end
        %scatter3(P_LCS{s}(1,:),P_LCS{s}(2,:),P_LCS{s}(3,:),30,...
        %    color_P_index(truth_P_index{s},:),'filled');
        %axis(1.5*[i_xmin i_xmax i_ymin i_ymax -10 10]);
        axis auto
        axisval = axis;
        %set(gca, 'xtick',
        xlabel('x Position relative to plot center [m]');
        ylabel('y Position relative to plot center [m]');
        zlabel('z Position relative to plot center [m]');    view(0,90);
        grid on
        %titlestr = sprintf('Scan %g',P_plot(s));
        %title(titlestr);
        legend_str{1} = sprintf('Scanner %g',P_plot(s));
        legend_str{2} = sprintf('Disjoint set 1');%, P_plot(s));
        legend_str{3} = sprintf('Disjoint set 2');%, P_plot(s));
        legend(legend_str,'location','northeast');
        legend boxoff       % Hides the legend's axes (legend border and background)
        set(gca, 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        filepath_save = sprintf('%sLCS_%02.0f.eps',path_save, P_plot(s));
        saveas(gcf,filepath_save,'psc2')
    end
    
    clear s hdummy t legend_str ix_randi mid
    
    % Outputs
    % P_ix1                 index to points in set 1
    % P_ix2                 index to points in set 2
    % n_ix1                 number of points in set 1
    % n_ix2                 number of points in set 2
    %}
    %% Find radius similarity
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
    %% Find connection network
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
    
    % Plot connection network
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
    %% Find combination network for the disjoint sets
    %{
    T1_combis = cell(n_S,1);
    T2_combis = cell(n_S,1);
    
    for s = 1:n_S;
        % Valid triangles are rows which do not include any trees from opposite
        % set
        isvalid = false(n_comb(s),n_ix2(s));
        for i = 1:n_ix2(s);
            isvalid(:,i) = any(T_comb{s}==P_ix2{s}(i),2);
        end
        T1_combis{s} = ~any(isvalid,2);
    end
    for s = 1:n_S;
        % Repeat
        isvalid = false(n_comb(s),n_ix1(s));
        for i = 1:n_ix1(s);
            isvalid(:,i) = any(T_comb{s}==P_ix1{s}(i),2);
        end
        T2_combis{s} = ~any(isvalid,2);
    end
    n_comb1 = cellfun(@sum,T1_combis);
    n_comb2 = cellfun(@sum,T2_combis);
    
    T1_comb = cell(n_S,1);
    T2_comb = cell(n_S,1);
    
    for s = 1:n_S;
        T1_comb{s} = T_comb{s}(T1_combis{s},:);
        T2_comb{s} = T_comb{s}(T2_combis{s},:);
    end
    
    %{
i = s;
hold on
for c = 1:200;% n_T1_comb;
    plot([P_LCS{s}(1,T1_comb{i}(c,1)) P_LCS{i}(1,T1_comb{i}(c,2))],...
         [P_LCS{s}(2,T1_comb{i}(c,1)) P_LCS{i}(2,T1_comb{i}(c,2))],'-r')
    plot([P_LCS{s}(1,T1_comb{i}(c,1)) P_LCS{i}(1,T1_comb{i}(c,3))],...
         [P_LCS{s}(2,T1_comb{i}(c,1)) P_LCS{i}(2,T1_comb{i}(c,3))],'-r')
    plot([P_LCS{s}(1,T1_comb{i}(c,2)) P_LCS{i}(1,T1_comb{i}(c,3))],...
         [P_LCS{s}(2,T1_comb{i}(c,2)) P_LCS{i}(2,T1_comb{i}(c,3))],'-r')
end
for c = 1:200;% n_T2_comb;
    plot([P_LCS{s}(1,T2_comb{i}(c,1)) P_LCS{s}(1,T2_comb{i}(c,2))],...
         [P_LCS{s}(2,T2_comb{i}(c,1)) P_LCS{s}(2,T2_comb{i}(c,2))],'-b')
    plot([P_LCS{s}(1,T2_comb{i}(c,1)) P_LCS{s}(1,T2_comb{i}(c,3))],...
         [P_LCS{s}(2,T2_comb{i}(c,1)) P_LCS{s}(2,T2_comb{i}(c,3))],'-b')
    plot([P_LCS{s}(1,T2_comb{i}(c,2)) P_LCS{s}(1,T2_comb{i}(c,3))],...
         [P_LCS{s}(2,T2_comb{i}(c,2)) P_LCS{s}(2,T2_comb{i}(c,3))],'-b')
end
    legend_str{4} = sprintf('Disjoint connections 1');
    legend_str{5} = sprintf('Disjoint connections 2');
    legend(legend_str,'location','northeast');

    %}
    
    clear i s c isvalid
    % Outputs:
    %       T1_combis{s}              - Logical array of valid triangles for set 1
    %       T1_comb{s}                - Triangle combinations of disjoint set 1
    %       n_comb1(s)               - number of triangles for disjoint set 1
    %       T2_combis{s}              - Logical array of valid triangles for set 2
    %       T2_comb{s}                - Triangle combinations of disjoint set 2
    %       n_comb2(s)               - number of triangles for disjoint set 2
    %}
    %% Sort triangle set combinations by decreasing radius
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
    %% Find radii for disjoint sets
    %{
    T1_LCS = cell(n_S,1);
    T1_R = cell(n_S,1);
    T2_LCS = cell(n_S,1);
    T2_R = cell(n_S,1);
    
    for s = 1:n_S;
        T1_LCS{s} = T_LCS{s}(T1_combis{s},:,:); %Already sorted by radius
        T1_R{s} = T_R{s}(T1_combis{s},:);
        % Repeat
        T2_LCS{s} = T_LCS{s}(T2_combis{s},:,:);
        T2_R{s} = T_R{s}(T2_combis{s},:);
    end
    clear s
    % Outputs:
    %       T1_LCS(i)               - Triangle points for disjoint set 1
    %       T1_R{i}                 - Triangle radii for disjoint set 1
    %       T2_LCS(i)               - Triangle points for disjoint set 2
    %       T2_R{i}                 - Triangle radii for disjoint set 2
    %}
    %% Find eigenvalues of each triangle
    
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
    % Histogram of normalized (smaller) eigenvalues
    %{
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
    %% Find eigenvalues for disjoint sets
    %{
    T1_eig = cell(n_S,1);
    T1_isncoll = cell(n_S,1);
    T2_eig = cell(n_S,1);
    T2_isncoll = cell(n_S,1);
    
    for s = 1:n_S;
        T1_eig{s} = T_eig{s}(T1_combis{s},:);
        T1_isncoll{s} = T_isncoll{s}(T1_combis{s});
        % Repeat
        T2_eig{s} = T_eig{s}(T2_combis{s},:);
        T2_isncoll{s} = T_isncoll{s}(T2_combis{s});
    end
    clear s
    % Outputs:
    %       T1_eig(i)               - Triangle eigenvalues for disjoint set 1
    %       T1_isncoll(i)           - Not collinear for disjoint set 1
    %       T2_eig(i)               - Triangle eigenvalues for disjoint set 2
    %       T2_isncoll(i)           - Not collinear for disjoint set 2
    
    %}
    %% Find likely RANSAC pairs
    % Create 1D array of valid triangles after filtering by collinearity
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
    %% Find likely RANSAC pairs by filtering eigenvalues, radius
    % Could also utilize 4-point combinations, graph theory, etc.
    % Also add radius information
    
    % More intuitive way: to check consistency
    %{
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
    clear match_israd match_rad_diff
    clear m_eig m_combix n_m
    clear M_eigaI M_eigaJ M_eigbI M_eigbJ M_error is_brute
    clear row col xx_combixI xx_combixJ xx_combI xx_combJ n_brute
    clear xx_israd linear xxM_error xxM_error_sort sortix r i j
    % M_bestixI{i,j}          - index of best i triangles for j-> i
    % M_bestixJ{i,j}          - index of best j triangles for j-> i
    % M_combI{i,j}            - triangle combinations
    % M_combJ{i,j}            - triangle combinations
    % M_eigerror{i,j}         - Eigenvalue error for each triangle match
    %% Construct ordered tri radius, locations for clarity
    M_triI = cell(n_S,n_S);
    M_triJ = cell(n_S,n_S);
    M_triradI = cell(n_S,n_S);
    M_triradJ = cell(n_S,n_S);
    for i = 1:n_S;
        for j = i + 1:n_S;
            M_triI{i,j} = zeros(n_M(i,j),3,3);
            M_triI{i,j}(:,1,:) = P_LCS{i}(:,M_combI{i,j}(:,1))';
            M_triI{i,j}(:,2,:) = P_LCS{i}(:,M_combI{i,j}(:,2))';
            M_triI{i,j}(:,3,:) = P_LCS{i}(:,M_combI{i,j}(:,3))';
            M_triradI{i,j} = P_rad{i}(M_combI{i,j});
            M_triJ{i,j} = zeros(n_M(i,j),3,3);
            M_triJ{i,j}(:,1,:) = P_LCS{j}(:,M_combJ{i,j}(:,1))';
            M_triJ{i,j}(:,2,:) = P_LCS{j}(:,M_combJ{i,j}(:,2))';
            M_triJ{i,j}(:,3,:) = P_LCS{j}(:,M_combJ{i,j}(:,3))';
            M_triradJ{i,j} = P_rad{j}(M_combJ{i,j});
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
    %% Unnecessary figures
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
    %% Permutations??  OLD CODE
    %{
perm_of_3 = perms(1:3);
n_perm_of_3 = size(perm_of_3,1);
perm_of_3a = perm_of_3;
perm_of_3b = repmat(perm_of_3(1,:),[n_perm_of_3,1]);
n_p = size(perm_of_3a,1);
    %}
    %% Find pairwise matches between full sets
    % fprintf('\nFind pairwise matches \n');
    
    % Initialize graph and corresponding functional relationships
    %{
    match_i = cell(n_S); %base
    match_j = cell(n_S); %mobile
    match_R = cell(n_S);
    match_t = cell(n_S);
    %}
    
    % Find rotation and tranlation from j to i along with point matching pairs
    for i = 1:n_S-1;
        for j = i+1:n_S; % Symmetric
            % fprintf('\n\tMatching %g to %g\n',j,i);
            
            % RANSAC setup
            n_match_best = 0;
            m1_best = [];
            m2_best = [];
            n_search = min(t_RANSAC_nsearch, n_M(i,j));
            
            for r = 1:n_search;
                
                triI = squeeze(M_triI{i,j}(r,:,:));
                triJ = squeeze(M_triJ{i,j}(r,:,:));
                
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
                dataJt = ((Rhat*P_LCS{j})+ repmat(that,1,n_tree(j)));
                % Add rotation information
                % Use n_match, round-trip or other method to determine fit
                P_LCSIrep = repmat(P_LCS{i}',[1,1,n_tree(j)]);
                P_LCSJtrep = zeros(n_tree(i),3,n_tree(j));
                for t = 1:n_tree(i);
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
                    %m1_best = m1;
                    %m2_best = m2;
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
                match_R = R_best;
                match_t = t_best;
                %match_i{i,j} = m1_best;
                %match_j{i,j} = m2_best;
            end
        end
    end
    
    % Declare match from i-i identity
    %{
    for i = 1:n_S;
        match_R{i,i} = eye(3);
        match_t{i,i} = zeros(3,1);
        match_i{i,i} = 1:numel(P_LCS{i});
        match_j{i,i} = 1:numel(P_LCS{i});
    end
    %}
    
    % Make R,t non-directed
    %{
    for i = 1:n_S;
        for j = i:n_S;
            if ~isempty(match_R{i,j});
                match_R{j,i} = match_R{i,j}';
                match_t{j,i} = -(match_R{i,j}')*match_t{i,j};
            end
        end
    end
    %}
    
    clear n_match_best m1_best m2_best n_search i j r Rhat that dataJt
    clear P_LCSIrep P_LCSJtrep t error_xyz m1 m2 n_match
    % Outputs
    % match_R               {i,j} is pairwise rotation from j into i
    % match_t               {i,j} is pairwise translation from j into i
    % match_i               {i,j} is pairwise matches from i
    % match_j               {i,j} is pairwise matches from j
    %% Update M array for disjoint sets
    %{
    M1_combis = cell(n_S,n_S);
    M2_combis = cell(n_S,n_S);
    M1_combixI = cell(n_S,n_S);
    M1_combixJ = cell(n_S,n_S);
    M2_combixI = cell(n_S,n_S);
    M2_combixJ = cell(n_S,n_S);
    
    M1_combI = cell(n_S,n_S);
    M1_combJ = cell(n_S,n_S);
    M2_combI = cell(n_S,n_S);
    M2_combJ = cell(n_S,n_S);
    M1_triI = cell(n_S,n_S);
    M1_triJ = cell(n_S,n_S);
    M1_triradI = cell(n_S,n_S);
    M1_triradJ = cell(n_S,n_S);
    M2_triI = cell(n_S,n_S);
    M2_triJ = cell(n_S,n_S);
    M2_triradI = cell(n_S,n_S);
    M2_triradJ = cell(n_S,n_S);
    n_M1 = zeros(n_S,n_S);
    n_M2 = zeros(n_S,n_S);
    
    T1_combix = cell(n_S,1);
    T2_combix = cell(n_S,1);
    
    for i = 1:n_S;
        T1_combix{i} = find(T1_combis{i});
        T2_combix{i} = find(T2_combis{i});
    end
    
    % Make M_bestixI and M_bestixJ symmetric
    for i = 1:n_S;
        for j = i+1:n_S;
            M_bestixJ{j,i} = M_bestixI{i,j};
            M_bestixI{j,i} = M_bestixJ{i,j};
        end
    end
    
    for i = 1:n_S;
        for j = 1:n_S;
            if i==j;
                continue
            end
            % Only one set is disjoint
            M1_combis{i,j} = ismember(M_bestixI{i,j}, T1_combix{i});
            M2_combis{i,j} = ismember(M_bestixJ{i,j}, T2_combix{j});
            % both are disjoint *Also change j:i+1
            %         M1_combis{i,j} = ismember(M_bestixI{i,j},T1_combix{i}) & ...
            %             ismember(M_bestixJ{i,j}, T1_combix{j});
            %         M2_combis{i,j} = ismember(M_bestixI{i,j},T2_combix{i}) & ...
            %             ismember(M_bestixJ{i,j}, T2_combix{j});
            M1_combixI{i,j} = M_bestixI{i,j}(M1_combis{i,j});
            M1_combixJ{i,j} = M_bestixJ{i,j}(M1_combis{i,j});
            M2_combixI{i,j} = M_bestixI{i,j}(M2_combis{i,j});
            M2_combixJ{i,j} = M_bestixJ{i,j}(M2_combis{i,j});
            M1_combI{i,j} = T_comb{i}(M1_combixI{i,j},:);
            M1_combJ{i,j} = T_comb{j}(M1_combixJ{i,j},:);
            M2_combI{i,j} = T_comb{i}(M2_combixI{i,j},:);
            M2_combJ{i,j} = T_comb{j}(M2_combixJ{i,j},:);
            n_M1(i,j) = size(M1_combI{i,j},1);
            n_M2(i,j) = size(M2_combI{i,j},1);
            M1_triI{i,j} = zeros(n_M1(i,j),3,3);
            M1_triI{i,j}(:,1,:) = P_LCS{i}(:,M1_combI{i,j}(:,1))';
            M1_triI{i,j}(:,2,:) = P_LCS{i}(:,M1_combI{i,j}(:,2))';
            M1_triI{i,j}(:,3,:) = P_LCS{i}(:,M1_combI{i,j}(:,3))';
            M1_triradI{i,j} = P_rad{i}(M1_combI{i,j});
            M1_triJ{i,j} = zeros(n_M1(i,j),3,3);
            M1_triJ{i,j}(:,1,:) = P_LCS{j}(:,M1_combJ{i,j}(:,1))';
            M1_triJ{i,j}(:,2,:) = P_LCS{j}(:,M1_combJ{i,j}(:,2))';
            M1_triJ{i,j}(:,3,:) = P_LCS{j}(:,M1_combJ{i,j}(:,3))';
            M1_triradJ{i,j} = P_rad{j}(M1_combJ{i,j});
            M2_triI{i,j} = zeros(n_M2(i,j),3,3);
            M2_triI{i,j}(:,1,:) = P_LCS{i}(:,M2_combI{i,j}(:,1))';
            M2_triI{i,j}(:,2,:) = P_LCS{i}(:,M2_combI{i,j}(:,2))';
            M2_triI{i,j}(:,3,:) = P_LCS{i}(:,M2_combI{i,j}(:,3))';
            M2_triradI{i,j} = P_rad{i}(M2_combI{i,j});
            M2_triJ{i,j} = zeros(n_M2(i,j),3,3);
            M2_triJ{i,j}(:,1,:) = P_LCS{j}(:,M2_combJ{i,j}(:,1))';
            M2_triJ{i,j}(:,2,:) = P_LCS{j}(:,M2_combJ{i,j}(:,2))';
            M2_triJ{i,j}(:,3,:) = P_LCS{j}(:,M2_combJ{i,j}(:,3))';
            M2_triradJ{i,j} = P_rad{j}(M2_combJ{i,j});
        end
    end
    
    clear i jM_isI M_isJ M1_combisI M1_combisJ M2_combisI M2_combisJ
    clear M1_combI M1_combJ M2_combI M2_combJ
    % Outputs:
    % M1_triI{i,j}             - triangle of points for set 1I
    % M1_triJ{i,j}             - triangle of points for set 1J
    % M1_triradI{i,j}          - radius of triangle of points for set 1I
    % M1_triradJ{i,j}          - radius of triangle of points for set 1J
    % M2_triI{i,j}             - triangle of points for set 1I
    % M2_triJ{i,j}             - triangle of points for set 1J
    % M2_triradI{i,j}          - radius of triangle of points for set 1I
    % M2_triradJ{i,j}          - radius of triangle of points for set 1J
    %}
    %% Find pairwise matches between disjoint sets for set 1
    %{
    % Set 1
    %save('line962.mat');
    %     else
    %         load('line962.mat');
    %     end
    
    fprintf('\nFind pairwise matches \n');
    
    % Initialize graph and corresponding functional relationships
    match1_i = cell(n_S); %base
    match1_j = cell(n_S); %mobile
    match1_R = cell(n_S);
    match1_t = cell(n_S);
    
    % Find rotation and tranlation from j to i along with point matching pairs
    for i = 1:n_S;
        for j = 1:n_S; % Not symmetric
            fprintf('\n\tMatching %g to %g\n',j,i);
            if i==j;
                continue
            end
            % RANSAC setup
            n_match_best = 0;
            m1_best = [];
            m2_best = [];
            
            for r = 1:n_M1(i,j);
                
                triI = squeeze(M1_triI{i,j}(r,:,:));
                triJ = squeeze(M1_triJ{i,j}(r,:,:));
                
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
                dataJt = ((Rhat*P_LCS{j})+ repmat(that,1,n_tree(j)));
                % Add rotation information
                % Use n_match, round-trip or other method to determine fit
                P_LCSIrep = repmat(P_LCS{i}',[1,1,n_tree(j)]);
                P_LCSJtrep = zeros(n_tree(i),3,n_tree(j));
                for t = 1:n_tree(i);
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
                end
                plot_triangle(triI);
                plot_triangle(triJt);
                plot_triangle_connection(triI,triJt);
                %}
            end
            
            if n_match_best > 0;
                match1_R{i,j} = R_best;
                match1_t{i,j} = t_best;
                match1_i{i,j} = m1_best;
                match1_j{i,j} = m2_best;
            end
        end
    end
    
    % Declare match from i-i identity
    for i = 1:n_S;
        match1_R{i,i} = eye(3);
        match1_t{i,i} = zeros(3,1);
        match1_i{i,i} = 1:numel(P_LCS{i});
        match1_j{i,i} = 1:numel(P_LCS{i});
    end
    
    % Make R,t non-directed
    %     for i = 1:n_S;
    %         for j = i:n_S;
    %             if ~isempty(match1_R{i,j});
    %                 match1_R{j,i} = match1_R{i,j}';
    %                 match1_t{j,i} = -(match1_R{i,j}')*match1_t{i,j};
    %             end
    %         end
    %     end
    
    clear filepath_match_R filepath_match_t filepath_match_i filepath_match_j
    clear n_match_best Rhat that m1_best m2_best n_search
    clear r i j h n_match_best triI triJ Rhat that triJt dataJt
    clear  P_LCSIrep t P_LCSJtrep  error_xyz m1 m2 n_match
    clear R_best t_best
    % Outputs
    % match1_R               (i,j) is pairwise rotation from j into i
    % match1_t               (i,j) is pairwise translation from j into i
    % match1_i               (i,j) is pairwise matches from i
    % match1_j               (i,j) is pairwise matches from j
    %}
    %% Find pairwise matches between disjoint sets for set 2
    %{
    % Set 1
    %     %save('line1100.mat');
    %      else
    %          load('line1100.mat');
    %      end
    
    options_rotation = false;
    fprintf('\nFind pairwise matches \n');
    
    % Initialize graph and corresponding functional relationships
    match2_i = cell(n_S); %base
    match2_j = cell(n_S); %mobile
    match2_R = cell(n_S);
    match2_t = cell(n_S);
    
    % Find rotation and tranlation from j to i along with point matching pairs
    for i = 1:n_S; %why -1 before?
        for j = 1:n_S; % Symmetric
            fprintf('\n\tMatching %g to %g\n',j,i);
            
            if i==j;
                continue
            end
            % RANSAC setup
            n_match_best = 0;
            m1_best = [];
            m2_best = [];
            
            for r = 1:n_M2(i,j);
                
                triI = squeeze(M2_triI{i,j}(r,:,:));
                triJ = squeeze(M2_triJ{i,j}(r,:,:));
                
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
                dataJt = ((Rhat*P_LCS{j})+ repmat(that,1,n_tree(j)));
                % Add rotation information
                % Use n_match, round-trip or other method to determine fit
                P_LCSIrep = repmat(P_LCS{i}',[1,1,n_tree(j)]);
                P_LCSJtrep = zeros(n_tree(i),3,n_tree(j));
                for t = 1:n_tree(i);
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
                end
                plot_triangle(triI);
                plot_triangle(triJt);
                plot_triangle_connection(triI,triJt);
                %}
            end
            
            if n_match_best > 0;
                match2_R{i,j} = R_best;
                match2_t{i,j} = t_best;
                match2_i{i,j} = m1_best;
                match2_j{i,j} = m2_best;
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
    
    % Flip R,t to go the opposite direction (J->I)
    match2_Rflip = cell(n_S,n_S);
    match2_tflip = cell(n_S,n_S);
    
    % Make R,t non-directed
    for i = 1:n_S;
        for j = 1:n_S;
            match2_Rflip{i,j} = match2_R{i,j}';
            match2_tflip{i,j} = -(match2_R{i,j}')*match2_t{i,j};
        end
    end
    match2_R = match2_Rflip;
    match2_t = match2_tflip;
    
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
    %
    %   save('line1238.mat');
    %  else
    %      load('line1238.mat');
    %  end
    %}
    %% Determine RMSE error of stem centers using circular path 
    %{
    match12_R = cellfun(@(x,y) x*y,match2_R, match1_R, 'uniformOutput', false);
    match12_t = cellfun(@(x,y,z) (x*y) + z,match2_R, match1_t, match2_t, 'uniformOutput', false);
    match12_rx = 180*cellfun(@decompose_rotation_rx,match12_R)/pi;
    match12_ry = 180*cellfun(@decompose_rotation_ry,match12_R)/pi;
    match12_rz = 180*cellfun(@decompose_rotation_rz,match12_R)/pi;
    match12_er = abs(match12_rx) + abs(match12_ry) + abs(match12_rz);
    %match12_et = cellfun(@sum,cellfun(@abs,match12_t,'uniformOutput',false));
    %match12_etxy = cellfun(@sum,cellfun(@abs,cellfun(@(x) x(1:2),match12_t, 'uniformoutput', false),'uniformOutput',false));
    
    % OR, could generate synthetic points about volume and R1t1R2t2 them
    match12_RMSE = zeros(n_S,n_S);
    for i = 1:n_S;
        for j = 1:n_S;
            P_LCSt =  (match12_R{i,j}*P_LCS{i})+ repmat(match12_t{i,j},1,n_tree(i));
            P_dist = sqrt(sum((P_LCSt - P_LCS{i}).^2));
            match12_RMSE(i,j) = mean(P_dist);
        end
    end
    clear i j P_dist P_LCSt
    % Outputs
    % match12_R               {i,j} Effective rotation
    % match12_t               {i,j} Effective translation
    %   save('line1309.mat');
    %else
    %    load('line1309.mat');
    %end
    %}
end

