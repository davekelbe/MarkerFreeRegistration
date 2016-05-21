function [ Rhat, that, m1_best, m2_best ] = toy_registrationfunction2(tree1_xyz, tree2_xyz, tree1_r, tree2_r, israd )
%UNTITLED8 Summary of this function goes here
%   Finds the transformation T which maps plot2 into plot1

%% Initialization

Rhat = [];
that = [];
m1_best = [];
m2_best = [];

% Set Internal Matlab Parameters
set(0,'DefaultFigureRenderer','OpenGL')
set(0,'defaultfigureposition',[966  322  557 472]);
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:legend:IgnoringExtraEntries');
set(0,'DefaultFigureColor', 'White');
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontName', 'Calibri')
set(0,'DefaultAxesFontSize', 14)
%%
% Radius difference
t_rad = 0.1;
% Eigenvalue difference
%t_eig = inf - just sorting now
% Distance between tree centers for match in RANSAC
t_xyz = 0.4.^2;
% Number of RANSAC iterations
n_search = 16;
options_rotation = false;


n_tree1 = numel(tree1_r);
n_tree2 = numel(tree2_r);

% Outputs:
%       tree1_r               - radius vector for set 1
%       tree1_xyz             - xyz vectors for set 1
%       tree2_r               - radius vector for set 2
%       tree2_xyz             - xyz vectors for set 2
clear tree1 tree2 tree i filepath_tree1 filepath_tree2
%% Plot DBH maps
% Future Work? Get initial transformation for plotting
%{
figure;
hold on
h1 = plot3(0,0,30,'ok','markerfacecolor',[.57 1 .53]);
h2 = plot3(0,0,30,'ok','markerfacecolor',[.96 .46 .50]);
set(h1, 'visible','off');
set(h2, 'visible','off');
for s = 1:n_tree1;
    h = filledCircle([tree1_xyz(s,1) tree1_xyz(s,2)],tree1_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_tree2;
    h = filledCircle([tree2_xyz(s,1) tree2_xyz(s,2)],tree2_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on
legend_str{1} = 'Tree i';
legend_str{2} = 'Tree j';
legend(legend_str);
%}
% Outputs:
clear h h1 h2 s
%% Find connection network

comb1 = combnk(1:n_tree1,3);
comb2 = combnk(1:n_tree2,3);
n_comb1 = size(comb1,1);
n_comb2 = size(comb2,1);

%{
for c = 1:200;% n_comb1; 
    plot([tree1_xyz(comb1(c,1),1) tree1_xyz(comb1(c,2),1)],...
         [tree1_xyz(comb1(c,1),2) tree1_xyz(comb1(c,2),2)],'-k')
    plot([tree1_xyz(comb1(c,1),1) tree1_xyz(comb1(c,3),1)],...
         [tree1_xyz(comb1(c,1),2) tree1_xyz(comb1(c,3),2)],'-k')
    plot([tree1_xyz(comb1(c,2),1) tree1_xyz(comb1(c,3),1)],...
         [tree1_xyz(comb1(c,2),2) tree1_xyz(comb1(c,3),2)],'-k')
end
for c = 1:200;% n_comb2; 
    plot([tree2_xyz(comb2(c,1),1) tree2_xyz(comb2(c,2),1)],...
         [tree2_xyz(comb2(c,1),2) tree2_xyz(comb2(c,2),2)],'-k')
    plot([tree2_xyz(comb2(c,1),1) tree2_xyz(comb2(c,3),1)],...
         [tree2_xyz(comb2(c,1),2) tree2_xyz(comb2(c,3),2)],'-k')
    plot([tree2_xyz(comb2(c,2),1) tree2_xyz(comb2(c,3),1)],...
         [tree2_xyz(comb2(c,2),2) tree2_xyz(comb2(c,3),2)],'-k')
end
%}
% Outputs:
%       comb1                 - Combinations of tree1's; n_comb1 x 3
%       bomb2                 - Combinations of tree2's; n_comb2 x 3
%       n_comb1               - number of combinations 1
%       n_comb2               - number of combinations 2
%% Filter combination network by radius similarity 


%% Find radius


rad1 = zeros(n_comb1,3);
for i = 1:n_comb1;
    rad1(i,:) = tree1_r(comb1(i,:))';
end
[rad1s,ixsort] = sort(rad1,2,'descend');
comb1s = zeros(n_comb1,3);
for i = 1:n_comb1;
    comb1s(i,:) = comb1(i,ixsort(i,:));
end

rad2 = zeros(n_comb2,3);
for i = 1:n_comb2;
    rad2(i,:) = tree2_r(comb2(i,:))';
end
[rad2s,ixsort] = sort(rad2,2,'descend');
comb2s = zeros(n_comb2,3);
for i = 1:n_comb2;
    comb2s(i,:) = comb2(i,ixsort(i,:));
end
% Outputs:
%       comb1s                - Sort Combinations based on radius
%       comb2s                - Sort Combinations based on radius
%       rad1s             - Radius for each
%       rad2s              - Radisu for each
clear i ixsort rad1 rad2
%% Find difference in radius
%{
rep_tree1r = repmat(tree1_r,[1,n_tree2]);
rep_tree2r = repmat(tree2_r',[n_tree1,1]);
rep_rdiff = abs((rep_tree1r - rep_tree2r)./(rep_tree1r + rep_tree2r));
rep_isrdiff = (rep_rdiff < t_rad);

is_valid = false(n_comb1, n_comb2);

for i = 1:n_comb1;
    rep_c1i = repmat(comb1s(i,:),[n_comb2,1]);
    comb1s_1d = reshape(rep_c1i,1,[]);
    comb2s_1d = reshape(comb2s,1,[]);
    linear_ix = sub2ind([n_tree1, n_tree2],comb1s_1d,comb2s_1d);
    ismatch_1d = rep_isrdiff(linear_ix);
    is_match = any(reshape(ismatch_1d,3,[]));
    is_valid(i,is_match) = true;
end
    %}

% Method 1 
%{
tic 
rep_ra1 = repmat(rad1s(:,1),[1,n_comb2]);
rep_ra2 = repmat(rad2s(:,1)',[n_comb1,1]);
error_ra = abs((rep_ra1-rep_ra2)./((rep_ra1 + rep_ra2)/2));
%figure; imagesc(error_ra<t_rad); title('A1');
rep_rb1 = repmat(rad1s(:,2),[1,n_comb2]);
rep_rb2 = repmat(rad2s(:,2)',[n_comb1,1]);
error_rb = abs((rep_rb1-rep_rb2)./((rep_rb1 + rep_rb2)/2));
%figure; imagesc((error_rb<t_rad)&(error_ra<t_rad)); title('A12');
rep_rc1 = repmat(rad1s(:,3),[1,n_comb2]);
rep_rc2 = repmat(rad2s(:,3)',[n_comb1,1]);
error_rc = abs((rep_rc1-rep_rc2)./((rep_rc1 + rep_rc2)/2));
%figure; imagesc((error_rc<t_rad)&(error_rb<t_rad)&(error_ra<t_rad)); title('A123');
is_valid1 = (error_ra<t_rad)&(error_rb<t_rad)&(error_rc<t_rad);
fprintf('\ntime = %3.2f sec\n',toc);
%clear rep_ra1 rep_ra2 error_ra rep_rb1 rep_rb2 error_rb rep_rc1 rep_rc2 error_rc
%}
%tic 
% Method 2 
% Error in the largest radius
rep_ra1 = repmat(rad1s(:,1),[1,n_comb2]);
rep_ra2 = repmat(rad2s(:,1)',[n_comb1,1]);
error_ra = abs((rep_ra1-rep_ra2)./((rep_ra1 + rep_ra2)/2));
is_ra = (error_ra<t_rad);
%figure; imagesc(is_ra); title('B1');
[row,col] = find(is_ra);
% Error in the middle radius
e_rb = (rad1s(row,2)-rad2s(col,2))./((rad1s(row,2)+rad2s(col,2))/2);
is_rb = (abs(e_rb)<t_rad);
row = row(is_rb);
col = col(is_rb);
%is_valid2 = false(n_comb1, n_comb2);
%ix = sub2ind([n_comb1 n_comb2], row,col);
%is_valid2(ix) = true;
%figure; imagesc(is_valid2); title('B12');
% Error in the smallest radius
e_rb = (rad1s(row,3)-rad2s(col,3))./((rad1s(row,3)+rad2s(col,3))/2);
is_rc = (abs(e_rb)<t_rad);
row = row(is_rc);
col = col(is_rc);
ix = sub2ind([n_comb1 n_comb2], row,col);
%fprintf('\ntime = %3.2f sec\n',toc);
is_valid = false(n_comb1, n_comb2);
is_valid(ix) = true;
%figure; imagesc(is_valid2); title('B123');

%error_r = sqrt(error_ra.^2 + error_rb.^2 + error_rc.^2);
clear error_ra error_rb error_rc 

clear rep_ra1 rep_ra2 error_ra rep_rb1 rep_rb2 error_rb rep_rc1 rep_rc2 error_rc
% Outputs:
%       is_valid              - valid combinations based on radius match
%       t_rad                 - radius threshold

%% Plotting
%{
error = error_r(is_valid);
[ix2,ix1] = meshgrid(1:n_comb2, 1:n_comb1);
index1 = ix1(is_valid);
index2 = ix2(is_valid);

[errors,Isort] = sort(error);
index1s = index1(Isort);
index2s = index2(Isort);
for i = 10001:20001
    t1 = comb1(index1s(i),:);
    t2 = comb2(index2s(i),:);
%{
    hcomba = plot([tree1_xyz(t1(1),1) tree2_xyz(t2(1),1)+ti(1)],...
        [tree1_xyz(t1(1),2) tree2_xyz(t2(1),2)+ti(2)],'-k');
    hcombb = plot([tree1_xyz(t1(2),1) tree2_xyz(t2(2),1)+ti(1)],...
        [tree1_xyz(t1(2),2) tree2_xyz(t2(2),2)+ti(2)],'-k');
    hcombc = plot([tree1_xyz(t1(3),1) tree2_xyz(t2(3),1)+ti(1)],...
        [tree1_xyz(t1(3),2) tree2_xyz(t2(3),2)+ti(2)],'-k');
%}
    htria = plot([tree1_xyz(t1(1),1) tree1_xyz(t1(2),1)],...
        [tree1_xyz(t1(1),2) tree1_xyz(t1(2),2)],'-g');
    htrib = plot([tree1_xyz(t1(2),1) tree1_xyz(t1(3),1)],...
        [tree1_xyz(t1(2),2) tree1_xyz(t1(3),2)],'-g');
    htric = plot([tree1_xyz(t1(3),1) tree1_xyz(t1(1),1)],...
        [tree1_xyz(t1(3),2) tree1_xyz(t1(1),2)],'-g');
    htrid = plot([tree2_xyz(t2(1),1) tree2_xyz(t2(2),1)] + ti(1),...
        [tree2_xyz(t2(1),2) tree2_xyz(t2(2),2)] + ti(2),'-r');
    htrie = plot([tree2_xyz(t2(2),1) tree2_xyz(t2(3),1)] + ti(1),...
        [tree2_xyz(t2(2),2) tree2_xyz(t2(3),2)] + ti(2),'-r');
    htrif = plot([tree2_xyz(t2(3),1) tree2_xyz(t2(1),1)] + ti(1),...
        [tree2_xyz(t2(3),2) tree2_xyz(t2(1),2)] + ti(2),'-r');
%{
    set(hcomba, 'visible', 'off');
    set(hcombb, 'visible', 'off');
    set(hcombc, 'visible', 'off');
        %
    set(htria, 'visible', 'off');
    set(htrib, 'visible', 'off');
    set(htric, 'visible', 'off');
    set(htrid, 'visible', 'off');
    set(htrie, 'visible', 'off');
    set(htrif, 'visible', 'off');
%}
    foo = 1;
end

%}

%% Find eigenvalues

% First and second eigenvalues for combinations 1
eigval1 = zeros(n_comb1,3);
for i = 1:n_comb1;
    data = tree1_xyz(comb1(i,:),:);
    eigval1(i,:) = eig(cov(data));
end
eigfirst1 = eigval1(:,2);
eigsecond1 = eigval1(:,3);

% First and second eigenvalues for combinations 2
eigval2 = zeros(n_comb2,3);
for i = 1:n_comb2
    data = tree2_xyz(comb2(i,:),:);
    eigval2(i,:) = eig(cov(data));
end
eigfirst2 = eigval2(:,2);
eigsecond2 = eigval2(:,3);

error = (eigfirst1 - eigfirst2).^2 +...
    (eigsecond1 - eigsecond2).^2;
[~,Isort] = sort(error);
index1s = row(Isort);
index2s = col(Isort);

%tic 
% New method 
error = (eigfirst1(row) - eigfirst2(col)).^2 +...
    (eigsecond1(row) - eigsecond2(col)).^2;
[~,Isort] = sort(error);
index1s = row(Isort);
index2s = col(Isort);
%fprintf('\n time new = %3.2f\n',toc);

%{
tic 
% Find total eigenvalue error
rep_ea1 = repmat(eigvala1, [1,n_comb2]);
rep_eb1 = repmat(eigvalb1, [1,n_comb2]);
rep_ea2 = repmat(eigvala2', [n_comb1,1]);
rep_eb2 = repmat(eigvalb2', [n_comb1,1]);
error_eigval = (rep_ea1 - rep_ea2).^2 + (rep_eb1 - rep_eb2).^2;

% Find the errors and combination indices which satisfy radius criteria
error = error_eigval(is_valid);
[ix2,ix1] = meshgrid(1:n_comb2, 1:n_comb1);
index1 = ix1(is_valid);
index2 = ix2(is_valid);

% Sort valid indexes by eigenvalue error 
[~,Isort] = sort(error);
index1s = index1(Isort);
index2s = index2(Isort);
%n_filt = numel(index1);
fprintf('\n time old = %3.2f\n',toc);
%}
% Outputs:tfa
%       errors                - Error vector for eigenvalues
%       index1s               - combination 1
%       index2s               - combination 2
%       n_filt                - Number of filtered combinations to test

clear eigval1 i data eigval2 eigvala1 eigvala2 eigbalb1 eigvalb2
clear rep_ea1 rep_eb1 rep_ea2 rep_eb2 error ix1 ix2 index1 index2

%% RANSAC

n_match_best = 0;
%R_best = eye(3);
%t_best = zeros(3,1);
m1_best = [];
m2_best = [];

perm_of_3 = perms(1:3);
n_perm_of_3 = size(perm_of_3,1);
perm_of_3a = perm_of_3;
perm_of_3b = repmat(perm_of_3(1,:),[n_perm_of_3,1]);
n_p = size(perm_of_3a,1);
n_search = min(n_search, size(index1s,1));

for i = 1:n_search;
    
    % fprintf('\n*************%g of %g**************\n', i,n_search)
    if i>numel(index1s);
        foo = 1;
    end
    if numel(index1s)<i;
        foo = 1;
    end
    t1_orig = comb1(index1s(i),:);
    t2_orig = comb2(index2s(i),:);
    
    if options_rotation
        T1 = zeros(n_p,3);
        T2 = zeros(n_p,3);
        for p = 1:n_p;
            T1(p,:) = t1_orig(perm_of_3a(p,:));
            T2(p,:) = t2_orig(perm_of_3b(p,:));
        end
    else
        perm_xyz1 = zeros(3,3,n_p);
        perm_xyz2 = zeros(3,3,n_p);
        for p = 1:n_p;
            perm_xyz1(:,:,p) = tree1_xyz(t1_orig(perm_of_3a(p,:)),:);
            perm_xyz2(:,:,p) = tree2_xyz(t2_orig(perm_of_3b(p,:)),:);
        end
        perm_error = sum((perm_xyz1-perm_xyz2).^2,2);
        perm_error = squeeze(sum(perm_error,1));
        ix_min = find(perm_error==min(perm_error));
        T1 = t1_orig(perm_of_3a(ix_min,:));
        T2 = t2_orig(perm_of_3b(ix_min,:));
    end
    
    % For each permutation, find tree1 and tree2 indexes
    for p = 1:size(T1,1);
        t1 = T1(p,:);
        t2 = T2(p,:);
        
        %{
        scatter(tree1_xyz(t1(1),1),tree1_xyz(t1(1),2),30,'r','filled')
        scatter(tree1_xyz(t1(2),1),tree1_xyz(t1(2),2),30,'g','filled')
        scatter(tree1_xyz(t1(3),1),tree1_xyz(t1(3),2),30,'b','filled')
        scatter(tree2_xyz(t2(1),1) + ti(1),tree2_xyz(t2(1),2)+ti(2),30,'r','filled')
        scatter(tree2_xyz(t2(2),1) + ti(1),tree2_xyz(t2(2),2)+ti(2),30,'g','filled')
        scatter(tree2_xyz(t2(3),1) + ti(1),tree2_xyz(t2(3),2)+ti(2),30,'b','filled')
        %}
        %{
        hcomba = plot([tree1_xyz(t1(1),1) tree2_xyz(t2(1),1)+ti(1)],...
            [tree1_xyz(t1(1),2) tree2_xyz(t2(1),2)+ti(2)],'-k');
        hcombb = plot([tree1_xyz(t1(2),1) tree2_xyz(t2(2),1)+ti(1)],...
            [tree1_xyz(t1(2),2) tree2_xyz(t2(2),2)+ti(2)],'-k');
        hcombc = plot([tree1_xyz(t1(3),1) tree2_xyz(t2(3),1)+ti(1)],...
            [tree1_xyz(t1(3),2) tree2_xyz(t2(3),2)+ti(2)],'-k');
        htria = plot([tree1_xyz(t1(1),1) tree1_xyz(t1(2),1)],...
            [tree1_xyz(t1(1),2) tree1_xyz(t1(2),2)],'-g');
        htrib = plot([tree1_xyz(t1(2),1) tree1_xyz(t1(3),1)],...
            [tree1_xyz(t1(2),2) tree1_xyz(t1(3),2)],'-g');
        htric = plot([tree1_xyz(t1(3),1) tree1_xyz(t1(1),1)],...
            [tree1_xyz(t1(3),2) tree1_xyz(t1(1),2)],'-g');
        htrid = plot([tree2_xyz(t2(1),1) tree2_xyz(t2(2),1)] + ti(1),...
            [tree2_xyz(t2(1),2) tree2_xyz(t2(2),2)] + ti(2),'-r');
        htrie = plot([tree2_xyz(t2(2),1) tree2_xyz(t2(3),1)] + ti(1),...
            [tree2_xyz(t2(2),2) tree2_xyz(t2(3),2)] + ti(2),'-r');
        htrif = plot([tree2_xyz(t2(3),1) tree2_xyz(t2(1),1)] + ti(1),...
            [tree2_xyz(t2(3),2) tree2_xyz(t2(1),2)] + ti(2),'-r');
        foo = 1;
        set(hcomba, 'visible', 'off');
        set(hcombb, 'visible', 'off');
        set(hcombc, 'visible', 'off');
        set(htria, 'visible', 'off');
        set(htrib, 'visible', 'off');
        set(htric, 'visible', 'off');
        set(htrid, 'visible', 'off');
        set(htrie, 'visible', 'off');
        set(htrif, 'visible', 'off');
        %}
        A = tree1_xyz(t1,:);
        B = tree2_xyz(t2,:);
        %[Rhat,that] = rigid_transform_3D(A,B);
        %tree1_xyzt = (Rhat*tree1_xyz') + repmat(that,1,n_tree1);
        %tree1_xyzt = tree1_xyzt';
        % Transform Tree2 into coordinate system defined by Tree1
        [Rhat,that] = rigid_transform_3D(B,A);
        tree2_xyzt = (Rhat*tree2_xyz') + repmat(that,1,n_tree2);
        tree2_xyzt = tree2_xyzt';
        
        match1 = repmat(tree1_xyz,[1,1,n_tree2]);
        match2 = zeros(n_tree1,3,n_tree2);
        for t = 1:n_tree2;
            match2(:,:,t) = repmat(tree2_xyzt(t,:),[n_tree1,1]);
        end
        error_xyz = squeeze(sum((match1-match2).^2,2));
        [m1,m2] = find(error_xyz<t_xyz);
        n_match = numel(m1);
        
        if n_match > n_match_best;
            n_match_best = n_match;
            % R_best = Rhat;
            % t_best = that;
            m1_best = m1;
            m2_best = m2;
        end
        
        %{
        ht = figure('position', [587 95 1026 832]);
        hold on
        h1 = plot3(0,0,30,'ok','markerfacecolor',[.57 1 .53]);
        h2 = scatter(0,0,30,'ok','markerfacecolor',[.96 .46 .50]);
        set(h1, 'visible','off');
        set(h2, 'visible','off');
        for s = 1:n_tree1;
            h = filledCircle([tree1_xyz(s,1) tree1_xyz(s,2)],tree1_r(s),1000,'g');
            set(h,'FaceAlpha',.5);
        end
        for s = 1:n_tree2;
            h = filledCircle([tree2_xyzt(s,1) tree2_xyzt(s,2)],tree2_r(s),1000,'r');
            set(h,'FaceAlpha',.5);
        end
        for m = 1:n_match;
            h = filledCircle([tree1_xyz(m1(m),1) tree1_xyz(m1(m),2)],tree1_r(m1(m)),1000,'b');
            set(h,'FaceAlpha',.5);
        end
        for m = 1:n_match;
            h = filledCircle([tree2_xyzt(m2(m),1) tree2_xyzt(m2(m),2)],tree2_r(m2(m)),1000,'c');
            set(h,'FaceAlpha',.5);
        end
        hcomba = plot([tree1_xyz(t1(1),1) tree2_xyzt(t2(1),1)],...
            [tree1_xyz(t1(1),2) tree2_xyzt(t2(1),2)],'-k');
        hcombb = plot([tree1_xyz(t1(2),1) tree2_xyzt(t2(2),1)],...
            [tree1_xyz(t1(2),2) tree2_xyzt(t2(2),2)],'-k');
        hcombc = plot([tree1_xyz(t1(3),1) tree2_xyzt(t2(3),1)],...
            [tree1_xyz(t1(3),2) tree2_xyzt(t2(3),2)],'-k');
        htria = plot([tree1_xyz(t1(1),1) tree1_xyz(t1(2),1)],...
            [tree1_xyz(t1(1),2) tree1_xyz(t1(2),2)],'-g');
        htrib = plot([tree1_xyz(t1(2),1) tree1_xyz(t1(3),1)],...
            [tree1_xyz(t1(2),2) tree1_xyz(t1(3),2)],'-g');
        htric = plot([tree1_xyz(t1(3),1) tree1_xyz(t1(1),1)],...
            [tree1_xyz(t1(3),2) tree1_xyz(t1(1),2)],'-g');
        htrid = plot([tree2_xyzt(t2(1),1) tree2_xyzt(t2(2),1)] ,...
            [tree2_xyzt(t2(1),2) tree2_xyzt(t2(2),2)] ,'-r');
        htrie = plot([tree2_xyzt(t2(2),1) tree2_xyzt(t2(3),1)] ,...
            [tree2_xyzt(t2(2),2) tree2_xyzt(t2(3),2)] ,'-r');
        htrif = plot([tree2_xyzt(t2(3),1) tree2_xyzt(t2(1),1)] ,...
            [tree2_xyzt(t2(3),2) tree2_xyzt(t2(1),2)] ,'-r');
        h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
        axis equal;
        axis([-16 16 -16 16]);
        xlabel('Easting');
        ylabel('Northing');
        grid on
        foo = 1;
        %}
        %{
        delete(ht);
        %}
    end
end

foo = 1;
if isempty(m1_best);
    Rhat = [];
    that = [];
    return
end

%% Re-estimation
%
A = tree1_xyz(m1_best,:);
B = tree2_xyz(m2_best,:);
if size(A,1) < 3;
    return
end
%[Rhat,that] = rigid_transform_3D(A,B);
%tree1_xyzt = (Rhat*tree1_xyz') + repmat(that,1,n_tree1);
%tree1_xyzt = tree1_xyzt';
[Rhat,that] = rigid_transform_3D(B,A);
foo1 = (Rhat*tree2_xyz');
foo2 = repmat(that,1,n_tree2);
if any(size(foo1) ~= size(foo2));
    foo = 1
end
tree2_xyzt = foo1 + foo2;
tree2_xyzt = tree2_xyzt';


match1 = repmat(tree1_xyz,[1,1,n_tree2]);
match2 = zeros(n_tree1,3,n_tree2);
for t = 1:n_tree2;
    match2(:,:,t) = repmat(tree2_xyzt(t,:),[n_tree1,1]);
end
error_xyz = squeeze(sum((match1-match2).^2,2));
[m1,m2] = find(error_xyz<t_xyz);
n_match = numel(m1);

if n_match > n_match_best;
    n_match_best = n_match;
    % R_best = Rhat;
    % t_best = that;
    m1_best = m1;
    m2_best = m2;
end


%{
ht = figure('position', [587 95 1026 832]);
hold on
h1 = scatter(0,0,30,[.57 1 .53],'filled');
h2 = scatter(0,0,30,[.96 .46 .50],'filled');
set(h1, 'visible','off');
set(h2, 'visible','off');
for s = 1:n_tree1;
    h = filledCircle([tree1_xyz(s,1) tree1_xyz(s,2)],tree1_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_tree2;
    h = filledCircle([tree2_xyzt(s,1) tree2_xyzt(s,2)],tree2_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
for m = 1:n_match_best
    plot([tree1_xyz(m1_best(m),1) tree2_xyzt(m2_best(m),1)],...
        [tree1_xyz(m1_best(m),2) tree2_xyzt(m2_best(m),2)],'-k');
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on
foo = 1;
%}

%{
ht = figure('position', [587 95 1026 832]);
hold on
h1 = scatter(0,0,30,[.57 1 .53],'filled');
h2 = scatter(0,0,30,[.96 .46 .50],'filled');
set(h1, 'visible','off');
set(h2, 'visible','off');
for s = 1:n_tree1;
    h = filledCircle([tree1_xyz(s,1) tree1_xyz(s,2)],tree1_r(s),1000,'g');
    set(h,'FaceAlpha',.5);
end
for s = 1:n_tree2;
    h = filledCircle([tree2_xyz(s,1) tree2_xyz(s,2)],tree2_r(s),1000,'r');
    set(h,'FaceAlpha',.5);
end
for m = 1:n_match_best
    plot([tree1_xyz(m1_best(m),1) tree2_xyz(m2_best(m),1)],...
        [tree1_xyz(m1_best(m),2) tree2_xyz(m2_best(m),2)],'-k');
end
h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
axis equal;
axis([-16 16 -16 16]);
xlabel('Easting');
ylabel('Northing');
grid on
foo = 1;
%}
end

