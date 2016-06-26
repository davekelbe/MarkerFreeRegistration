function [ seg_row, seg_col, seg_iter ] = andrieu_tree_edges(I12r, axis_a, axis_e, tree )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
%figure; imagesc(I12ieq);

% Find edges as viewed from sensor 
n_tree = numel(tree);
if n_tree == 1;
    n_tree = 0;
end
for t = 1:n_tree; 
    tree_n_s = numel(tree(t).r);
    temp_a = circshift(tree(t).loc,[0 -1]) - tree(t).loc;
    temp_a(:,end) = temp_a(:,end-1);
    temp_a = temp_a';
    edges_left = zeros(tree_n_s,3);
    edges_right = zeros(tree_n_s,3);
    for s = 1:tree_n_s;
    q1 = temp_a(s,:)';
    p2  = tree(t).loc(:,s);
        % Second vector is from sensor to object
        w2 = p2-dot(p2,q1)*q1;
        q2 = w2/norm(w2);
        % Third vector orthogonal
        %p3 = [0 1 1]';
        p3 = p2;
        p3(1) = -p3(1);
        p3(3) = -p3(3);
        w3 = p3 - dot(p3,q1)*q1 - dot(p3,q2)*q2;
        q3 = w3/norm(w3);
        if q3(1)<0;
            q3 = -q3;
        end
        M1 = [q1'; q2'; q3'];
        edges_left(s,:) = (tree(t).loc(:,s) - tree(t).r(s)*q3)';
        edges_right(s,:) = (tree(t).loc(:,s) + tree(t).r(s)*q3)';
    end
    tree(t).left = edges_left';
    tree(t).right = edges_right';
end

tree_n_s = zeros(n_tree,1);
for t = 1:n_tree;
    tree_n_s(t) = numel(tree(t).r);
end
tree_n_seg = tree_n_s - 1;

n_seg = sum(tree_n_seg);
seg_upper_left = zeros(n_seg,3);
seg_upper_right = zeros(n_seg,3);
seg_lower_left = zeros(n_seg,3);
seg_lower_right = zeros(n_seg,3);

ix_end = cumsum(tree_n_seg);
ix_start = ix_end-tree_n_seg+1;
for t = 1:n_tree;
    seg_lower_left(ix_start(t):ix_end(t),:) = tree(t).left(:,1:end-1)';
    seg_lower_right(ix_start(t):ix_end(t),:) = tree(t).right(:,1:end-1)';
    seg_upper_left(ix_start(t):ix_end(t),:) = tree(t).left(:,2:end)';
    seg_upper_right(ix_start(t):ix_end(t),:) = tree(t).right(:,2:end)';
end
seg_iter = (1:sum(tree_n_seg))';


[seg_col, seg_row, seg_iter] = andrieu_fix_wraparound(...
    I12r,axis_a,axis_e,seg_iter,...
    seg_upper_left, seg_lower_left, seg_upper_right, seg_lower_right);


is_finite = false(numel(seg_col),1);
for s = 1:numel(seg_col);
        is_finite(s) = ~(any(~isfinite(seg_col{s})) || any(~isfinite(seg_row{s})));
end
    seg_iter = seg_iter(is_finite);
    seg_row = seg_row(is_finite);
    seg_col = seg_col(is_finite);
    
    %{
    figure;
    imagesc(I12ieq);
    hold on
    for s = 1:numel(seg_col);
        plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-r')
    end
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    [n_row, n_col,~] = size(I12ieq);
    set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
%}


end

