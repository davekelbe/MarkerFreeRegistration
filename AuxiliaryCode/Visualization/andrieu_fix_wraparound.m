function [seg_col, seg_row, seg_iter] = andrieu_fix_wraparound(...
    I12r,axis_a,axis_e,seg_iter,...
    seg_xyz_upper_left, seg_xyz_lower_left, seg_xyz_upper_right, seg_xyz_lower_right)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

axis_c = 1:numel(axis_a);
axis_r = 1:numel(axis_e);

% Convert xyz to row,col
% Left side
if numel(seg_xyz_upper_left)==3;
    [seg_row_upper_left, seg_col_upper_left] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_upper_left(1), seg_xyz_upper_left(2),seg_xyz_upper_left(3));
[seg_row_lower_left, seg_col_lower_left] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_lower_left(1), seg_xyz_lower_left(2),seg_xyz_lower_left(3));
% Right side
[seg_row_upper_right, seg_col_upper_right] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_upper_right(1), seg_xyz_upper_right(2),seg_xyz_upper_right(3));
[seg_row_lower_right, seg_col_lower_right] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_lower_right(1), seg_xyz_lower_right(2),seg_xyz_lower_right(3));
n_seg = 1;
else
[seg_row_upper_left, seg_col_upper_left] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_upper_left(:,1), seg_xyz_upper_left(:,2),seg_xyz_upper_left(:,3));
[seg_row_lower_left, seg_col_lower_left] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_lower_left(:,1), seg_xyz_lower_left(:,2),seg_xyz_lower_left(:,3));
% Right side
[seg_row_upper_right, seg_col_upper_right] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_upper_right(:,1), seg_xyz_upper_right(:,2),seg_xyz_upper_right(:,3));
[seg_row_lower_right, seg_col_lower_right] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_xyz_lower_right(:,1), seg_xyz_lower_right(:,2),seg_xyz_lower_right(:,3));
n_seg = numel(seg_row_upper_left);
end


seg_x = cell(n_seg,1);
seg_y = cell(n_seg,1);

[n_row, n_col,~] = size(I12r);
for s = 1:n_seg
    is_intersect = false;
    x = [seg_col_upper_left(s) seg_col_lower_left(s) seg_col_lower_right(s) seg_col_upper_right(s)];
    y = [seg_row_upper_left(s) seg_row_lower_left(s),seg_row_lower_right(s) seg_row_upper_right(s)];
    % If the difference in column values is small, continue 
    if (max(x) - min(x))< n_col/2;
        seg_x{s} = x;
        seg_y{s} = y;
        continue
    end
    % Add n_col to wraparound
    x(x<n_col/2)= x(x<n_col/2) + n_col;
    % Define right edge of image
    e2 = [n_col n_col; 0 n_row];
    % Find intersection points of segment polyon with edge
    if numel(x) ~= numel(y) 
        seg_x{s} = [];
        seg_y{s} = [];
        continue
    end
    if ~isequal(isnan(x), isnan(y))
        seg_x{s} = [];
        seg_y{s} = [];
        continue
    end
    %{
    [xi,yi,ii] = polyxpoly([x x(1)],[y y(1)],e2(1,:),e2(2,:));
    %{
       figure;
       hold on
       plot(e2(1,:),e2(2,:),'-r');
       plot([x x(1)],[y y(1)],'-g')
       %scatter(xi,yi,'b');
    %}
    % Redefine new polygon vertices if intersection
    if s ==231 || s ==236 ;
        foo = 1;
    end
    ii = sort(ii(:,1));
    if ii(1) == 1 && ii(2) == 2;
        is_intersect = true; %Adjusted
        xn1 = [xi(1) xi(2) x(3) x(4) x(1)];
        yn1 = [yi(1) yi(2) y(3) y(4) y(1)];
        xn2 = [xi(1) x(2) xi(2)]-n_col;
        yn2 = [yi(1) y(2) yi(2)];
        %{
        figure
        hold on
        plot(xn1,yn1,'b')
        plot(xn2+n_col,yn2,'r')
        foo = 1;
        %}
    end
    if ii(1) == 2 && ii(2) == 4;
        is_intersect = true;
        xn1 = [xi(2) xi(1) x(3) x(4)]; %Adjusted
        yn1 = [yi(2) yi(1) y(3) y(4)];
        xn2 = [xi(2) x(1) x(2) xi(1)]-n_col;
        yn2 = [yi(2) y(1) y(2) yi(1)];
        %{
        figure
        hold on
        plot(xn1,yn1,'b')
        plot(xn2+n_col,yn2,'r')
        foo = 1;
        %}
    end
    if ii(1) == 2 && ii(2) == 3;
        is_intersect =true; % Adjusted 
        xn1 = [xi(2) xi(1) x(3)]; 
        yn1 = [yi(2) yi(1) y(3)];
        xn2 = [xi(2) x(4) x(1) x(2) xi(1)]-n_col;
        yn2 = [yi(2) y(4) y(1) y(2) yi(1)];
        %{
        figure
        hold on
        plot(xn1,yn1,'b')
        plot(xn2+n_col,yn2,'r')
        foo = 1;
        %}
    end
    if ii(1) == 3 && ii(2) == 4;
        is_intersect = true;
        xn1 = [xi(2) x(2) xi(1)]-n_col;
        yn1 = [yi(2) y(2) yi(1)];
        xn2 = [xi(2) xi(1) x(3) x(4) x(1)];
        yn2 = [yi(2) yi(1) y(3) y(4) y(1)];
        %{
        figure
        hold on
        plot(xn1,yn1,'b')
        plot(xn2+n_col,yn2,'r')
        foo = 1;
        %}
    end
    if ii(1) == 1 && ii(2) == 4;
        is_intersect = true; %Adjusted 
        xn1 = [xi(2) xi(1) x(2) x(3) x(4)];
        yn1 = [yi(2) yi(1) y(2) y(3) y(4)];
        xn2 = [xi(2) x(1) xi(1)]-n_col;
        yn2 = [yi(2) y(1) yi(1)];
        %{
        figure
        hold on
        plot(xn1,yn1,'b')
        plot(xn2+n_col,yn2,'r')
        foo = 1;
        %}
    end
    %{
       figure; imagesc(mask_inliers);
       hold on
       plot(x,y,'-r','linewidth',1)
       plot(xn1,yn1,'-r','linewidth',1)
       plot(xn2,yn2,'-r','linewidth',1)
    %}
    if is_intersect
        seg_x{s} = xn1;
        seg_y{s} = yn1;
        n_seg = numel(seg_x);
        seg_x{n_seg+1}= xn2;
        seg_y{n_seg+1} = yn2;
        seg_iter(n_seg+1) = s;
    else
        seg_x{s} = x;
        seg_y{s} = y;
    end
    %}
    % No mapping toolbox 
    seg_x{s} = x;
    seg_y{s} = y;
end

is_valid = ~cellfun(@isempty,seg_x) & ~cellfun(@isempty,seg_y);
seg_iter = seg_iter(is_valid);
seg_col = seg_x(is_valid);
seg_row = seg_y(is_valid); 

%{
    figure;
    imagesc(I12r);
    hold on
    for s = 1:numel(seg_col); 
        plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-r')
    end
%}


