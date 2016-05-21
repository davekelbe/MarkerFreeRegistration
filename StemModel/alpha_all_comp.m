function [ x_bnd, y_bnd, is_in_best, x_bnd_all, y_bnd_all, is_in_all ] = alpha_all_comp( bnd, x,y )
%ALPHA_MAXCOMP Determine the maximum alpha component by # inclusive points 
%   
%   Inputs:
%       bnd     - boundary indexes from alphavol()
%       x       - x values of all points from alphavol()
%       y       - y values of all points from alphavol()
%
%   Outptus:
%       x_bnd   - x values on the boundary
%       y_bnd   - y values on the boundary
%       is_in_best - inliers to points [x,y] which are inside largest shape
%
%   Example:
%   Random:
%   x = [rand(10,1) + 3;  rand(30,1)];
%   y = [rand(10,1) + 3;  rand(30,1)];
%   n_pts = numel(x);
%   [~,S] = alphavol([x,y],1,1);
%   [ x_bnd, y_bnd, is_in_best ] = alpha_maxcomp( S.bnd, x,y );
%   figure;
%   scatter(x(is_in_best),y(is_in_best),'b','filled')
%   hold on
%   scatter(x(~is_in_best),y(~is_in_best),'r','filled')
%   plot(x_bnd,y_bnd,'-r','LineWidth',2)
%   
%

n_b = size(bnd,1);
n_pts = numel(x);
g = 1;
p = 1;
n = 1;
start_new = false;
is_used = false(n_b,1);
vertex{g} = zeros(n_b,1);
i = 1;
notdone = true;
still_open = false;
counter = 0;


%hold on
%h = scatter(u(1),v(1));
while notdone
    if start_new;
        g = g+1; % Group number
        vertex{g} = zeros(n_b,1);
        i = 1; % Index in group
        start_new = false;
        ix_notused = find(~is_used);
        if isempty(ix_notused);
            notdone = false;
            continue;
        end
        p = ix_notused(1); % Index in points array
    end
    %   if ishandle(h);
    %       set(h,'visible', 'off');
    %   end
    %   h = scatter(x(bnd(p,1)),y(bnd(p,1)),70,'g','filled');
    if p>numel(is_used)
        x_bnd = [];
        y_bnd = [];
        is_in_best = [];
        return
    end
    if is_used(p);
        %       scatter(x(bnd(p,1)),y(bnd(p,1)),90,'r')
        p = p + 1;
        continue
    end
    if sum(is_used) == numel(is_used);
        notdone = false;
    end
    if numel(p) ~= 1;
        p = p(1);
        is_used(p) = true;
    else
        is_used(p) = true;
    end
    if i ==1
        is_used(p) = false;
        %       scatter(x(bnd(p,1)),y(bnd(p,1)),'b','filled');
    end
    %   scatter(x(bnd(p,1)),y(bnd(p,1)),'b','filled')
    vertex{g}(i) = bnd(p,1);
    %   if i > 1
    %       plot(x(vertex{g}(i-1:i)),y(vertex{g}(i-1:i)),'-g','linewidth',2);
    %   end
    if vertex{g}(1) == vertex{g}(i) && i~=1;
        start_new = true;
    end
    counter = counter + 1;
    
    p = find(bnd(:,1)==bnd(p,2)); %Update p
    if numel(p) > 1;
        ix_not_used = find(~is_used(p));
        if numel(ix_not_used) == 0;
            p = p(1);
        else
            p = p(ix_not_used(1));
        end
    end
    i = i + 1;
end
vertex = vertex(1:numel(vertex)-1);
n_g = size(vertex,2);
is_in = false(n_pts,n_g);
is_on = false(n_pts,n_g);
for g = 1:n_g;
    vertex{g} = vertex{g}((vertex{g})>0);
    [is_in(:,g) is_on(:,g)] = inpolygon(x,y,x(vertex{g}),y(vertex{g}));
    is_in(:,g) = is_in(:,g) | is_on(:,g);
end
n_in = sum(is_in,1);

[~, g_best] = max(n_in);    % Caveat that first is chosen if equal
bnd_best = [vertex{g_best}; vertex{g_best}(1)];
if n_g>1;
    is_in_best = is_in(:,g_best);
else
    is_in_best = is_in;
end
x_bnd = x(bnd_best);
y_bnd = y(bnd_best);

if nargout >3
bnd_all = cell(n_g,1);
x_bnd_all = cell(n_g,1);
y_bnd_all = cell(n_g,1);
for g = 1:n_g;
    bnd_all{g} = [vertex{g}; vertex{g}(1)];
    x_bnd_all{g} = x(bnd_all{g});
    y_bnd_all{g} = y(bnd_all{g});
end
is_in_all = is_in;
end


end

