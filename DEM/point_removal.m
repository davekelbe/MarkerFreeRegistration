function [ x, y, z_new ] = point_removal3( data_x, data_y, data_z, minindex, dt, t_angle_deviation, pngPath, plot1, verbose, saveim)
%POINTREMOVAL Remove outlier points from delaunay triangulation
%   Detailed explanation goes here
minindex = logical(minindex);

x = data_x(minindex);
y = data_y(minindex);
z = data_z(minindex);
z_new = z;

n_bad_simplex_previous = inf;

while true
    
    z = z_new;
    
    tr = TriRep(dt.Triangulation,x,y,z); % Create 2D triangulation representation
    P = incenters(tr);
    fn = faceNormals(tr);
    
    % figure;
    % trimesh(dt.Triangulation,x,y,z)%, ...
    %     % 'FaceColor', 'cyan', 'faceAlpha', 0.8);
    % axis equal;
    % hold on;
    % quiver3(P(:,1),P(:,2),P(:,3), ...
    %      fn(:,1),fn(:,2),fn(:,3),0.5, 'color','r');
    % hold off;
    
    if verbose;
        figure;
        scatter3(fn(:,1),fn(:,2),fn(:,3),20,fn(:,3), 'filled');
        title('Unit Normal to Triangulation Planes');
        xlabel('x', 'fontsize', 14);
        ylabel('y', 'fontsize', 14);
        zlabel('z', 'fontsize', 14);
    end
    if saveim;  saveas( gcf(), [ pngPath 'normals' num2str( plot1 ),'.png' ], 'png' ); end
    
    %Find the average normal vector
    fn_ave = mean(fn,1);
    %Find the angle between all normal vectors and average;
    diff_angle = acosd(dot(fn, repmat(fn_ave,[size(fn,1),1]),2));
    diff_angle_mean = mean(diff_angle);
    index = (diff_angle<diff_angle_mean + t_angle_deviation)&...
        (diff_angle>diff_angle_mean - t_angle_deviation); % "Good ones"
    [theta, ~]  = cart2pol(P(:,1),P(:,2));
    theta = rad2deg(theta)+180;
    quadrant = floor(theta/90);
    xsign = sign(fn(:,1));
    ysign = sign(fn(:,2));
    Q0 = (quadrant == 0) & (xsign == 1) & (ysign == 1);
    Q1 = (quadrant == 1) & (xsign == -1) & (ysign == 1);
    Q2 = (quadrant == 2) & (xsign == -1) & (ysign == -1);
    Q3 = (quadrant == 3) & (xsign == 1) & (ysign == -1);
    
    peaks = (Q0 | Q1 | Q2 | Q3);
    index1 = ~(~index & peaks);
    xy = sqrt(P(:,1).^2 + P(:,2).^2);
    nearcenter = (xy<10);
    index2 = ~(~index & nearcenter);
    
    index = ~(~index1|~index2);
    n_bad_simplex = sum(~index);
    %fprintf('\n%5.0f bad simplexes remain\n',n_bad_simplex); 
    if sum(~index)<1 || n_bad_simplex>n_bad_simplex_previous ; % No more "bad ones"
        break
    end
    n_bad_simplex_previous = n_bad_simplex;
    
    if verbose;
        figure;
        trimesh(dt.Triangulation,x,y,z)
        axis equal;
        hold on;
        quiver3(P(index,1),P(index,2),P(index,3), ...
            fn(index,1),fn(index,2),fn(index,3),0.5, 'color','r');
        quiver3(P(~index,1),P(~index,2),P(~index,3), ...
            fn(~index,1),fn(~index,2),fn(~index,3),2, 'color','k');
       %  quiver3(P(i,1),P(newindex,2),P(newindex,3), ...
       %     fn(newindex,1),fn(newindex,2),fn(newindex,3),0.5, 'color','r');
       % quiver3(P(~newindex,1),P(~newindex,2),P(~newindex,3), ...
       %     fn(~newindex,1),fn(~newindex,2),fn(~newindex,3),2, 'color','k');
        hold off;
    end;
    if saveim;  saveas( gcf(), [ pngPath 'quiver' num2str( plot1 ),'.png' ], 'png' ); end
    
    bad_simplex_verts = dt.Triangulation(~index,:);
    simplex_z = [z(bad_simplex_verts(:,1))...
        z(bad_simplex_verts(:,2))...
        z(bad_simplex_verts(:,3))];
    [Ymin,~] = min(simplex_z,[],2);
    [~,Imax] = max(simplex_z,[],2);
    
    n_bad_simplex = numel(Ymin);
    %I_simplex_verts_min = zeros(n_bad_simplex,1);
    I_simplex_verts_max = zeros(n_bad_simplex,1);
    for i = 1:n_bad_simplex;
        %I_simplex_verts_min(i) = bad_simplex_verts(i,Imin(i));
        I_simplex_verts_max(i) = bad_simplex_verts(i,Imax(i));
    end
    z_new = z;
    z_new(I_simplex_verts_max) = Ymin;
    
end



end


