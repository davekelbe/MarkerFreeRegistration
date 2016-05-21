function [ hfigrmse ] = label_plot_rmse( aux )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

info_site = aux.info_site;
i = aux.i;
match_R = aux.match_R;
match_t = aux.match_t;
match_rx = aux.match_rx;
match_ry = aux.match_ry;
match_rz = aux.match_rz;
match_tx = aux.match_tx;
match_ty = aux.match_ty;
match_tz = aux.match_tz;
matchret_R = aux.matchret_R;
matchret_t = aux.matchret_t;
matchret_rx = aux.matchret_rx;
matchret_ry = aux.matchret_ry;
matchret_rz = aux.matchret_rz;
matchret_tx = aux.matchret_tx;
matchret_ty = aux.matchret_ty;
matchret_tz = aux.matchret_tz;
match_tdiff = aux.match_tdiff;
match_txydiff = aux.match_txydiff; 
match_RMSE = aux.match_RMSE;
G_R = aux.G_R;
G_t = aux.G_t;
G_rx = aux.G_rx;
G_ry = aux.G_ry;
G_rz = aux.G_rz;
G_tx = aux.G_tx;
G_ty = aux.G_ty;
G_tz = aux.G_tz;
is_guess = aux.is_guess;
n_scan = aux.n_scan;
vec_plot = aux.vec_plot;
vec_x = aux.vec_x;
vec_y = aux.vec_y;
RT = aux.RT;
isempty0 = aux.isempty0;
i = aux.i;

n_S = size(match_R,1);

x_axis = [0 0 0; 1 0 0]';
y_axis = [0 0 0; 0 1 0]';
z_axis = [0 0 0; 0 0 1]';
[x_sph, y_sph, z_sph] = sphere;
exclude_site = [];

% For each i in a site
% Plot all #1
%err_alpha = match12_RMSE(i,:);
%err_color = (double(vec2cmap(err_alpha, 'jet', 0,2 )))./255;
%figure('position', [251         392        1156         559]);
hfigrmse = figure('position', [68         413        1593         538]);
subplot(1,3,1)
if RT==1;
    title(sprintf('t diff for Site %3.0f-%g', info_site,i));
elseif RT ==2;
    title(sprintf('t diff for Site %3.0f-%g', info_site,i));
end
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
err_color = vec2cmap(match_tdiff(i,:), 'jet');
% Plot all pairwise matches from j
for j = 1:n_S;
    if ~isempty0(i,j);
        x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
        y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
        z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
        x_axisrt = matchret_R{i,j}*x_axis + repmat(matchret_t{i,j},[1,2]);
        y_axisrt = matchret_R{i,j}*y_axis + repmat(matchret_t{i,j},[1,2]);
        z_axisrt = matchret_R{i,j}*z_axis + repmat(matchret_t{i,j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),.1+x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),.1+y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),.1+z_axist(3,:),'-b', 'linewidth',2)
        plot3(x_axisrt(1,:),x_axisrt(2,:),.1+x_axisrt(3,:),'-r', 'linewidth',2)
        plot3(y_axisrt(1,:),y_axisrt(2,:),.1+y_axisrt(3,:),'-g', 'linewidth',2)
        plot3(z_axisrt(1,:),z_axisrt(2,:),.1+z_axisrt(3,:),'-b', 'linewidth',2)
        plot3([x_axist(1,1) x_axisrt(1,1)],[x_axist(2,1) x_axisrt(2,1)],[x_axist(3,1) x_axisrt(3,1)],'-r')
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
        %hrmse = filledCircle(match_t{i,j}(1:2),1,1000,err_color(j,:));
         h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
             z_sph+match_t{i,j}(3));
         alpha(0.2)
         set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
    end
end
curr1_xlim = xlim;
curr1_xlim(1) = curr1_xlim(1) -2;
curr1_xlim(2) = curr1_xlim(2) +2;
curr1_ylim = ylim;
curr1_ylim(1) = curr1_ylim(1) -2;
curr1_ylim(2) = curr1_ylim(2) +2;
xrange = curr1_xlim(2) - curr1_xlim(1);
yrange = curr1_ylim(2) - curr1_ylim(1);
if xrange>yrange;
    diff = xrange-yrange;
    curr1_ylim(1) = curr1_ylim(1) - diff/2;
    curr1_ylim(2) = curr1_ylim(2) + diff/2;
elseif yrange>xrange;
    diff = yrange-xrange;
    curr1_xlim(1) = curr1_xlim(1) - diff/2;
    curr1_xlim(2) = curr1_xlim(2) + diff/2;
end

xlimMode = 'manual' ;
ylimMode = 'manual';
zlimMode = 'manual';

% Plot reference map
subplot(1,3,2)
if RT==1;
    title(sprintf('txy for Site %3.0f-%g', info_site,i));
elseif RT ==2;
    title(sprintf('txy for Site %3.0f-%g', info_site,i));
end
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
err_color = vec2cmap(match_txydiff(i,:), 'jet');
% Plot all pairwise matches from j
for j = 1:n_S;
    if ~isempty0(i,j);
        x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
        y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
        z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),.1+x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),.1+y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),.1+z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
        %hrmse = filledCircle(match_t{i,j}(1:2),1,1000,err_color(j,:));
         h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
             z_sph+match_t{i,j}(3));
         alpha(0.2)
         set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
    end
end


% Plot graph best guess
subplot(1,3,3)
if RT==1;
    title(sprintf('RMSE for Site %3.0f-%g', info_site,i));
elseif RT ==2;
    title(sprintf('RMSE for Site %3.0f-%g', info_site,i));
end
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
err_color = vec2cmap(match_RMSE(i,:), 'jet');
% Plot all pairwise matches from j
for j = 1:n_S;
    if ~isempty0(i,j);
        x_axist = match_R{i,j}*x_axis + repmat(match_t{i,j},[1,2]);
        y_axist = match_R{i,j}*y_axis + repmat(match_t{i,j},[1,2]);
        z_axist = match_R{i,j}*z_axis + repmat(match_t{i,j},[1,2]);
        plot3(x_axist(1,:),x_axist(2,:),.1+x_axist(3,:),'-r', 'linewidth',2)
        plot3(y_axist(1,:),y_axist(2,:),.1+y_axist(3,:),'-g', 'linewidth',2)
        plot3(z_axist(1,:),z_axist(2,:),.1+z_axist(3,:),'-b', 'linewidth',2)
        textloc = (x_axist(:,2) + y_axist(:,2))/2;
        textstr = sprintf('%g', j);
        text(textloc(1), textloc(2), textloc(3), textstr);
        %hrmse = filledCircle(match_t{i,j}(1:2),1,1000,err_color(j,:));
         h = surf(x_sph+match_t{i,j}(1), y_sph+match_t{i,j}(2), ...
             z_sph+match_t{i,j}(3));
         alpha(0.2)
         set(h, 'Facecolor',err_color(j,:)', 'edgecolor','none')
    end
end
curr3_xlim = xlim;
curr3_xlim(1) = curr3_xlim(1) -2;
curr3_xlim(2) = curr3_xlim(2) +2;
curr3_ylim = ylim;
curr3_ylim(1) = curr3_ylim(1) -2;
curr3_ylim(2) = curr3_ylim(2) +2;
xrange = curr3_xlim(2) - curr3_xlim(1);
yrange = curr3_ylim(2) - curr3_ylim(1);
if xrange>yrange;
    diff = xrange-yrange;
    curr3_ylim(1) = curr3_ylim(1) - diff/2;
    curr3_ylim(2) = curr3_ylim(2) + diff/2;
elseif yrange>xrange;
    diff = yrange-xrange;
    curr3_xlim(1) = curr3_xlim(1) - diff/2;
    curr3_xlim(2) = curr3_xlim(2) + diff/2;
end

% curr_xlim = [min([curr1_xlim(1),curr2_xlim(1),curr3_xlim(1)]) max([curr1_xlim(2),curr2_xlim(2),curr3_xlim(2)])];
% curr_ylim = [min([curr1_ylim(1),curr2_ylim(1),curr3_ylim(1)]) max([curr1_ylim(2),curr2_ylim(2),curr3_ylim(2)])];
% subplot(1,3,1);
% axis([curr_xlim curr_ylim]);
% subplot(1,3,2);
% axis([curr_xlim curr_ylim]);
% subplot(1,3,3);
% axis([curr_xlim curr_ylim]);


% for j = 1:n_S;
%     if isempty(match_t{i,j});
%         continue
%     end
%     if is_guess(i,j);
%         clr = 'blue';
%     else
%         clr = 'red';
%     end
%     subplot(1,3,1);
%     h5 = filledCircle(match_t{i,j}(1:2),1,1000,clr);
%     set(h5, 'FaceAlpha', 0.2)
%     subplot(1,3,2);
%     h6 = filledCircle([vec_x(j)-vec_x(i) vec_y(j)-vec_y(i)],1,1000,clr);
%     set(h6, 'FaceAlpha', 0.2)
%     subplot(1,3,3);
%     h3d = filledCircle(G_t{i,j}(1:2),1,1000,clr);
%     set(h3d, 'FaceAlpha', 0.2);
% end
    
end

