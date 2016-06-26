function [ I12ieq_all ] = andrieu_colorLUT_axes( data_x,data_y,data_z,I12ieq,colorLUT,...
    seg_index,axis_a,axis_e,linecolor,title,cbar_minmax,...
    seg_z_center,seg_y_center,seg_z_left, seg_y_left, seg_z_right, seg_y_right)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

I12ieq_all = I12ieq;
axis_c = 1:numel(axis_a);
axis_r = 1:numel(axis_e);

%points = [seg_z_center seg_y_center seg_z_left seg_y_left seg_z_right seg_y_right]

if nargin> 9 % Center
[seg_row1, seg_col1] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_z_center(:,1), seg_z_center(:,2),seg_z_center(:,3));
[seg_row2, seg_col2] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_y_center(:,1), seg_y_center(:,2),seg_y_center(:,3));
is_valid1 = (abs(seg_col1-seg_col2)<size(I12ieq,2)/3);
end
if nargin > 11 % Left side 
[seg_row3, seg_col3] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_z_left(:,1), seg_z_left(:,2),seg_z_left(:,3));
[seg_row4, seg_col4] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_y_left(:,1), seg_y_left(:,2),seg_y_left(:,3));
is_valid2 = (abs(seg_col3-seg_col4)<size(I12ieq,2)/3);
end
if nargin > 13 % Right side 
[seg_row5, seg_col5] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_z_right(:,1), seg_z_right(:,2),seg_z_right(:,3));
[seg_row6, seg_col6] = andrieu_xyz2rowcol(axis_a,axis_e, axis_c, axis_r, ...
    seg_y_right(:,1), seg_y_right(:,2),seg_y_right(:,3));
is_valid3 = (abs(seg_col5-seg_col6)<size(I12ieq,2)/3);
end

if nargin>13
is_valid = (is_valid1 & is_valid2 & is_valid3);
seg_row1 = seg_row1(is_valid);
seg_col1 = seg_col1(is_valid);
seg_row2 = seg_row2(is_valid);
seg_col2 = seg_col2(is_valid);
seg_row3 = seg_row3(is_valid);
seg_col3 = seg_col3(is_valid);
seg_row4 = seg_row4(is_valid);
seg_col4 = seg_col4(is_valid);
seg_row5 = seg_row5(is_valid);
seg_col5 = seg_col5(is_valid);
seg_row6 = seg_row6(is_valid);
seg_col6 = seg_col6(is_valid);
colorLUT  = colorLUT(is_valid);
seg_index = seg_index(is_valid);
end

if  nargin>13;
    n_sub = numel(seg_row1);
else
    n_sub = size(seg_z_center,1);
end
cmap = jet(100);
%color = zeros(n_sub,3);
if ~isempty(cbar_minmax)
    colorLUT(colorLUT<cbar_minmax(1)) = cbar_minmax(1);
    colorLUT(colorLUT>cbar_minmax(2)) = cbar_minmax(2);
end

colorLUT_norm = colorLUT - sign(min(colorLUT))*abs((min(colorLUT)));
colorLUT_norm = colorLUT_norm/ max(colorLUT_norm);
colorLUT_norm = round(colorLUT_norm*100);
colorLUT_norm(colorLUT_norm==0) = 1;
    
figure('visible','off');
%{
for s = 1:n_sub
    color(s,:) = 255*cmap(colorLUT_norm(s),:);
end
color = uint8(color);
%}
color = uint8(255*cmap(colorLUT_norm,:));
%{
figure;
x = round(linspace(min(colorLUT_norm),max(colorLUT_norm),50));
y = histc(colorLUT_norm,x);
for s = 1:50;
    hold on
    scatter(x(s),y(s),10,cmap(x(s),:),'filled');
end
%}
    %{
figure;
subplot(2,1,1)
scatter(1:n_sub,ones(n_sub,1),50,cmap(colorLUT_norm,:),'filled'); colorbar
subplot(2,1,2)
scatter(1:n_sub,ones(n_sub,1),50,color,'filled'); colorbar
%}

for s = 1:n_sub
    sub_x = data_x(seg_index{s});
    sub_y = data_y(seg_index{s});
    sub_z = data_z(seg_index{s});
    
    [sub_a, sub_e,~] = cart2sph(sub_x,sub_y,sub_z);
    sub_a = rad2deg(sub_a);
    sub_e = rad2deg(sub_e);
    is_a_neg = sub_a<0;
    sub_a(is_a_neg)=sub_a(is_a_neg)+360;
    sub_row = round(interp1(axis_e,axis_r,sub_e));
    sub_col = round(interp1(axis_a,axis_c,sub_a));
    is_remove = isnan(sub_row) | isnan(sub_col);
    sub_row = sub_row(~is_remove);
    sub_col = sub_col(~is_remove);
    for p = 1:numel(sub_row);
        %fprintf('pixel %g\n',p);
        val = I12ieq(sub_row(p),sub_col(p),3);
        I12ieq_all(sub_row(p),sub_col(p),1) = color(s,1);
        I12ieq_all(sub_row(p),sub_col(p),2) = color(s,2);
        I12ieq_all(sub_row(p),sub_col(p),3) = color(s,3);
    end
end

I12ieq_all(:,:,1) = fliplr(I12ieq_all(:,:,1));
I12ieq_all(:,:,2) = fliplr(I12ieq_all(:,:,2));
I12ieq_all(:,:,3) = fliplr(I12ieq_all(:,:,3));
%
figure;
imshow(I12ieq_all)
colormap('gray');
axis off; axis image
hold on
set(gca,'Units','normalized','Position',[0 0 1 1]);
[n_row, n_col,~] = size(I12ieq);
set(gcf,'Units','pixels','Position',[200 200 n_col n_row]);
for s = 1:n_sub
    plot([n_col-seg_col1(s),n_col-seg_col2(s)],[seg_row1(s) seg_row2(s)],linecolor,'linewidth', .75);
    plot([n_col-seg_col3(s),n_col-seg_col4(s)],[seg_row3(s) seg_row4(s)],linecolor,'linewidth', .75);
    plot([n_col-seg_col5(s),n_col-seg_col6(s)],[seg_row5(s) seg_row6(s)],linecolor,'linewidth', .75);
end

% Colorbar
y = [min(colorLUT) max(colorLUT)];
x = [0 100];
xq = 0:10:100;
yq = interp1(x,y,xq);
colormap(cmap);
hcolorbar = colorbar('peer',gca);
set(hcolorbar,'location','south');
set(hcolorbar,'xTick',0:10:100);
set(hcolorbar,'xTickLabel',sprintf('%2.2f |',yq));
set(get(hcolorbar,'xlabel'),'String',title);

%}
end

