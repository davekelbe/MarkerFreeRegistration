function [ I12ieq_all ] = andrieu_colorLUT_axes3( data_x,data_y,data_z,I12ieq,colorLUT,...
    seg_index,axis_a,axis_e,title,cbar_minmax,...
    seg_col,seg_row,is_valid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

is_finite = isfinite(colorLUT);
if sum(~is_finite)>0;
    fprintf('\nWarning: some cylinders removed due to NaN\n')
end
colorLUT = colorLUT(is_finite);
seg_col = seg_col(is_finite);
seg_row = seg_row(is_finite);
seg_index = seg_index(is_finite);
is_valid = is_valid(is_finite);

I12ieq_all = I12ieq;
axis_c = 1:numel(axis_a);
axis_r = 1:numel(axis_e);
n_sub = numel(seg_col);

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
    
%figure('visible','off');
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

%I12ieq_all(:,:,1) = fliplr(I12ieq_all(:,:,1));
%I12ieq_all(:,:,2) = fliplr(I12ieq_all(:,:,2));
%I12ieq_all(:,:,3) = fliplr(I12ieq_all(:,:,3));
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
        if is_valid(s);
        plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-b')
        else
        plot([seg_col{s} seg_col{s}(1)],[seg_row{s} seg_row{s}(1)],'-r')
        end
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

