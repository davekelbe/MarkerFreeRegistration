function [ I12ieq_c ] = color_andrieu( I12ieq,axis_a, axis_e, all_x, all_y, all_z, is_inlier_alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    axis_c = 1:numel(axis_a);
    axis_r = 1:numel(axis_e);
    [sub_a, sub_e,~] = cart2sph(all_x(is_inlier_alpha),all_y(is_inlier_alpha),all_z(is_inlier_alpha));
    sub_a = rad2deg(sub_a);
    sub_e = rad2deg(sub_e);
    is_a_neg = sub_a<0;
    sub_a(is_a_neg)=sub_a(is_a_neg)+360;
    sub_row= round(interp1(axis_e,axis_r,sub_e));
    sub_col = round(interp1(axis_a,axis_c,sub_a));
    I12ieq_c = I12ieq;
    is_remove = isnan(sub_row) | isnan(sub_col);
    sub_row = sub_row(~is_remove);
    sub_col = sub_col(~is_remove);
    for p = 1:numel(sub_row);
        val = I12ieq(sub_row(p),sub_col(p),3);
        I12ieq_c(sub_row(p),sub_col(p),1) = val/2;
        I12ieq_c(sub_row(p),sub_col(p),2) = val/2;
    end
    I12ieq_c(:,:,1) = fliplr(I12ieq_c(:,:,1));
    I12ieq_c(:,:,2) = fliplr(I12ieq_c(:,:,2));
    I12ieq_c(:,:,3) = fliplr(I12ieq_c(:,:,3));
end

