function [ match ] = match_stem( unique_SEG_x, unique_SEG_y, unique_SEG_r, unique_SEG_id,...
                             unique_ROI_x, unique_ROI_y, unique_ROI_r,unique_ROI_id,...
                             unique_NEON_x, unique_NEON_y, unique_NEON_r, ...
                             delta_xy )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%%
n_unique_SEG = numel(unique_SEG_x);
n_unique_ROI = numel(unique_ROI_x);
n_unique_NEON = numel(unique_NEON_x);

unique_SEG_xy = sqrt(unique_SEG_x.^2 + unique_SEG_y.^2);
unique_ROI_xy = sqrt(unique_ROI_x.^2 + unique_ROI_y.^2);
unique_NEON_xy = sqrt(unique_NEON_x.^2 + unique_NEON_y.^2);


match_ROI_ix = nan(n_unique_SEG,1);
match_NEON_ix = nan(n_unique_SEG,1);

for s = 1:n_unique_SEG;
    is_in_ROI = (sqrt((unique_ROI_x-unique_SEG_x(s)).^2 + (unique_ROI_y-unique_SEG_y(s)).^2)<=delta_xy);
    is_in_NEON = (sqrt((unique_NEON_x-unique_SEG_x(s)).^2 + (unique_NEON_y-unique_SEG_y(s)).^2)<=delta_xy);
    if sum(is_in_ROI)>0
        r_diff = abs(unique_ROI_r-unique_SEG_r(s));
        r_diff(~is_in_ROI) = inf;
        [~,imin] = min(r_diff);
        if (r_diff(imin)./unique_ROI_r(imin) < 1);
            match_ROI_ix(s) = imin;
        end
    end
    if sum(is_in_NEON)>0
        r_diff = abs(unique_NEON_r-unique_SEG_r(s));
        r_diff(~is_in_NEON) = inf;
        [~,imin] = min(r_diff);
        if (r_diff(imin)./unique_NEON_r(imin) < 1);
            match_NEON_ix(s) = imin;
        end
    end
end


isnotnan_ROI = ~isnan(match_ROI_ix);
isnotnan_NEON = ~isnan(match_NEON_ix);
match_ROI_roi_ix = match_ROI_ix(isnotnan_ROI);  
match_NEON_neon_ix = match_NEON_ix(isnotnan_NEON);
clear match_NEON_ix match_ROI_ix

unique_NEON_id = (1:n_unique_NEON)';

%match_ROI_index = (1:numel(match_ROI_roi_ix))';
match_ROI_roi_id = unique_ROI_id(match_ROI_roi_ix); % ID of roi's for ROI match
match_ROI_roi_r= unique_ROI_r(match_ROI_roi_ix);
match_ROI_roi_xy = unique_ROI_xy(match_ROI_roi_ix);
match_ROI_roi_x = unique_ROI_x(match_ROI_roi_ix);
match_ROI_roi_y = unique_ROI_y(match_ROI_roi_ix);
match_ROI_seg_id = unique_SEG_id(isnotnan_ROI);
match_ROI_seg_r = unique_SEG_r(isnotnan_ROI);
match_ROI_seg_xy = unique_SEG_xy(isnotnan_ROI);
match_ROI_seg_x = unique_SEG_x(isnotnan_ROI);
match_ROI_seg_y = unique_SEG_y(isnotnan_ROI);


%match_NEON_index = (1:numel(match_NEON_neon_ix))';
match_NEON_neon_id = unique_NEON_id(match_NEON_neon_ix); % ID of NEON for neon match
match_NEON_neon_r= unique_NEON_r(match_NEON_neon_ix);
match_NEON_neon_xy = unique_NEON_xy(match_NEON_neon_ix);
match_NEON_neon_x = unique_NEON_x(match_NEON_neon_ix);
match_NEON_neon_y = unique_NEON_y(match_NEON_neon_ix);
match_NEON_seg_id = unique_SEG_id(isnotnan_NEON);
match_NEON_seg_r = unique_SEG_r(isnotnan_NEON);
match_NEON_seg_xy = unique_SEG_xy(isnotnan_NEON);
match_NEON_seg_x = unique_SEG_x(isnotnan_NEON);
match_NEON_seg_y = unique_SEG_y(isnotnan_NEON);

% Check for duplicates 

[n, bin] = histc(match_ROI_roi_id,unique(match_ROI_roi_id) ); % Repeated seg
multiple = find(n > 1);
ROI_duplicate = find(ismember(bin, multiple));
%duplicate_index = match_ROI_index(ROI_duplicate);
duplicate_ROI_roi_ix = match_ROI_roi_ix(ROI_duplicate);
duplicate_ROI_roi_r = match_ROI_roi_r(ROI_duplicate);
duplicate_ROI_roi_id = match_ROI_roi_id(ROI_duplicate);
duplicate_ROI_seg_r = match_ROI_seg_r(ROI_duplicate);
duplicate_ROI_seg_id = match_ROI_seg_id(ROI_duplicate);
duplicate_ROI_seg_rdiff = abs(duplicate_ROI_roi_r-duplicate_ROI_seg_r);
unique_duplicate_ROI_id = unique(duplicate_ROI_roi_id);
n_unique_ROI = numel(unique_duplicate_ROI_id);
n_duplicate = numel(duplicate_ROI_roi_id);
is_keep = false(n_duplicate,1);
for u = 1:n_unique_ROI;
    is_u = duplicate_ROI_roi_id==unique_duplicate_ROI_id(u);
    temp = duplicate_ROI_seg_rdiff;
    temp(~is_u) = inf;
    [~,imin] = min(temp);
    is_keep(imin) = true;
end
index_keep = true(numel(match_ROI_roi_ix),1);
index_keep(ROI_duplicate(~is_keep)) = false;
match.ROI_roi_id = match_ROI_roi_id(index_keep);
match.ROI_roi_r = match_ROI_roi_r(index_keep);
match.ROI_roi_xy = match_ROI_roi_xy(index_keep);
match.ROI_roi_x = match_ROI_roi_x(index_keep);
match.ROI_roi_y = match_ROI_roi_y(index_keep);
match.ROI_seg_id = match_ROI_seg_id(index_keep);
match.ROI_seg_r = match_ROI_seg_r(index_keep);
match.ROI_seg_xy = match_ROI_seg_xy(index_keep);
match.ROI_seg_x = match_ROI_seg_x(index_keep);
match.ROI_seg_y = match_ROI_seg_y(index_keep);

[n, bin] = histc(match_NEON_neon_id,unique(match_NEON_neon_id) ); % Repeated seg
multiple = find(n > 1);
NEON_duplicate = find(ismember(bin, multiple));
%duplicate_index = match_ROI_index(ROI_duplicate);
duplicate_NEON_neon_ix = match_NEON_neon_ix(NEON_duplicate);
duplicate_NEON_neon_r = match_NEON_neon_r(NEON_duplicate);
duplicate_NEON_neon_id = match_NEON_neon_id(NEON_duplicate);
duplicate_NEON_seg_r = match_NEON_seg_r(NEON_duplicate);
duplicate_NEON_seg_id = match_NEON_seg_id(NEON_duplicate);
duplicate_NEON_seg_rdiff = abs(duplicate_NEON_neon_r-duplicate_NEON_seg_r);
unique_duplicate_NEON_id = unique(duplicate_NEON_neon_id);
n_unique_NEON = numel(unique_duplicate_NEON_id);
n_duplicate = numel(duplicate_NEON_neon_id);
is_keep = false(n_duplicate,1);
for u = 1:n_unique_NEON;
    is_u = duplicate_NEON_neon_id==unique_duplicate_NEON_id(u);
    temp = duplicate_NEON_seg_rdiff;
    temp(~is_u) = inf;
    [~,imin] = min(temp);
    is_keep(imin) = true;
end
index_keep = true(numel(match_NEON_neon_ix),1);
index_keep(NEON_duplicate(~is_keep)) = false;
match.NEON_neon_id = match_NEON_neon_id(index_keep);
match.NEON_neon_r = match_NEON_neon_r(index_keep);
match.NEON_neon_xy = match_NEON_neon_xy(index_keep);
match.NEON_neon_x = match_NEON_neon_x(index_keep);
match.NEON_neon_y = match_NEON_neon_y(index_keep);
match.NEON_seg_id = match_NEON_seg_id(index_keep);
match.NEON_seg_r = match_NEON_seg_r(index_keep);
match.NEON_seg_xy = match_NEON_seg_xy(index_keep);
match.NEON_seg_x = match_NEON_seg_x(index_keep);
match.NEON_seg_y = match_NEON_seg_y(index_keep);

end

