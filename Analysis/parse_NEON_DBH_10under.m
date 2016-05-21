function [ D ] = parse_NEON_DBH_10under( filepath_input,info_site )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


file.input = filepath_input;
fid  = fopen(file.input);
C = textscan(fid,'%s%u%s%f%f%f%f%s%s%s%s%s%s%s%u','HeaderLines',1, 'delimiter', ',');
%C = textscan(fid,'%u%u','HeaderLines',1, 'delimiter', ',');
C{end} = [C{end};0]; % Last entry missing
%is_empty = cellfun(@isempty,C{1});
n_i = numel(C{3});
is_site = false(n_i,1);
for i = 1:n_i;
    if strcmp(C{1}(i),sprintf('A%02.0f',info_site));
        is_site(i) = true;
    end
end

D.plot_id = C{1}(is_site);
D.module = C{2}(is_site);
D.nested_subplot = C{3}(is_site);
D.dbh = C{4}(is_site);
D.stem_height = C{5}(is_site);
D.ddh = C{6}(is_site);
D.canopy_diameter = C{7}(is_site);
D.vertical_position = C{8}(is_site);
D.common_species = C{9}(is_site);
D.other_species = C{10}(is_site);
D.stem_notes = C{11}(is_site);
D.stem_status = C{12}(is_site);
D.status_notes = C{13}(is_site);
D.shotgun_sampling = C{14}(is_site);
D.point_id = C{15}(is_site);

end
            