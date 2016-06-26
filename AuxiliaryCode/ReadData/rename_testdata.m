
dir_source = '/Volumes/Asterix/Test/Data/';
D = dir(dir_source);
dir_target = '/Volumes/Asterix/Test/DataRename/';

D = remove_hiddenfiles(D);
n_D = numel(D);
info_exp = 'test';
info_site = 0;
D_id = cell(n_D,1);

for d = 1:n_D;
    ix_dot = strfind(D{d},'.');
    ix_underscore = strfind(D{d},'_');
    id = D{d}(ix_underscore(3)+1: ix_underscore(3) + 6);
    D_id{d} = id;
end

unique_id = unique(D_id);
D_unique = zeros(n_D,1);

ix_match = zeros(n_D,1);
for d = 1:n_D;
    ix = strfind(unique_id, D{d}(ix_underscore(3)+1: ix_underscore(3) + 6));
    is = ~cellfun(@isempty, ix);
    D_unique(d) = find(is);
end

for d = 1:n_D;
    info_plot = D_unique(d);
    new_filename = sprintf('%s-%03.0f-%03.0f', 'TEST', info_site, info_plot);
    suffix =  D{d}( ix_underscore(3) + 6 + 1:end);
    filepath_new = sprintf('%s%s%s',dir_target, new_filename, suffix);
    filepath_old = sprintf('%s%s', dir_source, D{d});
    command = sprintf('cp %s %s', filepath_old, filepath_new);
    [~,~] = system(command);
end
