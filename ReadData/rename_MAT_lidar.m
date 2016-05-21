%% Rename MAT
info_experiment = 'MAT';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-SPE-800-934-9-11-15/SPE800-9-11-15/';
%site = '001';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-SPE-800-934-9-11-15/SPE934-9-11-15/';
%site = '002';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-SPE-1064-10-30-15/';
%site = '003';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-SPE-1116-9-29-15/';
%site = '004';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-WPL-1116-11-3-15/';
%site = '005';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-WPL-1204-1274-9-24-15/WPL1204-9-24-15/';
%site = '006';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-WPL-1204-1274-9-24-15/WPL1274-9-24-15/';
%site = '007';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-HAK-1468-1600-10-27-15/MATHAK1468-Lidar/';
%site = '008';

path_source = '/Volumes/Asterix/Original/CBL-MAT-HAK-1468-1600-10-27-15/MATHAK1600-Lidar/';
site = '009';

path_target = '/Volumes/Asterix/Renamed/MAT-lidar/';
if ~exist(path_target, 'dir');
    mkdir(path_target);
end


D = dir(path_source);
D = remove_hiddenfiles(D);


n_D = numel(D);
for d = 1:n_D;
    if ~isempty(strfind(D{d}, 'zip'));
        continue
    end
    if ~isempty(strfind(D{d}, 'junk'));
        continue
    end
    if ~isempty(strfind(D{d}, 'Junk'));
        continue
    end
    if ~isempty(strfind(D{d}, 'photos'));
        continue
    end
    if ~isempty(strfind(D{d}, 'Photos'));
        continue
    end
    if ~isempty(strfind(D{d}, 'Fun'));
        continue
    end
    if ~isempty(strfind(D{d}, 'jpg'));
        continue
    end
    if ~isempty(strfind(D{d}, 'Mek'));
        continue
    end
    if ~isempty(strfind(D{d}, 'FUN'));
        continue
    end
    ix_dash = strfind(D{d}, '-');
    ix_dot = strfind(D{d}, '.');
    ix_underscore = strfind(D{d}, '_');
    row = D{d}(ix_dash(1)+1:ix_dash(2)-1);
    if ~isempty(ix_underscore)
        col = D{d}(ix_dash(2)+1:ix_underscore-1);
    else
        col = D{d}(ix_dash(2)+1:ix_dot(1)-1);
    end
    row = str2double(row);
    col = str2double(col);
    plot = 5*(row-1) + col;
    filepath_old = sprintf('%s%s', path_source, D{d});
    if ~isempty(ix_underscore)
        extension = D{d}(ix_underscore:end);
    else
        extension = D{d}(ix_dot:end);
    end
    filepath_new = sprintf('%s%s-%s-%03.0f%s', path_target, 'MAT', site, plot,extension );
    command = sprintf('cp %s %s', filepath_old, filepath_new);
    [a,b,] = system(command);
    
end