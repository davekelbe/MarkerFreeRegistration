%% Rename BIRD
%info_experiment = 'HAKT1';
%path_source = '/Volumes/Asterix/Original/CBL-HAK-BIRDS-T1-9-30-15/';

%info_experiment = 'HAKT2';
%path_source = '/Volumes/Asterix/Original/CBL-HAK-BIRDS-T2-12-1-15/';

%path_source = '/Volumes/Asterix/Original/CBL-KIP-BIRDS-10-20-15/';

%path_source = '/Volumes/Asterix/Original/CBL-KIP-BIRDS-11-6-15/';

path_source = '/Volumes/Asterix/Original/CBL-KIP-BIRDS-11-17-15/';



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
    if ~isempty(strfind(D{d}, 'JPG'));
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
    info_experiment = D{d}(1:ix_dash(2)-1);
    info_experiment = strrep(info_experiment, '-','');
    
    
    path_target = sprintf('%s%s-lidar%s', '/Volumes/Asterix/Renamed/', info_experiment, '/');
    if ~exist(path_target, 'dir');
        mkdir(path_target);
    end
    
    site = D{d}(ix_dash(2)+1:ix_dash(3)-1);
    site = str2double(site);
    site = sprintf('%03.0f', site);
    if ~isempty(ix_underscore)
        heading = D{d}(ix_dash(3)+1:ix_underscore-1);
    else
        heading = D{d}(ix_dash(3)+1:ix_dot(1)-1);
    end
    switch heading;
        case 'C';
            plot = 0;
        case 'N';
            plot = 1;
        case 'E';
            plot = 2;
        case 'S';
            plot = 3;
        case 'W';
            plot = 4;
    end
    filepath_old = sprintf('%s%s', path_source, D{d});
    if ~isempty(ix_underscore)
        extension = D{d}(ix_underscore:end);
    else
        extension = D{d}(ix_dot:end);
    end
    filepath_new = sprintf('%s%s-%s-%03.0f%s', path_target, info_experiment, site, plot,extension );
    command = sprintf('cp %s %s', filepath_old, filepath_new);
    [a,b,] = system(command);
    
end