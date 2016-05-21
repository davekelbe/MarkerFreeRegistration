%% Rename HAVO 

%path_source = '/Volumes/Asterix/Original/CBL-MAT-SPE-1064-10-30-15/MAT-SPE1064-Photos-10-30-15/';
%site = '003';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-WPL-1116-11-3-15/MAT-WPL1116-Photos-11-3-15/';
%site = '005';

%path_source = '/Volumes/Asterix/Original/CBL-MAT-HAK-1468-1600-10-27-15/MAT-HAK1468-Photos-10-27-15/';
%site = '008';

path_source = '/Volumes/Asterix/Original/CBL-MAT-HAK-1468-1600-10-27-15/MAT-HAK1600-Photos-10-27-15/';
site = '009';

path_target = '/Volumes/Asterix/Renamed/MAT-photos/';

if ~exist(path_target);
    mkdir(path_target);
end

D = dir(path_source);
D = remove_hiddenfiles(D);


n_D = numel(D);
for d = 1:n_D;
    ix_dash = strfind(D{d}, '-');
    %ix_dot = strfind(D{d}, '.');
    %ix_underscore = strfind(D{d}, '_');
    row = D{d}(ix_dash(1)+1:ix_dash(2)-1);
  %  if ~isempty(ix_underscore)
  %      col = D{d}(ix_dash(2)+1:ix_underscore-1);
  %  else
        col = D{d}(ix_dash(2)+1:ix_dash(3)-1);
  %  end
    row = str2double(row);
    col = str2double(col);
    plot = 5*(row-1) + col;
    filepath_old = sprintf('%s%s', path_source, D{d});
   % if ~isempty(ix_underscore)
        extension = D{d}(ix_dash(3):end);
   % else
   %     extension = D{d}(ix_dot:end);
   % end
    filepath_new = sprintf('%s%s-%s-%03.0f%s', path_target, 'HAVO', site, plot,extension );
    command = sprintf('cp %s %s', filepath_old, filepath_new);
    [a,b,] = system(command);
    
end