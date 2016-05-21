%% Rename HAVO 

path_source = '/Volumes/Asterix/Original/CBL-HAVO-11-25-15/CBL-HAVO-Photos-11-25-15/';
path_target = '/Volumes/Asterix/Renamed/HAVO-photos/';

if ~exist(path_target);
    mkdir(path_target);
end

D = dir(path_source);
D = remove_hiddenfiles(D);

site = '000';

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
    plot = 11*row + col;
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