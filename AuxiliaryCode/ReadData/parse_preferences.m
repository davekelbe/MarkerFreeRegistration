function [ mypreferences ] = parse_preferences( filepath_preferences )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(filepath_preferences);
    mypreferences = [];
    return
end
fid = fopen(filepath_preferences);
if fid <0;
    mypreferences = [];
    return
end

preferences = textscan(fid, '%s', 'Delimiter', '\n');
preferences = preferences{1};

is_empty = cellfun(@isempty,preferences);
preferences = preferences(~is_empty);

search = 'path_up';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
path_up = preferences(is);
mypreferences.path_up = path_up{1}(numel(search)+2:end);

search = 'path_source';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
path_source = preferences(is);
mypreferences.path_source = path_source{1}(numel(search)+2:end);

search = 'info_experiment';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
info_experiment = preferences(is);
mypreferences.info_experiment = info_experiment{1}(numel(search)+2:end);

search = 'info_sites';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
info_sites = preferences(is);
mypreferences.info_sites = info_sites{1}(numel(search)+2:end);

search = 'info_plots';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
info_plots = preferences(is);
mypreferences.info_plots = info_plots{1}(numel(search)+2:end);

search = 'info_suffix';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
info_suffix = preferences(is);
mypreferences.info_suffix = info_suffix{1}(numel(search)+2:end);

search = 'options_skipseg';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
info_suffix = preferences(is);
mypreferences.options_skipseg = info_suffix{1}(numel(search)+2:end);

end

