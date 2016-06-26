
% Load Lookup file for TLS data
path_soaprootTLSlookup = 'D:\Users\djk2312\Documents\Paper01\Canopy\Data\California\TLS\';
filepath_soaprootTLSlookup = sprintf('%s%s',path_soaprootTLSlookup, 'SoaprootLiDAR.csv');
fid = fopen(filepath_soaprootTLSlookup);
TLSLookup = textscan(fid, '%d %d %s', 'HeaderLines', 1, 'Delimiter', ',');
TLS.site = TLSLookup{1};
TLS.plot = TLSLookup{2}; 
TLS.filename = TLSLookup{3};
fclose all

n_scans = numel(TLS.site);
warning('off', 'arguments:exteriordata');

for s = 1:n_scans;
    site = TLS.site(s);
    if site~=331;
        continue
    end
    plot = TLS.plot(s);
    fprintf('\nSite %d Plot %d\n', site, plot)
    lidar(site,plot,'tree', 'soaproot');
end





