savepath = '/Users/Kelbe/Documents/Research/Experiments/HarvardForest/Field Data/DBH Maps/';

file.input = '/Users/Kelbe/Documents/Research/Experiments/HarvardForest/Field Data/NEON_DBH_10over.csv';
fid  = fopen(file.input);
C = textscan(fid,'%s%s%s%d%f%s%f%f%f%f%s%s%s%s%s%s%s%s%f%f%f%f%d','HeaderLines',1, 'delimiter', ',');

D.easting = C{1};
D.northing = C{2};
D.plot_id = C{3};
D.module = C{4};
D.dbh = C{5};
D.trupulse_position = C{6};
D.stem_distance = C{7};
D.stem_angle = C{8};
D.stem_height = C{9};
D.canopy_diameter = C{10};
D.vertical_position = C{11};
D.common_species = C{12};
D.other_species = C{13};
D.stem_notes = C{14};
D.stem_status = C{15};
D.status_notes = C{16};
D.shotgun_sampling = C{17};
D.shotgun_id = C{18};
D.canopy_diam_90deg = C{19};
D.canopy_diam_150deg = C{20};
D.canopy_diam_210deg = C{21};
D.point_id = C{22};
D.dbhm = D.dbh/100;

n.entry = size(D.easting,1);

D.x = zeros(n.entry,1);
D.y = zeros(n.entry,1); 

D.offset_x = zeros(n.entry,1);
D.offset_y = zeros(n.entry,1);

ix = strcmp(D.trupulse_position,repmat('Plot center',n.entry,1));
D.offset_x(ix) = 0;
D.offset_y(ix) = 0;
ix = strcmp(D.trupulse_position,repmat('SE corner',n.entry,1));
D.offset_x(ix) = 10;
D.offset_y(ix) = -10;
ix = strcmp(D.trupulse_position,repmat('SW corner',n.entry,1));
D.offset_x(ix) = -10;
D.offset_y(ix) = -10;
ix = strcmp(D.trupulse_position,repmat('NW corner',n.entry,1));
D.offset_x(ix) = -10;
D.offset_y(ix) = 10;
ix = strcmp(D.trupulse_position,repmat('NE corner',n.entry,1));
D.offset_x(ix) = 10;
D.offset_y(ix) = 10; 

D.x = D.offset_x  -(D.stem_distance+D.dbhm/2).*sind(D.stem_angle);
D.y = D.offset_y + (D.stem_distance+D.dbhm/2).*cosd(D.stem_angle);

plot.id = unique(D.plot_id);
n.plot = size(plot.id,1);
plot.ntree = zeros(n.plot,1);
 
for p = 1:n.plot;
    ix = strcmp(D.plot_id,repmat(plot.id{p},n.entry,1));
    plot.ntree(p) = sum(ix);  
    
    ix = find(ix);
    figure;
    hold on
    titlestr = sprintf('DBH map - %s',plot.id{p});
    title(titlestr);
    for t = 1:plot.ntree(p)
        i = ix(t);
        rectangle('Position',[D.x(i)-D.dbhm(i)/2 D.y(i)-D.dbhm(i)/2 D.dbhm(i) D.dbhm(i)],...
            'Curvature',[1 1]);
    end
    axis equal;
    axis([-12 12 -12 12]);
    xlabel('Easting');
    ylabel('Northing');
    grid on
    name = sprintf('DBHMap_%s',plot.id{p});
    filename = sprintf('%sEPS/%s.eps',savepath,name);
    print(gcf, '-deps', filename);
    filename = sprintf('%sPNG/%s.png',savepath,name);
    print(gcf, '-dpng', filename);
    filename = sprintf('%sNEON_DBH_10over.mat',savepath);
    save(filename,'D');
end
            