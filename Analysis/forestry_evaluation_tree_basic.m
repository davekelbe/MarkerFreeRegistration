function [ n_stems, ba, atdbh, match ] = forestry_evaluation_tree_basic( tree, NEON_DBH_10over, NEON_DBH_10under )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Stem Map 

  figure('position', [587 95 1026 832]);
    hold on
    for s = 1:numel(atdbh.ROI_x);
        h = filledCircle([atdbh.ROI_x(s) atdbh.ROI_y(s)],atdbh.ROI_r(s),1000,'b');
        set(h,'FaceAlpha',.5);
    end
    for s = 1:numel(atdbh.SEG_x);
        h = filledCircle([atdbh.SEG_x(s) atdbh.SEG_y(s)],atdbh.SEG_r(s),1000,'r');
        set(h,'FaceAlpha',.5);
    end
    for s = 1:numel(atdbh.NEON_x)
        h = filledCircle([atdbh.NEON_x(s) atdbh.NEON_y(s)],atdbh.NEON_r(s),1000,'g');
        set(h,'FaceAlpha',.5);
    end
    h = rectangle('position',[-10 -10 20 20], 'Curvature',[0 0]);
    axis equal;
    axis([-16 16 -16 16]);
    xlabel('Easting');
    ylabel('Northing');
    grid on
    print(gcf,'-depsc','-opengl',filepath_dbhmap)



unique_NEON_x = NEON_DBH_10over.x;
unique_NEON_y = NEON_DBH_10over.y;
unique_NEON_xy = sqrt(unique_NEON_x.^2 + unique_NEON_y.^2);
unique_NEON_r = NEON_DBH_10over.dbh/(2*100);
n_unique_NEON = numel(unique_NEON_x);

%% Stem Count

r_step = 2;
rmax = max([unique_ROI_xy; unique_SEG_xy; 10]);
r_array = 0:r_step:rmax+r_step;
n_r = numel(r_array)-1;
n_stems.SEG_of_xy = zeros(n_r,1);
n_stems.ROI_of_xy = zeros(n_r,1);

%% Basal Area

%% DBH 


end

