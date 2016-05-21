% Author: Kevin Sacca
% Title: DEM_metrics2.m
% Purpose: Extract quantitative metrics from ground-based lidar point
%          clouds. First uses the DEM obtained from Dave Kelbe's code to
%          calculate the height of each point above the DEM. Then uses the
%          height of each point to segment layers of vegetation in order to
%          create separate data products to link with bird song
%          characteristics.

% Adjust these variables before running
nSites = 1; % Only have data for KIP3. Need other KIP datasets.
datestr = '03-8';

% Plot points up to 25 meters above DEM
z_height = 0:0.5:25;

% Set up empty arrays to be filled
z_freq = zeros(6,numel(z_height),3); % 5 captures. 6th = Average
plot_avgs = zeros(3,numel(z_height),nSites);

for siteIdx = [3] % Only processed site is KIP3 so far.
    for plotIdx = 1:3
        for scanIdx = 0:4
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/data_dem.mat'))
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/data_x.mat'))
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/data_y.mat'))
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/data_z.mat'))
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/dem_qx.mat'))
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/dem_qy.mat'))
            load(strcat('/Users/kevinsacca/Desktop/BIRDS/Processed/KIP',num2str(siteIdx),'/',datestr,'/00',num2str(plotIdx),'/00',num2str(scanIdx),'/mat/dem_qz.mat'))

            % Subtract DEM term from each height value (z-axis) to get height above DEM
            data_z0 = data_z(:) - data_dem(:);

            % Isolate all hits measured to be above the DEM
            data_x_ag = data_x(data_z0 >= 0);
            data_y_ag = data_y(data_z0 >= 0);
            data_z_ag = data_z0(data_z0 >= 0);

            % Isolate the hits measured to be the DEM
            data_z_dem = data_z0(abs(data_z0) < 0.25);
            % figure('Position',[100 400 560 420],'color', [1 1 1 ]);
            % mesh(dem_qx,dem_qy,dem_qz);
            % axis('equal');
            % xlabel('x, [meters]', 'fontsize', 14);
            % ylabel('y, [meters]', 'fontsize', 14);
            % zlabel('z, [meters]', 'fontsize', 14);
            % colormap('jet');
            % view(45,30);
            % grid off;
            % title('DEM', 'fontsize', 14)
            % colorbar
            
            % Normalize for oversampling 
            % Convert to spherical coordiantes
            % data_e is in [rad] from the xy plane
            [~,data_e,data_r] = cart2sph(data_x,data_y,data_z);
            % Take absolute value so don't have to worry about negative el
            data_e = abs(data_e);
            % Convert to degrees 
            data_e = rad2deg(data_e);
            % Compute number of samples at each elevation angle 
            % Using the form y = m*x + b
            data_nsamp = ((720-1)./90).*data_e + 1;
            % Take inverse to compute weight
            data_w = 1./data_nsamp;
            % Now multiply this weight in hist calculation (see line 84) 
        
            

            % Calculate statistics of hits above/below ground (ag / bg)
            n_hits = numel(data_z0);
            n_hits_ag = numel(data_z_ag);
            n_hits_bg = n_hits - n_hits_ag;

            % Plot above-DEM point cloud
            data_z_bin = floor(data_z0) + floor( (data_z0-floor(data_z0))/0.5) * 0.5;

            % Histogram of returns at 0.5m resolution height values.
            i = 1;
            for z = z_height
                tmp = data_z_bin;
                tmp(data_z_bin == z) = 1;
                tmp(data_z_bin ~= z) = 0;
                z_freq(scanIdx+1,i,plotIdx) = sum(tmp.*data_w);
                i = i+1;
            end
            % Find average height vs. freq
            z_freq(6,:,plotIdx) = mean(z_freq(:,:,plotIdx));

        end     % end of scan loop
        
        % Move averages to array
        plot_avgs(plotIdx,:,1) = z_freq(6,:,plotIdx);
        
    end     % end of plot loop
    
    % Calculate average height vs. freq for all plots for each site
    plot_avg = mean(plot_avgs(:,:,1));
    
end     % end of site loop


%% Plots

% Height vs. frequency of hits
for plotIdx = 1:3
    figure('Position', [100 400 860 620])
    hold on
    plot(z_freq(1,:,plotIdx), z_height, 'Color', [1,0.5,0])
    plot(z_freq(2,:,plotIdx), z_height, 'Color', [1,0.6,0])
    plot(z_freq(3,:,plotIdx), z_height, 'Color', [1,0.7,0])
    plot(z_freq(4,:,plotIdx), z_height, 'Color', [1,0.8,0])
    plot(z_freq(5,:,plotIdx), z_height, 'Color', [1,0.9,0])
    plot(z_freq(6,:,plotIdx), z_height, 'Color', [0,0.5,1])
    hold off
    title('Frequency of Hits at Height above DEM', 'fontSize', 14, 'fontWeight','bold')
    xlabel('Frequency of Hits', 'fontSize', 14)
    ylabel('Height above DEM [m]', 'fontSize', 14)
    legend('Center','North','East','South','West','Average','Location','NorthEast')
    set(gca,'fontName','Helvetica','fontSize',10,'TickDir','out')
    
end
    
% Plot average height vs. freq. hits for three sites
figure('Position', [100 400 860 620])
hold on
plot(z_freq(6,:,1), z_height, 'Color', [1,0.5,0])
plot(z_freq(6,:,2), z_height, 'Color', [1,0.7,0])
plot(z_freq(6,:,3), z_height, 'Color', [1,0.9,0])
plot(plot_avg, z_height, 'Color', [0,0.5,1])
hold off
title('Frequency of Hits at Height above DEM', 'fontSize', 14, 'fontWeight','bold')
xlabel('Frequency of Hits', 'fontSize', 14)
ylabel('Height above DEM [m]', 'fontSize', 14)
legend('Plot 001','Plot 002','Plot 003', 'Average','Location','NorthEast')
set(gca,'fontName','Helvetica','fontSize',10,'TickDir','out')













