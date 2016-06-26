function [ ] = fig_match_tree_neon_color_dbh( site_all, plot_all, match_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Plot results

figure;

n_site_all = numel(site_all);
color = jet(n_site_all);
%legend_text = cell(n_site_all,1);
hold on; 

x = [];
y = [];
for i = 1:n_site_all;
    x = [x ;2*match_all{i}.NEON_neon_r];
    y = [y ;2*match_all{i}.NEON_tree_r];
end

x_step = .1;
x_array = 0:x_step:max(x);
n_x = numel(x_array)-1;
color = jet(n_x);
counter = 1;
for i = 1:n_x;
    is_valid = x>x_array(i)&x<=x_array(i+1);
   if sum(is_valid)>0;
        legend_text{counter} = sprintf('%2.1f - %2.1f m',x_array(i),x_array(i+1));
        plot(x(is_valid),y(is_valid),'x',...
            'markersize',10,'linewidth',2,'markeredgecolor',color(i,:));
        counter = counter + 1;
   end
end

xmax = max([x;y]);
p = polyfit(x,y,1);
yfit  = polyval(p,x);
%yfit = p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;

legend_text{counter} = sprintf('Fit');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);
legend_text{counter+1} = sprintf('1:1');%y = %3.2fx + %3.2f, R^2 = %3.2f',p(2), p(1), rsq);

hold on;
plot([0 xmax],[0 xmax],'-k');
plot([0 xmax],[polyval(p,0) polyval(p,xmax)],':k')
htxt = annotation('textbox',[.6 .8 .1 .1], 'string', sprintf('R^2 = %4.2f',rsq));
set(htxt,'edgecolor','w');

title('DBH Retrieval ');
hleg = legend(legend_text);
set(hleg, 'edgecolor','w');
set(hleg, 'location', 'bestoutside');
hlegtitle = get(hleg,'title');
set(hlegtitle,'string', 'DBH');

xlim([0 xmax ]);
ylim([0 xmax]);
ylabel('Detected [m]');
xlabel('NEON truth data [m]');
axis equal 

end

