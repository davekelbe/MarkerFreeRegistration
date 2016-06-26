%% Allometry: Stem height vs dbh 

filepath_HF_fielddata = 'D:\Data\Harvard_Forest\Fieldwork_2012\Ben_VegStructure_2012_HF_Aplot_data.csv';
fid = fopen(filepath_HF_fielddata);
HF_field = textscan(fid, '%s%s%s%u%f%s%f%f%f%f%s%s%s%s%s%s%s%s%f%f%f%f%u', 'delimiter', ',', 'headerlines',1); 
fclose all

HF_true.dbh = HF_field{5};
HF_true.height = HF_field{9};
HF_true.abh = (pi.*HF_true.dbh.^2)./4;

isvalid = (HF_true.dbh~=99);
HF_true.dbh = HF_true.dbh(isvalid);
HF_true.height = HF_true.height(isvalid);
HF_true.abh = HF_true.abh(isvalid);

xstep = linspace(min(HF_true.dbh), max(HF_true.dbh), 100);
% DBH vs area
%{
figure; 
scatter(log(HF_true.dbh), HF_true.height, '.b');
xlabel('DBH [cm]');
ylabel('Height [m]');
%}
% Figure: Height vs area
%{
figure; 
scatter(HF_true.abh, HF_true.height, '.b');
xlabel('Area at BH [cm]');
ylabel('Height [m]');
%}

% Split training and testing data 
n_all = numel(HF_true.dbh);
ix = randperm(n_all);
p_train = 20;
n_train = floor((p_train./100)*n_all);
ix_valid_train = ix(1:n_train);
is_valid_train = false(n_all,1);
is_valid_train(ix_valid_train) = true;
is_valid_test = ~is_valid_train;
% Figure: test and train 
%{
figure; 
scatter(log(HF_true.dbh(is_valid_train)), HF_true.height(is_valid_train), '.b');
hold on
scatter(log(HF_true.dbh(is_valid_test)), HF_true.height(is_valid_test), '.r');
xlabel('DBH [cm]');
ylabel('Height [m]');
%}

%% Training - Log 
x_train = HF_true.dbh(is_valid_train);
y_train = HF_true.height(is_valid_train);
x_test = HF_true.dbh(is_valid_test);
y_test = HF_true.height(is_valid_test);
btest = [1 0]; 
[p, stats, t0slope, t0int, tcrit ] = regress_with_test( 1./(x_train),log(y_train),btest );
if abs(t0int) > tcrit;
    fprintf('\n MST Reject null hypothesis intercept\n');
else
    fprintf('\n MST Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n MST Reject null hypothesis slope\n');
else
    fprintf('\n MST Fail to reject null hypothesis slope\n');
end

yfit_train = exp(polyval(p,1./xstep));

figure
ms = 10;
clear legend_str
legend_ctr = 1;
plot(1./x_train,log(y_train),'o', 'markersize', ms/2, 'color', 'b', 'markerfacecolor', [.85 .85 .85]);
hold on
legend_str{legend_ctr} = 'Train';
legend_ctr = legend_ctr + 1;
plot(1./x_test,log(y_test),'o', 'markersize', ms/2, 'color', 'r', 'markerfacecolor', [.85 .85 .85]);
legend_str{legend_ctr} = 'Test';
legend_ctr = legend_ctr + 1;
%plot(xstep,yfit_train,'o-', 'markersize', ms/2, 'color', 'b', 'markerfacecolor', [.85 .85 .85]);
%legend_str{legend_ctr} = 'Train Fit';
%legend_ctr = legend_ctr + 1;
xlabel('1 / DBH [cm]');
ylabel('log(Height [m])');

figure
ms = 10;
clear legend_str
legend_ctr = 1;
plot(x_train,y_train,'o', 'markersize', ms/2, 'color', 'b', 'markerfacecolor', [.85 .85 .85]);
hold on
legend_str{legend_ctr} = 'Train';
legend_ctr = legend_ctr + 1;
plot(x_test,y_test,'o', 'markersize', ms/2, 'color', 'r', 'markerfacecolor', [.85 .85 .85]);
legend_str{legend_ctr} = 'Test';
legend_ctr = legend_ctr + 1;
plot(xstep,yfit_train,'o-', 'markersize', ms/2, 'color', 'b', 'markerfacecolor', [.85 .85 .85]);
legend_str{legend_ctr} = 'Train Fit';
legend_ctr = legend_ctr + 1;
xlabel('DBH [cm]');
ylabel('Height [m]');
%% SQRT model 
%
x_train = HF_true.dbh(is_valid_train);
y_train = HF_true.height(is_valid_train);
x_test = HF_true.dbh(is_valid_test);
y_test = HF_true.height(is_valid_test);
btest = [1 0]; 
[p, stats, t0slope, t0int, tcrit ] = regress_with_test( sqrt(x_train),y_train,btest );
if abs(t0int) > tcrit;
    fprintf('\n MST Reject null hypothesis intercept\n');
else
    fprintf('\n MST Fail to reject null hypothesis intercept\n');
end
if abs(t0slope) > tcrit;
    fprintf('\n MST Reject null hypothesis slope\n');
else
    fprintf('\n MST Fail to reject null hypothesis slope\n');
end

yfit_train = polyval(p,sqrt(xstep));

figure
ms = 10;
clear legend_str
legend_ctr = 1;
plot(x_train,y_train,'o', 'markersize', ms/2, 'color', 'b', 'markerfacecolor', [.85 .85 .85]);
hold on
legend_str{legend_ctr} = 'Train';
legend_ctr = legend_ctr + 1;
plot(x_test,y_test,'o', 'markersize', ms/2, 'color', 'r', 'markerfacecolor', [.85 .85 .85]);
legend_str{legend_ctr} = 'Test';
legend_ctr = legend_ctr + 1;
plot(xstep,yfit_train,'o-', 'markersize', ms/2, 'color', 'b', 'markerfacecolor', [.85 .85 .85]);
legend_str{legend_ctr} = 'Train Fit';
legend_ctr = legend_ctr + 1;
xlabel('DBH [cm]');
ylabel('Height [m]');
%}
%% Testing 
