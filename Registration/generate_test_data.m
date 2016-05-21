function [ P_LCS,P_rad, P_plot,P_n_tree, n_S ] = generate_test_data(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_S = 2;
P_LCS = cell(n_S,1);
P_rad = cell(n_S,1);

P_n_tree = [5,4]';
P_LCS{1} = zeros(3,P_n_tree(1));
a = -10; b = 10;
P_LCS{1}(1,:) = a+ (b-a).*rand(P_n_tree(1),1);
P_LCS{1}(2,:) = a+ (b-a).*rand(P_n_tree(1),1);
a = -.5; b = .5;
P_LCS{1}(3,:) = a+ (b-a).*rand(P_n_tree(1),1);
a = 0.1; b = 0.6;
P_rad{1} = a+ (b-a).*rand(P_n_tree(1),1);

P_LCS{2} = zeros(3,P_n_tree(2));
a = -10; b = 10;
P_LCS{2}(1,:) = a+ (b-a).*rand(P_n_tree(2),1);
P_LCS{2}(2,:) = a+ (b-a).*rand(P_n_tree(2),1);
a = -.5; b = .5;
P_LCS{2}(3,:) = a+ (b-a).*rand(P_n_tree(2),1);
a = 0.1; b = 0.6;
P_rad{2} = a+ (b-a).*rand(P_n_tree(2),1);

ix_match1 = [1 3 4 5];
ix_match2 = [4 3 2 1];

n_match = numel(ix_match1);

t = [-5 0 0]';
R = eye(3);

P_LCS{2}(:,ix_match2) = R*P_LCS{1}(:,ix_match1) + repmat(t, [1,n_match]);

% Add noise 
a = -.1; b =.1;
noise = [a+ (b-a).*rand(P_n_tree(2),1)';...
    a+ (b-a).*rand(P_n_tree(2),1)';...
    a+ (b-a).*rand(P_n_tree(2),1)'];
    
P_LCS{2} = P_LCS{2} + noise;

P_rad{2}(ix_match2) = P_rad{1}(ix_match1);
a = -.01; b =.01;
noise = a+ (b-a).*rand(P_n_tree(2),1);
P_rad{2} = P_rad{2} + noise;

P_plot =[8 13]';

end

