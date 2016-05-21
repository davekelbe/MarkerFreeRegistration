%% Levenberg Marquardt 

%Load data 
R_13_13 = eye(3);
t_13_13 = zeros(3,1);
filename_R_15_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\15\mat\R_031-15-13.mat';
load(filename_R_15_13);
R_15_13 = R;
filename_t_15_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\15\mat\t_031-15-13.mat';
load(filename_t_15_13);
t_15_13 = t;
filename_R_08_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\08\mat\R_031-08-13.mat';
load(filename_R_15_13);
R_08_13 = R;
filename_t_08_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\08\mat\t_031-08-13.mat';
load(filename_t_08_13);
t_08_13 = t;
filename_m1_15_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\15\mat\m1_best_031-15-13.mat';
load(filename_m1_15_13);
m_13_15 = m1_best;
filename_m2_15_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\15\mat\m2_best_031-15-13.mat';
load(filename_m2_15_13);
m_15_13 = m2_best;
filename_m1_08_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\08\mat\m1_best_031-08-13.mat';
load(filename_m1_08_13);
m_13_08 = m1_best;
filename_m2_08_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\08\mat\m2_best_031-08-13.mat';
load(filename_m2_08_13);
m_08_13 = m2_best;
filename_tree_13 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\13\mat\tree.mat';
load(filename_tree_13);
tree_13 = tree;
filename_tree_08 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\08\mat\tree.mat';
load(filename_tree_08);
tree_08 = tree;
filename_tree_15 = 'D:\Users\djk2312\Documents\Harvard\03-01\031\15\mat\tree.mat';
load(filename_tree_15);
tree_15 = tree;
clear R t m1_best m2_best tree
%%
n_m_15 = numel(m_13_15);
n_m_08 = numel(m_13_08);

T_13_15 = tree_13(m_13_15);
X_13_15 = zeros(3,n_m_15);
D_13_15 = zeros(3,n_m_15);
for i = 1:n_m_15;
    X_13_15(:,i) = T_13_15(i).loc(:,1);
    D_13_15(:,i) = T_13_15(i).r(1);
end 

T_15_13 = tree_15(m_15_13);
X_15_13 = zeros(3,n_m_15);
D_15_13 = zeros(3,n_m_15);
for i = 1:n_m_15;
    X_15_13(:,i) = T_15_13(i).loc(:,1);
    D_15_13(:,i) = T_15_13(i).r(1);
end

T_13_08 = tree_13(m_13_08);
X_13_08 = zeros(3,n_m_08);
D_13_08 = zeros(3,n_m_08);
for i = 1:n_m_08;
    X_13_08(:,i) = T_13_08(i).loc(:,1);
    D_13_08(:,i) = T_13_08(i).r(1);
end 

T_08_13 = tree_08(m_08_13);
X_08_13 = zeros(3,n_m_08);
D_08_13 = zeros(3,n_m_08);
for i = 1:n_m_08;
    X_08_13(:,i) = T_08_13(i).loc(:,1);
    D_08_13(:,i) = T_08_13(i).r(1);
end

%r_08_13 = vrrotmat2vec(R_08_13);
[rx_08_13 ry_08_13, rz_08_13] = decompose_rotation(R_08_13);
[rx_15_13 ry_15_13, rz_15_13] = decompose_rotation(R_15_13);

%% Levenberg Marquardt

x0 = [X_13_15(1,:),X_13_15(2,:),X_13_15(3,:),...
    X_13_08(1,:),X_13_08(2,:),X_13_08(3,:),...
    X_15_13(1,:),X_15_13(2,:),X_15_13(3,:),...
    X_08_13(1,:),X_08_13(2,:),X_08_13(3,:),...
    rx_15_13 ry_15_13 rz_15_13,...
    rx_08_13 ry_08_13 rz_08_13,...
    t_15_13',...
    t_08_13'];
x = lsqnonlin(@LM_demo,x0);

X_13_15_LM = [x(1:20);x(21:40);x(41:60)];
X_13_08_LM = [x(61:88);x(89:116);x(117:144)];
X_15_13_LM = [x(145:164);x(165:184);x(185:204)];
X_08_13_LM = [x(205:232);x(233:260);x(261:288)];
rx_15_13_LM = x(289);
ry_15_13_LM = x(290);
rz_15_13_LM = x(291);
rx_08_13_LM = x(292);
ry_08_13_LM = x(293);
rz_08_13_LM = x(294);
t_15_13_LM = x(295:297)';
t_08_13_LM = x(298:300)';
R_15_13_LM = compose_rotation(rx_15_13_LM,ry_15_13_LM,rz_15_13_LM);
R_08_13_LM = compose_rotation(rx_08_13_LM,ry_08_13_LM,rz_08_13_LM);

Y_15_13 = R_15_13 * X_15_13 + repmat(t_15_13,[1,n_m_15]);
Y_08_13 = R_08_13 * X_08_13 + repmat(t_08_13,[1,n_m_08]);
Y_15_13_LM = R_15_13_LM * X_15_13_LM + repmat(t_15_13_LM,[1,n_m_15]);
Y_08_13_LM = R_08_13_LM * X_08_13_LM + repmat(t_08_13_LM,[1,n_m_08]);

%% Visualization
[Sx, Sy, Sz] = sphere(5);
figure;
hold on
% Original matches 
for i = 1:n_m_15;
    Px = Sx.*D_13_15(i)+X_13_15(1,i);
    Py = Sy.*D_13_15(i)+X_13_15(2,i);
    Pz = Sz.*D_13_15(i)+X_13_15(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[1 0 0],'edgecolor',[.5 0 0],'facealpha',.25);
end
for i = 1:n_m_15;
    Px = Sx.*D_13_15(i)+X_13_15_LM(1,i);
    Py = Sy.*D_13_15(i)+X_13_15_LM(2,i);
    Pz = Sz.*D_13_15(i)+X_13_15_LM(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[0 1 0],'edgecolor',[0 .5 0],'facealpha',.25);
end

for i = 1:n_m_08;
    Px = Sx.*D_13_08(i)+X_13_08(1,i);
    Py = Sy.*D_13_08(i)+X_13_08(2,i);
    Pz = Sz.*D_13_08(i)+X_13_08(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[0 0 1],'edgecolor',[0 0 .5],'facealpha',.25);
end
for i = 1:n_m_08;
    Px = Sx.*D_13_08(i)+X_13_08_LM(1,i);
    Py = Sy.*D_13_08(i)+X_13_08_LM(2,i);
    Pz = Sz.*D_13_08(i)+X_13_08_LM(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[0 1 1],'edgecolor',[0 .5 .5],'facealpha',.25);
end

% Transformed points 
for i = 1:n_m_15;
    Px = Sx.*D_13_15(i)+X_15_13(1,i);
    Py = Sy.*D_13_15(i)+X_15_13(2,i);
    Pz = Sz.*D_13_15(i)+X_15_13(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[1 1 0],'edgecolor',[.5 .5 0],'facealpha',.25);
end
for i = 1:n_m_15;
    Px = Sx.*D_13_15(i)+X_15_13_LM(1,i);
    Py = Sy.*D_13_15(i)+X_15_13_LM(2,i);
    Pz = Sz.*D_13_15(i)+X_15_13_LM(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[1 0 1],'edgecolor',[.5 0 .5],'facealpha',.25);
end

for i = 1:n_m_08;
    Px = Sx.*D_13_08(i)+X_08_13(1,i);
    Py = Sy.*D_13_08(i)+X_08_13(2,i);
    Pz = Sz.*D_13_08(i)+X_08_13(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[.25 0 .75],'edgecolor',[.125 0 .375],'facealpha',.25);
end
for i = 1:n_m_08;
    Px = Sx.*D_13_08(i)+X_08_13_LM(1,i);
    Py = Sy.*D_13_08(i)+X_08_13_LM(2,i);
    Pz = Sz.*D_13_08(i)+X_08_13_LM(3,i);
    h = surf(Px,Py,Pz);
    set(h,'facecolor',[0.75 0 .25],'edgecolor',[.375 0 .25],'facealpha',.25);
end

axis auto 
axis equal
foo = 1;









