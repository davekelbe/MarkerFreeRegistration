% Point cloud figures for paper 


% Initial transformation of PLY

load('D:\Users\djk2312\Documents\Harvard\reg\031\25\mat\match12_RMSE.mat');
load('D:\Users\djk2312\Documents\Harvard\reg\031\25\mat\match1_R.mat');
load('D:\Users\djk2312\Documents\Harvard\reg\031\25\mat\match1_t.mat');

filepath_plyJ = 'D:\Users\djk2312\Documents\Harvard\reg\031\13\ply\points_full_031-13.ply';
%[vertexI, ~] = read_ply(filepath_plyI);
[vertexJ, ~] = read_ply(filepath_plyJ);
    data_xJ = vertexJ(:,1);
    data_yJ = vertexJ(:,2);
    data_zJ = vertexJ(:,3);
    colorJ = vertexJ(:,4:6);
    
filepath_ply_reg = 'D:\Users\djk2312\Documents\Harvard\reg\031\25\ply\points_full_031-13r.ply';
i = 25;
j = 13;
xyz2t = (match1_R{i,j}*[data_xJ data_yJ data_zJ]') + ...
        repmat(match1_t{i,j},1,numel(data_xJ));
write2ply(filepath_ply_reg,xyz2t', colorJ);


% color red and blue 

colorJ_base = [239 206 203]./255;
colorI_base = [204 221 238]./255;


filepath_ply_reg_red = 'D:\Users\djk2312\Documents\Harvard\reg\031\25\ply\points_full_031-13r_red.ply';
color_red = colorJ.*repmat(colorJ_base, [numel(data_xJ), 1]);
write2ply(filepath_ply_reg_red,xyz2t', color_red);

filepath_plyI = 'D:\Users\djk2312\Documents\Harvard\reg\031\25\ply\points_full_031-25.ply';
[vertexI, ~] = read_ply(filepath_plyI);
    data_xI = vertexI(:,1);
    data_yI = vertexI(:,2);
    data_zI = vertexI(:,3);
    colorI = vertexI(:,4:6);
filepath_ply_reg_blue = 'D:\Users\djk2312\Documents\Harvard\reg\031\25\ply\points_full_031-25_blue.ply';
color_blue = colorI.*repmat(colorI_base, [numel(data_xI), 1]);
xyz = [data_xI data_yI data_zI]';
write2ply(filepath_ply_reg_blue,xyz', color_blue);
