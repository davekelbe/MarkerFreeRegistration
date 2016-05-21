function [I_block_xyz, xyz_b ] = andrieu_block_xyz( data_x,data_y,data_z,data_a, data_e, data_n,t_r )
%UNTITLED Creates an Andrieu Image from data showing blocks
%   
%   Outputs:
%   I_block_xyz is the Andrieu image block, colored according to block #
%   xyz_b  is a LUT giving the xb,yb,zb block numbers for the image (pixel
%       value) b
%
%
%   (C) David Kelbe, Rochester Institute of Technology 

block_x_lo = floor(min(data_x)):t_r:max(data_x);
block_x_hi = circshift(block_x_lo,[0,-1]); 
block_x_hi(end) = block_x_lo(end) + t_r;
block_y_lo = floor(min(data_y)):t_r:max(data_y);
block_y_hi = circshift(block_y_lo,[0,-1]); 
block_y_hi(end) = block_y_lo(end) + t_r;
block_z_lo = floor(min(data_z)):t_r:max(data_z);
block_z_hi = circshift(block_z_lo,[0,-1]); 
block_z_hi(end) = block_z_lo(end) + t_r;
n_x_block = numel(block_x_lo);
n_y_block = numel(block_y_lo);
n_z_block = numel(block_z_lo);
n_xyz_block = n_x_block*n_y_block*n_z_block;

n_data = numel(data_x);
%data_xb = zeros(n_data,1);
%data_yb = zeros(n_data,1);
data_xyzb = zeros(n_data,1);

%xb_array = 1:n_x_block;
%yb_array = 1:n_y_block;
xyzb_array = 1:n_xyz_block;

%xb_array_rand = xb_array(randperm(n_x_block)); %Random colors
%yb_array_rand = yb_array(randperm(n_y_block)); %Random colors
ix_rand = randperm(n_xyz_block);
xyzb_array_rand = xyzb_array(ix_rand);
xyz_b = zeros(n_xyz_block,3);

b = 1;
for xb = 1:n_x_block;
    for yb = 1:n_y_block;
        for zb = 1:n_z_block
          % fprintf('\n%05.0f of %05.0f\n',b,n_xyz_block);
        x_lo = block_x_lo(xb);
        x_hi = block_x_hi(xb);
        y_lo = block_y_lo(yb);
        y_hi = block_y_hi(yb);
        z_lo = block_z_lo(zb);
        z_hi = block_z_hi(zb);
        is_block = data_x>=x_lo & data_x<x_hi&...
                data_y>=y_lo & data_y<y_hi&...
                data_z>=z_lo & data_z<z_hi;
       % data_xb(is_block) = xb_array_rand(xb);
        %data_yb(is_block) = yb_array_rand(yb);
        data_xyzb(is_block) = xyzb_array_rand(b);
        xyz_b(xyzb_array_rand(b),1) = xb;
        xyz_b(xyzb_array_rand(b),2) = yb;
        xyz_b(xyzb_array_rand(b),3) = zb;
        b = b + 1;
        end
    end
end


%[I_block_x,~,~] = andrieu_image_from_points2(data_a,data_e,data_n,data_xb);
%[I_block_y,~,~] = andrieu_image_from_points2(data_a,data_e,data_n,data_yb);
[I_block_xyz,~,~] = andrieu_image_from_points(data_a,data_e,data_n,data_xyzb);
I_block_xyz = fliplr(I_block_xyz);

