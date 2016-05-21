function [I, I1, I2 ] = andrieu_image_from_points(data_a,data_e, data_n, data_i )
%ANDRIEU_IMAGE_FROM_POINTS Make and Andrieu Projection Image from point cloud
%
%   Example:
%   [I,I1,I2] = andrieu_image_from_points(data_a,data_e,data_n,data_i);
%   Andrieu image will be color coded according to the intensity 
%
%
%   Outputs:
%   I  - image with all returns plotted 
%   I1 - image with just first returns
%   I2 - image with just second returns 
%
%   (C) David Kelbe, Rochester Institute of Technology 

a_unique = sort(unique(data_a));
e_unique = sort(unique(data_e));

n_i = numel(data_a);
n_c = numel(a_unique);
n_r = numel(e_unique);

I = zeros( n_r, n_c);

data_r = 4*(data_e + 45)+1;

c_unique = (1:n_c)';
data_c = zeros(n_i,1);
for c = 1:n_c;
    ixc = (data_a==a_unique(c));
    data_c(ixc) = c_unique(c);
end

for i = 1:n_i
   c = data_c(i);
   r = data_r(i);
   I(r,c) = data_i(i);
end

I = flipud(I);
%figure; imagesc(I); colormap jet; axis image;
if nargout>=2;
I1 = zeros(n_r, n_c);
ix_1=(data_n==1);
n_i1 = sum(ix_1);
data_c1 = data_c(ix_1);
data_r1 = data_r(ix_1);
data_i1 = data_i(ix_1);
for i = 1:n_i1
   c = data_c1(i);
   r = data_r1(i);
   I1(r,c) = data_i1(i);
end
I1 = flipud(I1);
%figure; imagesc(I1); colormap jet; axis image;
end
if nargout==3;
I2 = zeros(n_r, n_c);
ix_2=(data_n==2);
n_i2 = sum(ix_2);
data_c2 = data_c(ix_2);
data_r2 = data_r(ix_2);
data_i2 = data_i(ix_2);
for i = 1:n_i2
   c = data_c2(i);
   r = data_r2(i);
   I2(r,c) = data_i2(i);
end
I2 = flipud(I2);
%figure; imagesc(I2); colormap jet; axis image;
end

end

