function  [data_a2, data_e]  = az_to_360( data_a, data_e, data_z )
%180AZTO360 Convert 180 azimuth range to 360 based on elevation
%   0 begins at -x axis and goes counter clockwise
%
% (c) David Kelbe, Rochester Institute of Technology

index = logical(data_e>90);
data_a(index) = data_a(index) + 180;
data_e(index) = abs(data_e(index) - 180);
index2 = logical(data_z<0&index);
data_e(index2) = -data_e(index2);
data_a2 = (360-data_a)-270;
ix_negative = (data_a2<0);
data_a2(ix_negative) = data_a2(ix_negative)+360;

end

