function [ minindex_filtered, minindex_bad ] = radial_subsample3( data_xy, data_a, data_e, data_z,...
    abins,...
    a_blocked,xy_blocked)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

verbose =0;
thetar = (360/abins)*pi/180;
thetad = (360/abins);

minindex= data_xy*0;

r=1.5;
rincrement = r;
while r<max(data_xy)
    r = r+r*thetar;
    rincrement = [rincrement,r];
end

if verbose; figure; hold on; end;


for a = min(data_a):thetad:max(data_a)
    aindex = logical(data_a>a&data_a<a+thetad);
    for i=1:numel(rincrement)-1;
        rindex = logical(data_xy>rincrement(i)&data_xy<rincrement(i+1));
        index =  aindex.*rindex;
        if verbose; scatter(5*rand(1)*rincrement(i), sum(index), 'r','filled'); end;
        if numel(index>0)
        temp = data_z;
        temp(~index)=inf;
        if min(temp)~=inf
        minimum = logical(temp==min(temp));
        minindex = minindex+minimum;
        end
        end
    end
end

minindex = logical(minindex);
minindex_filtered = false(numel(minindex),1);
minindex_bad = false(numel(minindex),1);

% Remove points within lidar plate shadow
Imin = find(minindex);
for i = 1:numel(Imin);
    a_temp = data_a(Imin(i));
    xy_temp = data_xy(Imin(i));
    [~,idx] = min(abs(a_blocked - a_temp));
    ab_temp = a_blocked(idx);
    xyb_temp = xy_blocked(idx);
    if xy_temp > xyb_temp;
        minindex_filtered(Imin(i)) = true;
    %    fprintf('\nPoint (%f,%f) was accepted beacuse limit is (%f,%f)\n',...
    %        a_temp, xy_temp, ab_temp, xyb_temp)
    else
     %  fprintf('\nPoint (%f,%f) was rejected beacuse limit is (%f,%f)\n',...
     %       a_temp, xy_temp, ab_temp, xyb_temp)
            minindex_bad(Imin(i)) = true;
    end
end

   

end

