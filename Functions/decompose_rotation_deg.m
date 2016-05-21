function [x,y,z] = decompose_rotation_deg(R)
    x = [];
    y = [];
    z = [];
    if isempty(R);
        return
    end
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
    x = rad2deg(x);
    y = rad2deg(y);
    z = rad2deg(z);
end