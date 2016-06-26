function [y] = decompose_rotation_ry(R)
    if isempty(R);
        y = nan;
        return
    end
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
end