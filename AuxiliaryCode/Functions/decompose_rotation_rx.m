function [x] = decompose_rotation_rx(R)
    if isempty(R);
        x = nan;
        return
    end
	x = atan2(R(3,2), R(3,3));
end