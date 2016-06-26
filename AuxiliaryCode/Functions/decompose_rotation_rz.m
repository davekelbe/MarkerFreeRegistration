function [z] = decompose_rotation_rz(R)
    if isempty(R);
        z = nan;
        return
    end
	z = atan2(R(2,1), R(1,1));
end