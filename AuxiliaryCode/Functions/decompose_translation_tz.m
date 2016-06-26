function [x] = decompose_translation_tz(t)
    if isempty(t);
        x = nan;
        return
    end
	x = t(3);
end