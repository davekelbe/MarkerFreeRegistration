function [x] = decompose_translation_ty(t)
    if isempty(t);
        x = nan;
        return
    end
	x = t(2);
end