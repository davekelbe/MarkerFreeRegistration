function [x] = decompose_translation_tx(t)
    if isempty(t);
        x = nan;
        return
    end
	x = t(1);
end