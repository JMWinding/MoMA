function rval = GenerateOversampleChips(chip, oversampling)
assert(size(chip,2) == 1);
rval = reshape([chip.'; zeros(oversampling-1,length(chip))],[],1);
end