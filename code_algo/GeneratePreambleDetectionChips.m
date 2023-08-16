function rval = GeneratePreambleDetectionChips(pchip, L)
%% 
% oversampling code are composed of {1,-1,0}
% we need to change all '0's to its preceding non-zero value
rval = pchip;
lastIdx = 1;
for k = 2:length(pchip)
    if pchip(k) == 0
        if k - lastIdx < L
            rval(k) = pchip(lastIdx);
        end
    else
        lastIdx = k;
    end
end

%%

end