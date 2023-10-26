function rval = GenerateDataChips(xbit, xchip)
%% input check
if size(xbit,1) == 1
    xbit = xbit.';
elseif size(xbit,2) ~= 1
    error("input is not a vector");
end
if size(xchip,1) == 1
    xchip = xchip.';
elseif size(xchip,2) ~= 1
    error("input is not a vector");
end

%% -1 x -1 -> 1
rval = reshape(xchip.*xbit.',[],1);
end