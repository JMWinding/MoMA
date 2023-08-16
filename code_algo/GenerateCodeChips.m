function rval = GenerateCodeChips(xchip,code)
%%
idx = find(xchip~=0,2);
if length(idx) == 2
    oversampling = idx(2) - idx(1);
else
    oversampling = length(xchip);
end

%%
switch code
    case {"goldman", "gold", "ooc"}
        rval = xchip;
    case {"goldman0", "ooc0"}
        rval = xchip>0;
    case {"man"}
        rval = zeros(length(xchip),1);
        rval(1:oversampling:end/2) = 1;
        rval(end/2+1:oversampling:end) = -1;
    case {"plain"}
        rval = zeros(length(xchip),1);
        rval(1:oversampling:end/2) = 1;
    case {"plain0"}
        rval = zeros(2*oversampling,1);
        rval(1) = 1;
    otherwise
        error("code not defined");
end
end