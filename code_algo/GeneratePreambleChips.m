function rval = GeneratePreambleChips(xchip, Lp)
%% vector check
if size(xchip,1) == 1
    xchip = xchip.';
elseif size(xchip,2) ~= 1
    error("input is not a vector");
end

%% 
% type 1: Lp bits (modulted by) code, -1 x -1 -> 1
% xbit = GeneratePreambleBits(Lp);
% rval = reshape(xchip.*xbit.',[],1);

% type 2: code (modulated by) Lp ones, -1 x -1 -> 1
% rval = reshape(ones(Lp,1).*xchip.',[],1);

% type 3: code (modulated by) PPM
% ppm = zeros(Lp,1);
% ppm(1) = 1;
% ppm(Lp/2+1) = -1;
% rval = reshape(ppm.*xchip.',[],1);

% type 4: for gold + manchester code
% TESTBED
ppm = zeros(2*Lp,1);
ppm(1:Lp) = 1;
ppm(Lp+1:end) = -1;
rval = reshape(ppm.*xchip(1:2:end).',[],1);
% rval = [rval; zeros(length(xchip)*2,1)];
end