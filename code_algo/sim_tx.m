function txOut = sim_tx(params)
%% input parameters
try noiseb = params.xChannel.noiseb; catch, noiseb = 0; end
try noisen = params.xChannel.noisen; catch, noisen = 0; end
try noisep = params.xChannel.noisep; catch, noisep = 0; end
try T = params.T; catch, T = 100; end
try nDegree = params.nDegree; catch, nDegree = 5; end
try nTx = params.nTx; catch, nTx = 2; end
try nBit = params.nBit; catch, nBit = 1000; end
try Lp = params.Lp; catch, Lp = 4; end

try txOffset = params.txOffset; catch 
    txOffset = cell(1,nTx); txOffset(:) = {0}; 
end
try xCIR = params.xChannel.xCIR; catch 
    xCIR = cell(1,nTx); xCIR(:) = sim_mc_cir3(struct("T", T));
end
try xChip = params.xChip; catch 
    xChip = get_gold_code(nTx, nDegree);
    xChip = mat2cell(xChip, size(xChip,1), ones(size(xChip,2),1));
end
try xBit = params.xBit; catch
    xBit = randi([0,1], [nBit,nTx]);
    xBit = 2 * xBit - 1;
    xBit = mat2cell(xBit, nBit, ones(nTx,1));
end
try pChip = params.pChip; catch
    pChip = cell(1,nTx);
    for i = 1:nTx
        pChip(i) = {GeneratePreambleChips(xChip{i}, Lp)};
    end
end

%% useful variables
nChip = length(xChip{1});
plen = nChip * Lp;

%% simulate tx
xTx = cell(1,nTx);
for i = 1:nTx
    xTx{i} = ToPos([pChip{i}; GenerateDataChips(xBit{i},xChip{i})]);
end

% conv and add pre packet samples

nPad = 5 * plen;
for i = 1:nTx
    txOffset{i} = txOffset{i} + nPad;
end

yRxLen = 0;
yTx = cell(1,nTx);
for i = 1:nTx
    yTx{i} = [zeros(txOffset{i},1); conv(xTx{i}, xCIR{i})];
    if length(yTx{i}) > yRxLen
        yRxLen = length(yTx{i});
    end
end

yRx = zeros(yRxLen,1);
for i = 1:nTx
    yRx(1:length(yTx{i})) = yRx(1:length(yTx{i})) + yTx{i};
end
yRx = [yRx; zeros(nPad,1)];

% add noise
% background noise
yRx = yRx + noiseb;
% signal-dependent normal noise
yRx = yRx + randn(size(yRx)) .* noisep .* sqrt(yRx);
yRx = max(yRx, 0);
% AWGN
yRx = yRx + randn(size(yRx)) .* noisen;

% compute SNR
SINR = cell(1,nTx);
pSig = rssq(yRx(:));
for i = 1:nTx
    yIN = yRx;
    yIN(1:length(yTx{i})) = yIN(1:length(yTx{i})) - yTx{i};
    pIN = rssq(yIN(:));
    SINR{i} = mag2db(pSig/pIN);
end

%% return txOut
xChannel = struct("xCIR", xCIR, "noiseb", noiseb, ...
    "noisen", noisen, "noisep", noisep);
txOut = struct( ...
    "T", T, ...
    "yRx", {yRx}, ...
    "txOffset", {txOffset}, ...
    "xChannel", {xChannel}, ...
    "xChip", {xChip}, ...
    "xBit", {xBit}, ...
    "Lp", Lp, ...
    "pChip", {pChip}, ...
    "yTx", {yTx}, ...
    "SINR", {SINR});

end