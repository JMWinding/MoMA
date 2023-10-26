function rxOut = decode_mmo_noncoherent_MMoNTx(params)
% channel estimation in PD phase for two consecutive blocks
% and compare to decide true/false
try debug_pd = params.debug_pd; catch, debug_pd = false; end
try debug_ce = params.debug_ce; catch, debug_ce = false; end

%% rxIn parameters
try
    yRx = params.yRx;
    xChannel = params.xChannel;
    xChip = params.xChip;
    xBit = params.xBit;
    Lp = params.Lp;
    pChip = params.pChip;
    txOffset = params.txOffset;
catch
    error("missing fields");
end

try isrepeat = params.isrepeat; catch, isrepeat = false; end
try noisemodel = params.noisemodel; catch, noisemodel = "pois"; end
try algoPD = params.algoPD; catch, algoPD = "sc"; end
try algoCE = params.algoCE; catch, algoCE = "ls"; end
try nTracks = params.nTracks; catch, nTracks = 2^8; end
try sync = params.sync; catch, sync = -1; end
try mode = params.mode; catch, mode = "dc"; end
try code = params.code; catch, code = "goldman"; end

%% related variables
nMo = size(xChip,1);
nTx = size(xChip,2);
nChip = length(xChip{1});
nBit = length(xBit{1});

try moOffset = params.moOffset; catch
    if nMo == 1
        moOffset = num2cell(zeros(1,nTx));
    else
        error("missing fields");
    end
end
moOffsetMax = max(cell2mat(moOffset),[],"all");

debug_pd = debug_pd & ~isequal(algoPD,"gt");

%%
yr = yRx;
ylen = inf;
for j = 1:nMo
    ylen = min(ylen, length(yr{j}));
end

%% sliding window
plen = nChip * Lp;

nChip = length(GenerateCodeChips(xChip{1},code));
swAdv = plen/2;
swSize = 3*plen + swAdv + moOffsetMax;
swStart = 1;
swEnd = swStart + swSize-1;

% packet detection variables
lags = inf(1,nTx);
lags_temp = inf(1,nTx);
pkOffset = zeros(1,nTx);
for i = 1:nTx
    pkLoc = zeros(nMo,1);
    for j = 1:nMo
        [~, pkLoc(j)] = max(xChannel.xCIR{j,i});
    end
    pkOffset(i) = round(median(pkLoc))-1;
end
% channel estimation variables
try
    hPre = params.hPre;
    hPost = params.hPost;
catch
%     hPre = ceil(0.6/T); hPost = ceil(1.2/T);
    hPre = 7; hPost = 10;
end
chan = struct("hp", {cell(nMo,nTx)}, ...
    "hpre", {cell(nMo,nTx)}, ...
    "nb", {cell(nMo,1)}, ...
    "nn", {cell(nMo,1)}, ...
    "np", {cell(nMo,1)});
chan.hpre(:) = {-1};
switch algoCE
    case "gt"
        for i = 1:nTx
            chan.hp(:,i) = xChannel.xCIR(:,i);
            chan.hpre(:,i) = {pkOffset(i)};
        end
        chan.nb = xChannel.noiseb;
        chan.nn = xChannel.noisen;
        chan.np = xChannel.noisep;
    case "af0"
        % "0" assumes known CIR estimation (possibly from previos packet)
        for i = 1:nTx
            for j = 1:nMo
                chan.hp(j,i) = {zeros(hPre+1+hPost,1)};
                for k = 1:length(xChannel.xCIR{j,i})
                    if k-pkOffset(i)+hPre < 1 || k-pkOffset(i) > hPost+1
                        continue;
                    end
                    chan.hp{j,i}(k-pkOffset(i)+hPre) = xChannel.xCIR{j,i}(k);
                end
            end
            chan.hpre(:,i) = {hPre};
        end
        chan.nb = xChannel.noiseb;
        chan.nn = xChannel.noisen;
        chan.np = xChannel.noisep;
        % remove "0" to perform normal channel estimation
        algoCE = "af";
end
chan_gt = chan;
% decoding variables
dBit = xBit;
for j = 1:size(dBit,1)
    for i = 1:nTx
        dBit{j,i}(:) = 0;
    end
end

% lag_ce, lags after this value need channel estimation
% lag_ce = hPre - swAdv;
lag_ce = -nChip*nBit;

if debug_pd
    labels_pd = zeros(0,4);
    errors_pd = nan(0,2);
    metrics_pd = zeros(0,2,nMo,nTx);
end

%%
lags_gt = cell2mat(txOffset) + cell2mat(moOffset) + pkOffset;
for j = 1:nMo
for i = 1:nTx
    temp = yRx{j,1}(lags_gt(j,i)+plen+(0:nBit*nChip-1));
    temp = reshape(temp, nChip, []).' * (xChip{j,i}>0);
    thrd = (min(temp(xBit{j,i}==1)) + max(temp(xBit{j,i}==-1))) / 2;
    dBit{j,i} = (xBit{j,i} >= thrd) * 2 - 1;
end
end

%% return
PDOff = lags_gt - cell2mat(txOffset) - pkOffset;
rxOut.PDOff = PDOff;

rxOut.chan = chan;

rxOut.dBit = dBit;
BER = cell(size(xBit));
for j = 1:size(xBit,1)
    for i = 1:nTx
        BER{j,i} = mean(xBit{j,i}~=dBit{j,i}(1:length(xBit{j,i})));
    end
end
rxOut.BER = BER;

if debug_pd
    rxOut.debug_pd = struct( ...
        "labels_pd", labels_pd, ...
        "errors_pd", errors_pd, ...
        "metrics_pd", metrics_pd);
end

end

%% secondary functions
function rval = EvaluateEstimate(chan, esti, lags)
nMo = size(chan.hp,1);
nTx = size(chan.hp,2);

assert(length(lags) == nTx);

%% Tx score (less than 1)
score = 1;
rcorr = nan(nMo,nTx);
for i = 1:nTx
    if isinf(lags(i)), continue; end
    
    % compute correlation for each molecule signal
    xcorrMo = cell(nMo,4);
    for j = 1:nMo
        sig1 = chan.hp{j,i};
        sig2 = esti.hp{j,i};
        power = 1 - sum(abs(sig1 - sig2)) / sqrt(sum(sig1)*sum(sig2));
        scorr = dot(sig1,sig2) / (norm(sig1)*norm(sig2));

        xcorrMo(j,2) = {power};
        xcorrMo(j,3) = {scorr};
    end

    % sum over all molecule correlation
    xcorrMo = cell2mat(xcorrMo(:,3).');
    xcorrTx = min(xcorrMo,[],2);

    % find the max correlation
    scoreTx = xcorrTx;
    rcorr(:,i) = xcorrMo.';
    score = min(score, scoreTx);
end

%% Tx power ratio
ratio = nan(nMo,nTx);
for i = 1:nTx
    if isinf(lags(i)), continue; end
    for j = 1:nMo
        ratio(j,i) = sum(esti.hp{j,i}) / sum(chan.hp{j,i});
    end
end

%% noise score
% for j = 1:nMo
%     if chan.nb{j} > 0
%         score = min(score, -abs(esti.nb{j}/chan.nb{j}-1));
%     end
%     if chan.nn{j} > 0
%         score = min(score, -abs(esti.nn{j}/chan.nn{j}-1));
%     end
%     if chan.np{j} > 0
%         score = min(score, -abs(esti.np{j}/chan.np{j}-1));
%     end
% end

%% return
rval.score = score;
rval.corr = rcorr;
rval.ratio = ratio;
end

function ret = relu(input)
ret = input .* (input > 0);
end