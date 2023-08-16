function txOut = sim_mmo_tx(params)
%% model parameters
try isrepeat = params.isrepeat; catch, isrepeat = false; end

%% txIn parameters
try nMo = params.nMo; catch, nMo = 1; end
try T = params.T; catch, T = 1e-1; end
try nDegree = params.nDegree; catch, nDegree = 5; end
try nTx = params.nTx; catch, nTx = 2; end
try nBit = params.nBit; catch, nBit = 100; end
try Lp = params.Lp; catch, Lp = 2; end

try txOffset = params.txOffset; catch 
    txOffset = cell(1,nTx); txOffset(:) = {0}; 
end
try moOffset = params.moOffset; catch
    moOffset = cell(nMo,nTx); moOffset(:) = {0};
end
try noiseb = params.xChannel.noiseb; assert(iscell(noiseb)); catch 
    noiseb = cell(nMo,1); noiseb(:) = {0}; 
end
try noisen = params.xChannel.noisen; assert(iscell(noisen)); catch 
    noisen = cell(nMo,1); noisen(:) = {0}; 
end
try noisep = params.xChannel.noisep; assert(iscell(noisep)); catch 
    noisep = cell(nMo,1); noisep(:) = {0}; 
end
try xCIR = params.xChannel.xCIR; assert(iscell(xCIR)); catch
    xCIR = cell(nMo,nTx); 
    for i = 1:nTx
        xCIR(:,i) = sim_mc_cir3(struct('T', T, 'nMo', nMo));
    end
end
try xChip = params.xChip; catch
    % for one Tx, each Mo uses different code
    % each Tx uses different code combination
    xChip = cell(nMo,nTx);
    xChipAll = get_gold_code2(0, nDegree);
%     xChipAll = [xChipAll; -xChipAll];
    ind = randperm(nchoosek(size(xChipAll,2),nMo)*factorial(nMo), nTx) - 1;
    for i = 1:nTx
        indtemp = ind(i);
        indrem = 1:size(xChipAll,2);
        for j = 1:nMo
            indthis = mod(indtemp,length(indrem))+1;
            xChip(j,i) = {xChipAll(:,indrem(indthis))};
            indtemp = floor(indtemp/length(indrem));
            indrem(indrem==indthis) = [];
        end
    end
end
try xBit = params.xBit; catch
    if isrepeat
        xBit = cell(1,nTx);
    else
        xBit = cell(nMo,nTx);
    end
end

%%
assert(sum(sum(cell2mat(moOffset)==0)==0)==0);

%%
yRx = cell(nMo,1);
pChip = cell(nMo,nTx);
yTx = cell(nMo,nTx);

SINR = cell(nMo,nTx);
for j = 1:nMo
    xChannel2 = struct( ...
        'noiseb', noiseb{j}, ...
        'noisen', noisen{j}, ...
        'noisep', noisep{j}, ...
        'xCIR', {xCIR(j,:)});
    txOffsetMo = txOffset;
    for i = 1:nTx
        txOffsetMo{1,i} = txOffset{1,i} + moOffset{j,i};
    end
    
    txIn2 = struct( ...
        'T', T, ...
        'nDegree', nDegree, ...
        'nTx', nTx, ...
        'nBit', nBit, ...
        'Lp', Lp, ...
        'txOffset', {txOffsetMo}, ...
        'xChannel', {xChannel2}, ...
        'xChip', {xChip(j,:)});

    if isrepeat
        assert(size(xBit,1) == 1);
        if ~isempty(xBit{1,1})
            txIn2.xBit = xBit;
        end
    else
        if ~isempty(xBit{j,1})
            txIn2.xBit = xBit(j,:);
        end
    end
    
    txOut2 = sim_tx(txIn2);
    
    yRx(j) = {txOut2.yRx};
    if ~isrepeat || j == 1
        xBit(j,:) = txOut2.xBit;
    end
    pChip(j,:) = txOut2.pChip;
    yTx(j,:) = txOut2.yTx;
    SINR(j,:) = txOut2.SINR;
    
    if j == 1
        txOffset2 = txOut2.txOffset;
        for i = 1:nTx
            txOffset2{1,i} = txOut2.txOffset{1,i} - moOffset{j,i};
        end
    else
        assert(isequal(cell2mat(txOffset2), cell2mat(txOut2.txOffset) - cell2mat(moOffset(j,:))));
    end
end

%% return txOut
xChannel = struct( ...
    'noiseb', {noiseb}, ...
    'noisen', {noisen}, ...
    'noisep', {noisep}, ...
    'xCIR', {xCIR});
txOut = struct( ...
    'yRx', {yRx}, ...
    'txOffset', {txOffset2}, ...
    'moOffset', {moOffset}, ...
    'xChannel', xChannel, ...
    'xChip', {xChip}, ...
    'xBit', {xBit}, ...
    'pChip', {pChip}, ...
    'yTx', {yTx}, ...
    'SINR', {SINR});
end