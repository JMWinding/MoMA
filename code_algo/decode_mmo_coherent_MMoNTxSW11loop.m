function rxOut = decode_mmo_coherent_MMoNTxSW11loop(params)
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
    error('missing fields');
end

try isrepeat = params.isrepeat; catch, isrepeat = false; end
try noisemodel = params.noisemodel; catch, noisemodel = 'pois'; end
try algoPD = params.algoPD; catch, algoPD = 'sc'; end
try algoCE = params.algoCE; catch, algoCE = 'ls'; end
try nTracks = params.nTracks; catch, nTracks = 2^8; end
try sync = params.sync; catch, sync = -1; end
try mode = params.mode; catch, mode = 'dc'; end
try code = params.code; catch, code = 'goldman'; end
try sameMo = params.sameMo; catch, sameMo = false; end

%% related variables
nMo = size(xChip,1);
nTx = size(xChip,2);
nChip = length(xChip{1});
nBit = length(xBit{1});

try moOffset = params.moOffset; catch
    if nMo == 1
        moOffset = num2cell(zeros(1,nTx));
    else
        error('missing fields');
    end
end
moOffsetMax = max(cell2mat(moOffset),[],'all');

try
    weights_ce = params.weights_ce;
catch
    weights_ce = struct( ...
        "afloss", "mean", ...
        "pos", 1, ...
        "posy", 0, ...
        "simTx", 1, ...
        "simMo", 0, ...
        "smth", 0, ...
        "cntr", 1);
end
try
    weights_pd = params.weights_pd;
catch
    weights_pd = weights_ce;
end
try thrd_pd_corr = params.thrd_pd.corr; catch
    thrd_pd_corr = 0.8;
end
try thrd_pd_ratio = params.thrd_pd.ratio; catch
    thrd_pd_ratio = 0.5;
end

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
    hPre = ceil(0.6/T); hPost = ceil(1.2/T);
end
chan = struct('hp', {cell(nMo,nTx)}, ...
    'hpre', {cell(nMo,nTx)}, ...
    'nb', {cell(nMo,1)}, ...
    'nn', {cell(nMo,1)}, ...
    'np', {cell(nMo,1)});
chan.hpre(:) = {-1};
switch algoCE
    case 'gt'
        for i = 1:nTx
            chan.hp(:,i) = xChannel.xCIR(:,i);
            chan.hpre(:,i) = {pkOffset(i)};
        end
        chan.nb = xChannel.noiseb;
        chan.nn = xChannel.noisen;
        chan.np = xChannel.noisep;
    case 'af0'
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
        algoCE = 'af';
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
    corr_pd = zeros(0,nMo+1,nTx);
    ratio_pd = zeros(0,nMo,nTx);
    ratio2_pd = zeros(0,nMo,nTx);
end

%%
isFirstPacket = true;
while true
    %% clip windowed signals
    yr2 = cell(nMo,1);
    for j = 1:nMo
        yr2{j} = yr{j}(swStart:swEnd);
    end

    %% range of different "lags"
    % lags \in (-inf, plen), previously confirmed
    % lags_new \in [plen, plen + 2*swAdv), currently confirmed

    % lags_old \in (-inf, lag_ce), no need for channel estimation
    lags_old = lags; lags_old(lags>=lag_ce) = inf;
    % lags_ce \in [lag_ce, plen + swAdv), need channel estimation

    % lags_gt \in [plen, plen + 2*swAdv), ground truth
    lags_gt = inf(1,nTx);
    for i = 1:nTx
        temp2 = txOffset{i} - (swStart-1) + pkOffset(i);
        if temp2 >= plen && temp2 < plen + 2*swAdv
            lags_gt(i) = temp2;
        end
    end

    % lags_min/lags_max, rough arriving range
    pdDo = true;
    if sync >= 0
        if ~exist('lags_min', 'var')
            lags_min = zeros(1,nTx);
            for i = 1:nTx
                lags_min(i) = txOffset{i} - (swStart-1) + pkOffset(i) - sync;
            end
            lags_max = lags_min + 2*sync;
        else
            lags_min = lags_min - swAdv;
            lags_max = lags_max - swAdv;
        end

        if sum(isinf(lags)) == 0 ...
                || min(lags_min(isinf(lags))) >= plen + 2*swAdv ...
                || max(lags_max(isinf(lags))) < plen
            pdDo = false;
        end
    end
    
    %% packet detection + temporal channel estimation
    chan_backup_gt = chan;
    for algoPD2 = [algoPD, "gt"]
        chan = chan_backup_gt;
        lags_new = lags_temp;

    while true
        %% decode (old Tx + confirmed new Tx)
        lags2 = lags; lags2(~isinf(lags_new)) = lags_new(~isinf(lags_new));
        if sum(~isinf(lags2)) > 0
            forwardIn = struct( ...
                'yd', {yr2}, 'xBit', {xBit}, ...
                'xChip', {xChip}, 'pChip', {pChip}, ...
                'lags', lags2, ...
                'moOffset', {moOffset}, 'hPres', {chan.hpre}, ...
                'chan', chan, 'mode', 'pd', ...
                'endIdx', min(3.5*plen+moOffsetMax,ylen-swStart+1), ...
                'nTracks', nTracks, 'code', code);
            if exist('checkpoint', 'var')
                forwardIn.checkpoint = checkpoint;
            end
            forward = DecodeSequence_ViterbiForward4(forwardIn);
        
            backwardIn = struct( ...
                'viterbi', forward.viterbi, ...
                'dBit', {dBit}, 'lags', lags2);
            backward = DecodeSequence_ViterbiBackward(backwardIn);
            dBit2 = backward.dBit;
        else
            dBit2 = dBit;
        end

        %% reconstruct decoded signal
        yo2 = yr2;
        for j = 1:nMo
            yo2{j}(:) = chan.nb{j};
        end
        for i = 1:nTx
            if isinf(lags(i)) && isinf(lags_new(i)), continue; end
            for j = 1:nMo
                dbitTemp = dBit2{j,i};
                xchip = GenerateCodeChips(xChip{j,i},code);
                dchipTemp = ToPos([pChip{j,i}; GenerateDataChips(dbitTemp,xchip)]);
                if ~isinf(lags(i))
                    delta = lags(i) + moOffset{j,i} - chan.hpre{j,i};
                else
                    delta = lags_new(i) + moOffset{j,i} - chan.hpre{j,i};
                end
                yo2{j} = RebuildKnownPacket(yo2{j},dchipTemp,chan.hp{j,i},delta);
            end
        end

        %% compute residual
        resRange = plen + swAdv + (1:1.5*plen);
        respTemp = zeros(length(resRange),nMo);
        for j = 1:nMo
            respTemp(:,j) = (yr2{j}(resRange) - yo2{j}(resRange)) ...
                ./ yr2{j}(resRange);
        end
        respTemp = reshape(respTemp, plen/2, [], nMo);
    
        % decide if potential new packets
        if mean(respTemp(:,2,:), 'all') < 0 ...
                || mean(respTemp(:,3,:), 'all') < 0.1
            pdDo = false;
        end

        %% packet detection
        pdIn = struct( ...
            'isrepeat', isrepeat, ...
            'yp', {yr2}, 'xChip', {xChip}, 'moOffset', {moOffset}, ...
            'xBit', {xBit}, 'dBit', {dBit2}, ...
            'Lp', Lp, 'pChip', {pChip}, ...
            'hPre', hPre, 'hPost', hPost, 'hp', {chan.hp}, ...
            'lags', lags2, 'hPres', {chan.hpre}, ...
            'algo', algoPD2, 'algoCE', algoCE, 'code', code);

        % other pdIn
        if sync >=0
            pdIn.lags_min = lags_min;
            pdIn.lags_max = lags_max;
        end
        pdIn.lags_gt = lags_gt;
        
        if isequal(algoPD2,"gt")
%             if ~debug_pd
%                 lags_new = lags_gt;
%                 break;
%             else
                lags_all = lags;
                lags_all(~isinf(lags_gt)) = lags_gt(~isinf(lags_gt));
%             end
        elseif pdDo
            pd = PreambleDetectionSW8(pdIn);
            lags_all = pd.lags;
            chan.hpre = pd.hPres;
            chan.hp = pd.hp;
        else
            lags_all = lags;
        end

        %% any potential new packet?
        lags_newp = lags_all;
        lags_newp(~isinf(lags)) = inf;
        lags_newp(~isinf(lags_new)) = inf;
        if sum(~isinf(lags_newp)) == 0
            % no new packet detected, so exit the while
            break;
        end

        %% check each potential new packet
        new_tx_confirmed = false;
        [~, idx_newp] = sort(lags_newp, 'ascend');
        for idx_newp2 = idx_newp
            if isinf(lags_newp(idx_newp2))
                break;
            end

            % estimate the channel
            lags_ce2 = lags; lags_ce2(lags<lag_ce) = inf;
            lags_ce2(~isinf(lags_new)) = lags_new(~isinf(lags_new));
            lags_ce2(idx_newp2) = lags_newp(idx_newp2);

            if isequal(algoCE,'gt')
                moOffset2 = cell2mat(moOffset)-cell2mat(chan.hpre);
                moOffset2(:,~isinf(lags_ce2)) = ...
                    cell2mat(moOffset(:,~isinf(lags_ce2))) - hPre;
                moOffset2 = num2cell(moOffset2);
            else
                moOffset2 = num2cell(cell2mat(moOffset)-cell2mat(chan.hpre));
            end

            if isFirstPacket
                nbPD = cell(nMo,1);
                for j = 1:nMo
                    nbPD{j} = mean(yr2{j}(1:min(lags_ce2)));
                end
            else
                nbPD = chan.nb;
            end

            % SIC
            dBit3_old = dBit2;
            lags3 = lags; lags3(~isinf(lags_ce2)) = lags_ce2(~isinf(lags_ce2));

            for niter = 0:nTx+1
                if niter == nTx+1
                    break;
                end

                %% estimate channel
                endIdxCE = max(lags_ce2(~isinf(lags_ce2))) + plen - hPre;
                estimateIn2 = struct( ...
                    'noisemodel', noisemodel, 'isrepeat', isrepeat, ...
                    'yr', {yr2}, 'lags_ce', lags_ce2, 'lags_old', lags_old, ...
                    'moOffset', {moOffset2}, 'xChip', {xChip}, ...
                    'dBit', {dBit3_old}, 'pChip', {pChip}, ...
                    'hPre', hPre, 'hPost', hPost, 'hp', {chan.hp}, 'nb', {nbPD}, ...
                    'algo', 'af', 'weights', weights_pd, 'code', code, ...
                    'sameMo', sameMo, ...
                    'endIdx', endIdxCE);
                estimate2 = EstimateChannelMMoSW(estimateIn2);

                %% decode
                forwardIn = struct( ...
                    'yd', {yr2}, 'xBit', {xBit}, ...
                    'xChip', {xChip}, 'pChip', {pChip}, ...
                    'lags', lags3, ...
                    'moOffset', {moOffset}, 'hPres', {chan.hpre}, ...
                    'chan', estimate2, 'mode', 'pd', ...
                    'endIdx', max(lags_ce2(~isinf(lags_ce2))) + plen, ...
                    'nTracks', nTracks, 'code', code);
                if exist('checkpoint', 'var')
                    forwardIn.checkpoint = checkpoint;
                end
                forward3 = DecodeSequence_ViterbiForward4(forwardIn);
            
                backwardIn = struct( ...
                    'viterbi', forward3.viterbi, ...
                    'dBit', {dBit}, 'lags', lags3);
                backward3 = DecodeSequence_ViterbiBackward(backwardIn);
                dBit3 = backward3.dBit;

                %% compare
                if niter == 0
                    dBit3_old = dBit3;
                elseif isequal(dBit3, dBit3_old)
                    break;
                else
                    dBit3_old = dBit3;
                end
            end
            
            if niter == nTx+1
                % iteration does not converge in expected steps,
                % so we assume false positive detection
                continue;
            else
                dBit2 = dBit3;
            end

            % debug: compare channel estimation
            while debug_ce
                f1 = figure("Units","Normalized","Outerposition",[0 0 1 1], ...
                    "Name", "PD");
                nTx2 = sum(~isinf(lags_ce2));

                %%
                i2 = 1;
                for i = 1:nTx
                    if isinf(lags_ce2(i)), continue; end
                    for j = 1:nMo
                        temp1 = floor((1 - lags_ce2(i) - plen) / nChip);
                        temp2 = floor((endIdxCE - lags_ce2(i) - plen) / nChip);

                        % bits
                        subplot(nTx2+2,2*nMo,(i2-1)*2*nMo+j*2-1);
                        hold on; box on; grid on;
                        stem(xBit{j,i});
                        stem(dBit2{j,i});
                        if temp1 > 0 && temp1 <= nBit
                            stem(temp1,xBit{j,i}(temp1),'k','LineWidth',2);
                        end
                        if temp2 > 0 && temp2 <= nBit
                            stem(temp2,xBit{j,i}(temp2),'--k','LineWidth',2);
                        end
    
                        % CIR
                        subplot(nTx2+2,2*nMo,(i2-1)*2*nMo+j*2);
                        hold on; box on; grid on;
                        title("Tx"+num2str(i)+", Mo"+num2str(j)+", pd");
                        plot((1:length(chan_gt.hp{j,i}))-chan_gt.hpre{j,i}-1, chan_gt.hp{j,i}, '-o');
                        plot((1:length(chan.hp{j,i}))-chan.hpre{j,i}-1, chan.hp{j,i}, '-+');
                        if isinf(lags_ce2(i)), continue; end
                        plot((1:length(estimate2.hp{j,i}))-hPre-1, estimate2.hp{j,i}, 'LineWidth', 2);
                    end
                    i2 = i2 + 1;
                end

                %%
                yo0 = yr2; % confirmed packet, xbit, chan_gt
                yo = yr2; % confirmed packet, xbit, chan
                yo2 = yr2; % confirmed packet, dbit2, chan
                yo3 = yr2; % confirmed + temp packet, dbit2, estimate
                for j = 1:nMo
                    yo0{j}(:) = chan_gt.nb{j};
                    yo{j}(:) = chan.nb{j};
                    yo2{j}(:) = chan.nb{j};
                    yo3{j}(:) = estimate2.nb{j};
                end
                for i = 1:nTx
                    if isinf(lags(i)) && isinf(lags_ce2(i)), continue; end
                    for j = 1:nMo
                        xbit = xBit{j,i};
                        dbit = dBit2{j,i};
                        xchip = GenerateCodeChips(xChip{j,i},code);
                        xchips = ToPos([pChip{j,i}; GenerateDataChips(xbit,xchip)]);
                        dchips = ToPos([pChip{j,i}; GenerateDataChips(dbit,xchip)]);
                        if ~isinf(lags(i))
                            delta = lags(i) + moOffset{j,i} - chan.hpre{j,i};
                        else
                            delta = lags_ce2(i) + moOffset{j,i} - chan.hpre{j,i};
                        end
                        if ~isinf(lags(i)) || ~isinf(lags_new(i))
                            yo0{j} = RebuildKnownPacket(yo0{j},xchips,chan_gt.hp{j,i},delta);
                            yo{j} = RebuildKnownPacket(yo{j},xchips,chan.hp{j,i},delta);
                            yo2{j} = RebuildKnownPacket(yo2{j},dchips,chan.hp{j,i},delta);
                        end
                        yo3{j} = RebuildKnownPacket(yo3{j},dchips,estimate2.hp{j,i},delta);
                    end
                end
        
                for j = 1:nMo
                    % total
                    subplot(nTx2+2,nMo,nTx2*nMo+j); hold on; box on; grid on;
                    plot(yr2{j});
                    plot(yo3{j});
                    plot([endIdxCE,endIdxCE],ylim,'--k','LineWidth',2);
                    xticks(0:plen/2:3.5*plen); xlim([0,3.5*plen]);
                    legend('rx signal', 'temp signal');
                    title("Mol "+string(j)+" total");

                    % subtract
                    subplot(nTx2+2,nMo,(nTx2+1)*nMo+j); hold on; box on; grid on;
                    plot(yr2{j}-yo0{j});
                    plot(yr2{j}-yo{j});
                    plot(yr2{j}-yo2{j});
                    plot([endIdxCE,endIdxCE],ylim,'--k','LineWidth',2);
                    xticks(0:plen/2:3.5*plen); xlim([0,3.5*plen]);
                    legend('gt residual', 'expected residual', 'actual residual');
                    title("Mol "+string(j)+" SIC");
                end

                %%
%                 close(f1);
                break;
            end

            % compare channel of two consecutive blocks
            estimateInHalf = struct( ...
                'noisemodel', noisemodel, 'isrepeat', isrepeat, ...
                'yr', {yr2}, 'lags_ce', lags_ce2, 'lags_old', lags_old, ...
                'moOffset', {moOffset2}, 'xChip', {xChip}, ...
                'dBit', {dBit3_old}, 'pChip', {pChip}, ...
                'hPre', hPre, 'hPost', hPost, 'hp', {estimate2.hp}, 'nb', {estimate2.nb}, ...
                'algo', 'af', 'weights', weights_pd, 'code', code, ...
                'sameMo', sameMo, ...
                'endIdx', max(lags_ce2(~isinf(lags_ce2))) + round((plen-hPre)/2));
            estimateHalf1 = EstimateChannelMMoSW(estimateInHalf);

            estimateInHalf.sttIdx = max(lags_ce2(~isinf(lags_ce2))) + round((plen-hPre)/2);
            estimateInHalf.endIdx = max(lags_ce2(~isinf(lags_ce2))) + (plen-hPre);
            estimateHalf2 = EstimateChannelMMoSW(estimateInHalf);

            % compare channel
            eval = EvaluateEstimate(estimateHalf1, estimateHalf2, lags_ce2);
            
            % DEBUG METRIC
            while debug_pd
                %
                ind_ahead = lags_gt<lags_gt(idx_newp2);
                cnt_ahead = sum(isinf(min(lags(ind_ahead),lags_new(ind_ahead))));
                ind_after = lags_gt>lags_gt(idx_newp2) & ~isinf(lags_gt);
                cnt_after = sum(isinf(min(lags(ind_after),lags_new(ind_after))));
                if abs(lags_newp(idx_newp2)-lags_gt(idx_newp2)) <= hPre
                    % 1 is a hit
                    labels_pd = [labels_pd; 1, isequal(algoPD2,"gt"), ...
                        cnt_ahead, cnt_after];
                else
                    % 0 is a wrong hit
                    labels_pd = [labels_pd; 0, isequal(algoPD2,"gt"), ...
                        cnt_ahead, cnt_after];
                end
                %
                errors_pd = [errors_pd; idx_newp2, lags_newp(idx_newp2)-lags_gt(idx_newp2)];
                %
                corr_temp = zeros([1,size(corr_pd,[2,3])]);
                corr_temp(1,:,:) = [eval.corr; eval.wcorr];
                corr_pd = [corr_pd; corr_temp];
                ratio_temp = zeros([1,size(ratio_pd,[2,3])]);
                ratio_temp(1,:,:) = eval.ratio;
                ratio_pd = [ratio_pd; ratio_temp];
                ratio2_temp = zeros([1,size(ratio2_pd,[2,3])]);
                ratio2_temp(1,:,:) = eval.ratio2;
                ratio2_pd = [ratio2_pd; ratio2_temp];

                break;
            end

            eval.ratio(eval.ratio>1) = 1./eval.ratio(eval.ratio>1);
            eval.ratio2(eval.ratio2>1) = 1./eval.ratio2(eval.ratio2>1);

            if isequal(algoPD2,"gt") ...
               || (min(eval.wcorr,[],2,"omitnan") > thrd_pd_corr ...
                   && min(mean(eval.ratio,1,"omitnan"),[],2,"omitnan") > thrd_pd_ratio)
                new_tx_confirmed = true;
                lags_new(idx_newp2) = lags_newp(idx_newp2);

                chan.hp = estimate2.hp;
                chan.nb = estimate2.nb;
                isFirstPacket = false;
                break;
            end
        end

        %% no confirmed new packet
        if ~new_tx_confirmed
            % no confirmed new packet
            if debug_pd
                if sum(isinf(lags_new(~isinf(lags_gt)))) > 0
                    % -1 is a miss
                    labels_pd = [labels_pd; -1, isequal(algoPD2,"gt"), ...
                        sum(isinf(lags_new(~isinf(lags_gt)))), nan];
                    %
                    errors_pd = [errors_pd; mean(respTemp(:,2,:), 'all'), ...
                        mean(respTemp(:,3,:), 'all')];
                    %
                    corr_temp = nan([1,size(corr_pd,[2,3])]);
                    corr_pd = [corr_pd; corr_temp];
                    ratio_temp = nan([1,size(ratio_pd,[2,3])]);
                    ratio_pd = [ratio_pd; ratio_temp];
                    ratio2_temp = nan([1,size(ratio2_pd,[2,3])]);
                    ratio2_pd = [ratio2_pd; ratio2_temp];
                end
            end
            % exit the while
            break;
        end

    end
    
    if isequal(algoPD,'gt') || ~debug_pd || isequal(lags_new,lags_gt)
        break;
    end
    end

    if ~isequal(mode,'dc') && (sum(isinf(lags)) == 0 || max(lags_gt) < 0)
        break;
    end
        
    %% final channel estimation
    lags_ce = lags; lags_ce(lags<lag_ce) = inf;
    lags_ce(lags_new<plen+swAdv) = lags_new(lags_new<plen+swAdv);
    if ~isequal(algoCE,'gt') && sum(~isinf(lags_ce)) > 0
        lags_ce_last = max(lags_ce(~isinf(lags_ce)));
        lags_ce_next = min(lags_new(lags_new>lags_ce_last));
        if numel(lags_ce_next) == 0, lags_ce_next = inf; end
        endIdxCE = max(lags_ce_last + hPost+1, ...
            min(lags_ce_next, max(plen, lags_ce_last) + plen) - hPre);
        lags_ce(lags_new<endIdxCE) = lags_new(lags_new<endIdxCE);

        moOffset2 = num2cell(cell2mat(moOffset)-cell2mat(chan.hpre));
        estimateIn = struct( ...
            'noisemodel', noisemodel, 'isrepeat', isrepeat, ...
            'yr', {yr2}, 'lags_ce', lags_ce, 'lags_old', lags_old, ...
            'moOffset', {moOffset2}, 'xChip', {xChip}, ...
            'dBit', {dBit2}, 'pChip', {pChip}, ...
            'hPre', hPre, 'hPost', hPost, 'hp', {chan.hp}, 'nb', {chan.nb}, ...
            'algo', algoCE, 'weights', weights_ce, 'code', code, ...
            'sameMo', sameMo, ...
            'endIdx', endIdxCE);

        estimateIn.hp = chan.hp;
        estimateIn.algo = algoCE;
        estimate = EstimateChannelMMoSW(estimateIn);

        % debug: compare channel estimation
        while debug_ce && sum(isinf(lags_ce))
            f2 = figure("Units","Normalized","Outerposition",[0 0 1 1], ...
                "Name", "CE");
            nTx2 = sum(~isinf(lags_ce));

            %%
            i2 = 1;
            for i = 1:nTx
                if isinf(lags_ce(i)), continue; end
                for j = 1:nMo
                    temp1 = floor((1 - lags_ce(i) - plen) / nChip);
                    temp2 = floor((endIdxCE - lags_ce(i) - plen) / nChip);
                    % bits
                    subplot(nTx2+2,2*nMo,(i2-1)*2*nMo+j*2-1);
                    hold on; box on; grid on;
                    stem(xBit{j,i});
                    stem(dBit2{j,i});
                    if temp1 > 0 && temp1 <= nBit
                        stem(temp1,xBit{j,i}(temp1),'k','LineWidth',2);
                    end
                    if temp2 > 0 && temp2 <= nBit
                        stem(temp2,xBit{j,i}(temp2),'--k','LineWidth',2);
                    end

                    % CIR
                    subplot(nTx2+2,2*nMo,(i2-1)*2*nMo+j*2);
                    hold on; box on; grid on;
                    title("Tx"+num2str(i)+", Mo"+num2str(j)+", ce");
                    plot((1:length(chan_gt.hp{j,i}))-chan_gt.hpre{j,i}-1, chan_gt.hp{j,i}, '-o');
                    plot((1:length(chan.hp{j,i}))-chan.hpre{j,i}-1, chan.hp{j,i}, '-+');
                    if isinf(lags_ce(i)), continue; end
                    plot((1:length(estimate.hp{j,i}))-hPre-1, estimate.hp{j,i}, 'LineWidth', 2);
                end
                i2 = i2 + 1;
            end

            %%
            yo = yr2;
            yo1 = yr2;
            yo2 = yr2;
            for j = 1:nMo
                yo{j}(:) = chan_gt.nb{j};
                yo1{j}(:) = chan.nb{j};
                yo2{j}(:) = estimate.nb{j};
            end
            for i = 1:nTx
                if isinf(lags(i)) && isinf(lags_ce(i)), continue; end
                for j = 1:nMo
                    xbit = xBit{j,i};
                    dbit = dBit2{j,i};
                    xchip = GenerateCodeChips(xChip{j,i},code);
                    dchip = ToPos([pChip{j,i}; GenerateDataChips(xbit,xchip)]);
                    dchipTemp = ToPos([pChip{j,i}; GenerateDataChips(dbit,xchip)]);
                    if ~isinf(lags(i))
                        delta = lags(i) + moOffset{j,i} - chan_gt.hpre{j,i};
                    else
                        delta = lags_ce(i) + moOffset{j,i} - chan_gt.hpre{j,i};
                    end
                    yo{j} = RebuildKnownPacket(yo{j},dchip,chan_gt.hp{j,i},delta);
                    yo1{j} = RebuildKnownPacket(yo1{j},dchipTemp,chan.hp{j,i},delta);
                    yo2{j} = RebuildKnownPacket(yo2{j},dchipTemp,estimate.hp{j,i},delta);
                end
            end
    
            for j = 1:nMo
                subplot(nTx2+2,nMo,[nTx2*nMo+j,(nTx2+1)*nMo+j]); hold on; box on; grid on;
                plot(yr2{j});
                plot(yo{j});
                plot(yo1{j});
                plot(yo2{j});
                plot([endIdxCE,endIdxCE],ylim,'--k','LineWidth',2);
                xlabel('sample index');
                xticks(0:plen/2:swSize);
                ylabel('signal');
                legend('rx', 'gt', 'ce old', 'ce');
            end

            %%
            if exist("f1","var") && isvalid(f1), close(f1); end
            if isvalid(f2), close(f2); end
            break;
        end

        chan.hp = estimate.hp;
        chan.nb = estimate.nb;
    end
    
    %% decode all Tx
    lags_all = lags; lags_all(~isinf(lags_new)) = lags_new((~isinf(lags_new)));
    if min(lags_all) < plen
        if swEnd == ylen
            endIdxDC = swEnd-swStart+1;
        elseif ~isequal(algoCE,'gt') && sum(~isinf(lags_ce)) > 0
            endIdxDC = endIdxCE;
        else
            endIdxDC = min(plen + 2*swAdv+moOffsetMax,swEnd-swStart+1);
        end
        forwardIn = struct( ...
            'yd', {yr2}, 'xBit', {xBit}, ...
            'xChip', {xChip}, 'pChip', {pChip}, ...
            'lags', lags_all, ...
            'moOffset', {moOffset}, 'hPres', {chan.hpre}, ...
            'chan', chan, 'mode', 'de', ...
            'cpIdx', swAdv, 'nTracks', nTracks, ...
            'endIdx', endIdxDC, 'code', code);
        if exist('checkpoint', 'var')
            forwardIn.checkpoint = checkpoint;
        end
        forward = DecodeSequence_ViterbiForward4(forwardIn);
        checkpoint = forward.checkpoint;

        backwardIn = struct( ...
            'viterbi', forward.viterbi, ...
            'dBit', {dBit}, 'lags', lags_all);
        backward = DecodeSequence_ViterbiBackward(backwardIn);
        dBit = backward.dBit;
    else
        clear('checkpoint');
    end
    
    %% next window
    if swEnd == ylen
        break;
    end
    swStart = swStart + swAdv;
    swEnd = min(ylen, swStart + swSize-1);

    lags = lags_all - swAdv;
    lags(lags >= plen) = inf;
    lags_temp = lags_all - swAdv;
    lags_temp(lags_temp < plen) = inf;
end

%% return
PDOff = lags + swStart-1 - cell2mat(txOffset) - pkOffset;
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
        'labels', labels_pd, ...
        'errors', errors_pd, ...
        'corr', corr_pd, ...
        'ratio', ratio_pd, ...
        'ratio2', ratio2_pd);
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
wcorr = nan(1,nTx);
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
    
    % weighted corr
    sig1 = cell2mat(chan.hp(:,i));
    sig2 = cell2mat(esti.hp(:,i));
    wcorr(1,i) = dot(sig1,sig2) / (norm(sig1)*norm(sig2));
end

%% Tx power ratio
ratio = nan(nMo,nTx);
ratio2 = nan(nMo,nTx);
for i = 1:nTx
    if isinf(lags(i)), continue; end
    for j = 1:nMo
        ratio(j,i) = sum(esti.hp{j,i}) / sum(chan.hp{j,i});
        ratio2(j,i) = norm(esti.hp{j,i}) / norm(chan.hp{j,i});
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
rval.ratio2 = ratio2;
rval.wcorr = wcorr;
end