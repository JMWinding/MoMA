function rval = PreambleDetectionSW8(params)
% Double correlation

%% inputs
try 
    yp = params.yp;
    xChip = params.xChip;
    pChip = params.pChip;
    Lp = params.Lp;
    dBit = params.dBit;
catch
    error("preamble detection missing inputs");
end

try isrepeat = params.isrepeat; catch, isrepeat = false; end
try algo = params.algo; catch, algo = "sc"; end
try algoCE = params.algoCE; catch, algoCE = "gt"; end
try code = params.code; catch, code = "goldman"; end

%% related variables
nMo = size(xChip,1);
nTx = size(xChip,2);
nChip = length(xChip{1});
assert(size(yp,1) == nMo);

try moOffset = params.moOffset; catch
    if nMo == 1
        moOffset = num2cell(zeros(1,nTx));
    else
        error("missing fields");
    end
end
assert(sum(sum(cell2mat(moOffset)==0)<1)==0);

try
    lags = params.lags; % lags relative to current windowed signal
    hPres = params.hPres;
    hPre = params.hPre;
    hPost = params.hPost;
    hp = params.hp;
    
    assert(size(hp,1) == nMo);
    assert(size(hp,2) == nTx);
catch
    lags = inf(1,nTx);
    hPres = cell(nMo,nTx);
    hPre = 0;
    hPost = 0;
    hp = cell(nMo,nTx);
end

try 
    lags_min = params.lags_min;
    lags_max = params.lags_max;
catch
    lags_min = zeros(1,nTx);
    lags_max = inf(1,nTx);
end
try lags_gt = params.lags_gt; catch, lags_gt = inf(1,nTx); end

PDVal = inf(1,nTx);

%% packet detection
plen = nChip*Lp;
mlen = 2.5*plen; % if swAdv changes, this should also be changed

if mlen > plen && sum(isinf(lags)) > 0
    %% remove signal from known packets
    yp2 = yp;
    for i = 1:nTx
        if ~isinf(lags(i))
            for j = 1:nMo
                if isrepeat
                    dbit = dBit{1,i};
                else
                    dbit = dBit{j,i};
                end
                xchip = GenerateCodeChips(xChip{j,i},code);
                dchip = ToPos([pChip{j,i}; GenerateDataChips(dbit,xchip)]);
                delta = lags(i) + moOffset{j,i} - hPres{j,i};
                yp2{j} = RemoveKnownPacket(yp2{j},dchip,hp{j,i},delta);
            end
        end
    end

    %% compute cross-correlation
    pcorr = nan(mlen,nMo,nTx);
    pcorr2 = nan(mlen,nMo,nTx);
    for i = 1:nTx
        if ~isinf(lags(i)), continue; end
        
        for j = 1:nMo
            pchip = pChip{j,i};
            % corr
            if isequal(pChip{1}(1:Lp), xChip{1})
                siga = GeneratePreambleDetectionChips(pchip,nChip);
            else
%                 siga = GeneratePreambleDetectionChips(pchip,Lp/2);
                siga = GeneratePreambleDetectionChips(pchip,1);
                siga(plen+1:end) = [];
            end
            temp = conv(yp2{j,1}(moOffset{j,i}+1:end),flip(siga)) / norm(siga);
            pcorr(:,j,i) = temp(length(siga)+(0:size(pcorr,1)-1));
            % corr of corr (corr2)
            sigb = conv(siga,flip(siga)); sigb(end-length(siga)+2:end) = [];
            for k = 1:size(pcorr,1)
                sig2 = pcorr(max(0,k-length(sigb))+1:k,j,i);
                pcorr2(k,j,i) = dot(sig2,sigb(end-length(sig2)+1:end)) ...
                    / (norm(sig2)*norm(sigb(end-length(sig2)+1:end)));
            end
        end
    end

    %% compute metric
    pcorrTx = squeeze(sum(pcorr,2));
    pcorr2Tx = squeeze(sum(pcorr2,2));
    metric = pcorrTx .* pcorr2Tx .* (1-2*(pcorrTx<0 & pcorr2Tx<0));
    metric_final = metric;
    
    %% search for potential packets
    lags_detected = inf(1,nTx);
    switch algo
        case "gt"
            for i = 1:nTx
                if isinf(lags_gt(i)), continue; end

                if lags_gt(i) >= plen && lags_gt(i) < 2*plen
                    lags_detected(i) = lags_gt(i);
                end
            end

        case "sc"
            pksTemp = -inf(1,nTx);
            lagsTemp = inf(1,nTx);
            for i = 1:nTx
                if ~isinf(lags(i)) || lags_min(i)>=2*plen || lags_max(i)<plen
                    continue;
                end

                lagMin = max(lags_min(i),plen);
                lagMax = min(lags_max(i),2*plen-1);
                [pkTemp, lagTemp] = max(metric_final(lagMin+1:lagMax+1,i),[],"omitnan");
                lagTemp = lagTemp + lagMin;

                if pkTemp < max(metric_final(lagTemp+1:min(lagTemp+2*Lp*2,size(metric_final,1)),i))
                    continue;
                end

                pksTemp(i) = pkTemp;
                lagsTemp(i) = lagTemp;
            end

            [~, maxInd] = sort(pksTemp,"descend");
            for i = 1:nTx
                if isinf(pksTemp(maxInd(i))), continue; end
                ind = find(lagsTemp>lagsTemp(maxInd(i))-plen ...
                    & lagsTemp<lagsTemp(maxInd(i))+plen);
                lagsTemp(ind(pksTemp(ind)<pksTemp(maxInd(i))*0.6)) = inf;
            end
            ind = find(lagsTemp<=2*plen);
            lags_detected(ind) = lagsTemp(ind);
            PDVal(ind) = pksTemp(ind);
            lags_detected = lags_detected - 1;
    end

    %% debug
    while (~isequal(lags_detected, lags_gt) ... 
            || sum(~isinf(lags_gt)) > 0 ...
          ) && false
        if isempty(groot().Children)
            figure("units","normalized","outerposition",[0 0 1 1]);
        end
        tiledlayout(3,1,"TileSpacing","compact");
        nexttile;
        plot(pcorrTx); title("corr"); box on; grid on;
        xticks(0:plen/2:mlen); xlim([0,mlen]);
        nexttile;
        plot(pcorr2Tx); title("corr2"); box on; grid on;
        xticks(0:plen/2:mlen); xlim([0,mlen]);
        nexttile;
        plot(metric); title("metric"); box on; grid on;
        xticks(0:plen/2:mlen); xlim([0,mlen]);
        lg = legend(cellstr(num2str([lags; lags_detected; lags_gt].')));
        lg.Layout.Tile = "east";

        %%
%         lags_detected = lags_new;
        break;
    end
    
    %% compute related values
    lags(~isinf(lags_detected)) = lags_detected(~isinf(lags_detected));
    if sum(~isinf(lags_detected)) > 0 && ~isequal(algoCE, "gt")
        for j = 1:nMo
            temp = mean(cell2mat(hp(j,:)),2);
            for i = 1:nTx
                if isinf(lags_detected(i)), continue; end
                if ~isempty(hp{j,i}), continue; end

                hPres{j,i} = hPre;
                if isempty(temp)
                    hp{j,i} = [zeros(hPre,1); pksTemp(i)/nMo; zeros(hPost,1)];
                else
                    hp{j,i} = temp;
                end
            end
        end
    end
        
end

%% return
rval.lags = lags;
rval.hPres = hPres;
rval.hp = hp;
rval.PDVal = PDVal;

end