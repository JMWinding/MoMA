function rval = DecodeSequence_ViterbiForward4(params)
%% inputs
try 
    yd = params.yd;
    xBit = params.xBit;
    xChip = params.xChip;
    pChip = params.pChip;
    lags = params.lags;
    moOffset = params.moOffset;
    hPres = params.hPres;
    chan = params.chan;
catch
    error("decoding missing variables");
end

%%
nMo = size(xChip,1);
nTx = size(xChip,2);
nBit = length(xBit{1});

% lags_new means new Tx and not ground truth PD
try lags_new = params.lags_new; catch, lags_new = inf(1,nTx); end
try cpIdx = params.cpIdx; catch, cpIdx = 0; end
try endIdx = params.endIdx; catch, endIdx = length(yd{1}); end
try mode = params.mode; catch, mode = "dc"; end

try code = params.code; catch, code = "goldman"; end
for j = 1:nMo
    for i = 1:nTx
        xChip{j,i} = GenerateCodeChips(xChip{j,i},code);
    end
end
nChip = length(xChip{1});

%% inherit states from last window
try nTracks = params.nTracks; catch, nTracks = 2^8; end
try 
    checkpoint = params.checkpoint;
    tracedBits = checkpoint.tracedBits;
    VTraceLast = checkpoint.VTraceLast;
    VStateLast = checkpoint.VStateLast;
catch
    % the 2nd index is for signal index
    tracedBits = cell(nMo,1);
    tracedBits(:) = {-inf(1,nTx)};
    VTraceLast = cell(nMo,1);
    VTraceLast(:) = {[nan(1,nTx), 0, -1]};
    VStateLast = cell(nMo,1);
    VStateLast(:) = {[cell(1,nTx), {0}]};
end

%% compute channel
hpre = cell(nMo,nTx);
hpcs = cell(nMo,nTx,nChip); % CIR for code "c" (bit 1)
hprs = cell(nMo,nTx,nChip); % CIR for code "1-c" (bit -1)
for j = 1:nMo
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        hpre(j,i) = {conv(chan.hp{j,i}, ToPos(pChip{j,i}))};
        xchip = xChip{j,i};
        hpcFull = conv(chan.hp{j,i}, ToPos(xchip));
        hprFull = conv(chan.hp{j,i}, ToPos(-xchip));
        for k = 1:nChip
            hpcs{j,i,k} = hpcFull(mod(k-1,nChip)+1:nChip:end);
            hprs{j,i,k} = hprFull(mod(k-1,nChip)+1:nChip:end);
        end
    end
end

%% initialize states
VTraces = cell(nMo,endIdx+1);
VTraces(:,1) = VTraceLast;        
VStates = VStateLast;
for j = 1:nMo
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        if ~isinf(tracedBits{j}(i)), continue; end

        tracedBits{j}(i) = 0;

        hLen = 0;
        for k = 1:nChip
            hLen = max(hLen, length(hpcs{j,i,k}));
            hLen = max(hLen, length(hprs{j,i,k}));
        end
        for k = 1:size(VStates{j},1)
            VStates{j}{k,i} = zeros(1,hLen);
        end
    end
end

%% forward

% forward
for jy = 1:endIdx
    for j = 1:nMo
        %% init
        VStatesTemp = VStates{j};
        VTraces{j,jy+1} = nan(size(VStates{j},1), 2+nTx);
        VTraces{j,jy+1}(:,end) = (1:size(VTraces{j,jy+1},1)).';

        %% add new bit
        for i = 1:nTx
            if isinf(lags(i)), continue; end
            % mark the start of the preamble
            if lags(i)+moOffset{j,i}-hPres{j,i} + 1 == jy
                if isinf(lags_new(i))
                    % gt PD, or confirmed packet
                    for k = 1:size(VStatesTemp,1)
                        % this is a flag when the packet has not gone to
                        % its data part
                        VStatesTemp{k,i}(end) = 1;
                    end
                else
                    % unconfirmed packet
                    VStatesTemp = repmat(VStatesTemp, [2,1]);
                    VTraces{j,jy+1} = repmat(VTraces{j,jy+1}, [2,1]);
                    for k = size(VStatesTemp,1)/2+1:size(VStatesTemp,1)
                        VStatesTemp{k,i}(end) = 1;
                    end
                end
                continue;
            end

            % add a new bit for the data part
            if lags(i)+moOffset{j,i}-hPres{j,i} + length(pChip{j,i})+1 > jy
                % data part has not started
                continue;
            end
            if mod(jy-(lags(i)+moOffset{j,i}-hPres{j,i} + length(pChip{j,i})+1),nChip) ~= 0
                % new bit is not transmitted (still in a chip)
                continue;
            end

            tracedBits{j}(i) = tracedBits{j}(i) + 1;
            assert(tracedBits{j}(i) > 0);
            % split the state
            if tracedBits{j}(i) <= nBit
                VStatesTemp = repmat(VStatesTemp, [2,1]);
                VTraces{j,jy+1} = repmat(VTraces{j,jy+1}, [2,1]);
            end
            % add the new bit to each state
            for k = 1:size(VStatesTemp,1)
                if tracedBits{j}(i) <= nBit
                    if VStatesTemp{k,i}(1) == 0 && VStatesTemp{k,i}(end) == 0
                        % assume that the packet does not exist in this state
                        newBit = 0;
                    elseif k > size(VStatesTemp,1)/2
                        newBit = 1;
                    else
                        newBit = -1;
                    end
                    VTraces{j,jy+1}(k,i) = newBit;
                else
                    newBit = 0;
                end
                VStatesTemp{k,i}(2:end) = VStatesTemp{k,i}(1:end-1);
                VStatesTemp{k,i}(1) = newBit;
            end
        end

        %% remove excessive states (and update probability)
        % set repeated states with lower probability to inf
        % (to be deleted)
        VStatesTemp2 = cell2mat(VStatesTemp(:,1:end-1));
        [~, ind] = sortrows(VStatesTemp2);
        for k = length(ind):-1:2
            if isequal(VStatesTemp2(ind(k),:), VStatesTemp2(ind(k-1),:))
                if VStatesTemp{ind(k),end} > VStatesTemp{ind(k-1),end}
                    VStatesTemp{ind(k-1),end} = -inf;
                    ind(k-1) = ind(k);
                else
                    VStatesTemp{ind(k),end} = -inf;
                end
            end
        end
        clear("VStatesTemp2");

        % delete inf
        ind = isinf(cell2mat(VStatesTemp(:,end)));
        VStatesTemp(ind,:) = [];
        VTraces{j,jy+1}(ind,:) = [];
        if size(VStatesTemp,1) == 0
            error("viterbi reaches all states with zero probability");
        end

        % probability of current state
        ViterbiLogProbCurrent = nan(size(VStatesTemp,1),1);
        for k = 1:size(VStatesTemp,1)
            yh = chan.nb{j};

            for i = 1:nTx
                if isinf(lags(i)), continue; end
                if VStatesTemp{k,i}(1) == 0 && VStatesTemp{k,i}(end) == 0
                    % assume that the packet does not exist in this state
                    continue;
                end
                
                % preamble                
                pind = jy - (lags(i)+moOffset{j,i}-hPres{j,i});
                if pind > 0 && pind <= length(hpre{j,i})
                    yh = yh + hpre{j,i}(pind);
                end

                % data bits
                if tracedBits{j}(i) > 0
                    % bit 1 correpsonds to "hpc"
                    hpc = hpcs{j,i,mod(jy-(lags(i)+moOffset{j,i}-hPres{j,i})-1,nChip)+1};
                    yh = yh + dot(hpc, double(VStatesTemp{k,i}(1:length(hpc))>0));
                    % bit -1 in preamble corresponds to "hpr"
                    hpr = hprs{j,i,mod(jy-(lags(i)+moOffset{j,i}-hPres{j,i})-1,nChip)+1};
                    yh = yh + dot(hpr, double(VStatesTemp{k,i}(1:length(hpr))<0));
                    % bit 0 in data correspond to nothing
                end
            end

            switch mode
                case "pd"
                    ViterbiLogProbCurrent(k) = -abs(yd{j}(jy) - yh);
                case "ce"
                    ViterbiLogProbCurrent(k) =  ComputeViterbiLogProb(yd{j}(jy), yh, ...
                          10*chan.nn{j}, 10*chan.np{j});
                case "dc"
                    ViterbiLogProbCurrent(k) = ComputeViterbiLogProb(yd{j}(jy), yh, ...
                          chan.nn{j}, chan.np{j});
            end
        end

        % handle -inf
        ind = isinf(ViterbiLogProbCurrent);
        if sum(~ind) == 0
            % if all probabilities are too small, skip this sample
        else
            % update joint probability
            for k = 1:size(VStatesTemp,1)
                VStatesTemp{k,end} = VStatesTemp{k,end} + ViterbiLogProbCurrent(k);
            end
            VStatesTemp(ind,:) = [];
            VTraces{j,jy+1}(ind,:) = [];
        end

        % sort and truncate lower probability
        [~, ind] = sort(cell2mat(VStatesTemp(:,end)), "descend");
        if length(ind) > nTracks
            ind(nTracks+1:end) = [];
        end
        VStatesTemp = VStatesTemp(ind,:);
        VStates{j} = VStatesTemp;
        VTraces{j,jy+1} = VTraces{j,jy+1}(ind,:);
        VTraces{j,jy+1}(:,end-1) = cell2mat(VStates{j}(:,end));
    end

    %% values for next window (only once)
    if jy == cpIdx
        checkpoint = struct( ...
            "tracedBits", {tracedBits}, ...
            "VTraceLast", {VTraces(:,jy+1)}, ...
            "VStateLast", {VStates});
    end
end

%%
if exist("checkpoint", "var")
    rval.checkpoint = checkpoint;
end
rval.viterbi = struct( ...
    "tracedBits", {tracedBits}, ...
    "VStates", {VStates}, ...
    "VTraces", {VTraces});

end