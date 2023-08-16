function rval = DecodeSequence_ViterbiBackward(params)
%% inputs
try
    viterbi = params.viterbi;
    tracedBits = viterbi.tracedBits;
    VStates = viterbi.VStates;
    VTraces = viterbi.VTraces;
    dBit = params.dBit;
    lags = params.lags;
catch
    error('decoding missing variables');
end

%%
nMo = length(tracedBits);
nTx = length(tracedBits{1});
nBit = length(dBit{1});

try
    stateInds = params.stateInds;
catch
    stateInds = cell(nMo,1);
    for j = 1:nMo
        [~, stateInds{j}] = max(cell2mat(VStates{j}(:,end)));
    end
end

%% backward
for j = 1:nMo
    backBits = tracedBits{j};
    backBits = min(nBit, backBits);

    stateInd = stateInds{j};
    for jy = size(VTraces,2)-1:-1:1
        trace = VTraces{j,jy+1}(stateInd,:);
        stateInd = trace(end);
        for i = 1:nTx
            if isinf(lags(i)), continue; end
            if ~isnan(trace(i))
                dBit{j,i}(backBits(i)) = trace(i);
                backBits(i) = backBits(i)-1;
            end
        end
    end
end

%%
rval.dBit = dBit;

end