function rval = EstimateChannelMMoSW(params)

%% model parameters
try noisemodel = params.noisemodel; catch, noisemodel = 'pois'; end

%% inputs
try yr = params.yr; catch, error('channel estimation missing variables'); end
try xChip = params.xChip; catch, error('channel estimation missing variables'); end
try dBit = params.dBit; catch, error('channel estimation missing variables'); end
try pChip = params.pChip; catch, error('channel estimation missing variables'); end
try lags_ce = params.lags_ce; catch, error('channel estimation missing variables'); end
try algo = params.algo; catch, algo = 'af'; end
try code = params.code; catch, code = 'goldman'; end

%%
nMo = size(xChip,1);
nTx = size(xChip,2);
nChip = length(xChip{1});
plen = length(pChip{1});

%% other inputs
try sameMo = params.sameMo; catch, sameMo = false; end
try hp = params.hp; catch, hp = cell(nMo,nTx); end
try nb = params.nb; catch, nb = cell(nMo,1); nb(:) = {0}; end
try hPre = params.hPre; catch, hPre = nChip; end
try hPost = params.hPost; catch, hPost = 2*nChip; end
hTotal = hPre + hPost + 1;

try moOffset = params.moOffset; catch
    if nMo == 1
        moOffset = num2cell(zeros(1,nTx));
    else
        error('channel estimation missing fields');
    end
end

try lags_old = params.lags_old; catch, lags_old = inf(1,nTx); end
try endIdx = params.endIdx; catch, endIdx = 1.5*plen; end
% try sttIdx = params.sttIdx; catch, sttIdx = 1; end
try sttIdx = params.sttIdx; catch
    if isinf(lags_ce(1))
        sttIdx = 1;
    else
        sttIdx = max(1,lags_ce(1)+min(cell2mat(moOffset(:,1))));
    end
end

try afloss = params.weights.afloss; catch, afloss = "mean"; end
try nbest = params.nbest; catch, nbest = true; end
if isfield(params, "weights")
    try weight_pos = params.weights.pos; catch, weight_pos = 1; end
    try weight_posy = params.weights.posy; catch, weight_posy = 0; end
    try weight_simTx = params.weights.simTx; catch, weight_simTx = 10; end
    try weight_simMo = params.weights.simMo; catch, weight_simMo = 0; end
    try weight_smth = params.weights.smth; catch, weight_smth = 1; end
    try weight_cntr = params.weights.cntr; catch, weight_cntr = 1; end
end

%% rebuild signal from old Tx
yo = yr;
for j = 1:nMo
    yo{j}(:) = 0;
end
for i = 1:nTx
    if isinf(lags_old(i)), continue; end
    for j = 1:nMo
        dbit = dBit{j,i};
        xchip = GenerateCodeChips(xChip{j,i},code);
        dchip = ToPos([pChip{j,i}; GenerateDataChips(dbit,xchip)]);
        delta = lags_old(i) + moOffset{j,i};
        yo{j} = RebuildKnownPacket(yo{j},dchip,hp{j,i},delta);
    end
end

%% compute Toeplitz matrices
Xp = cell(nMo,nTx);
yp = cell(nMo,1);

lens_ce = zeros(nMo,1);
for j = 1:nMo
    lags_ce2 = lags_ce + cell2mat(moOffset(j,:));
    lens_ce(j) = endIdx - sttIdx + 1;
    
    for i = 1:nTx
        if isinf(lags_ce(i)), continue; end
        dbit = dBit{j,i};
        xchip = GenerateCodeChips(xChip{j,i},code);
        pchip = ToPos([pChip{j,i}; GenerateDataChips(dbit,xchip)]);
        Xp{j,i} = get_decode_matrix(pchip, lens_ce(j), hTotal, lags_ce2(i)-sttIdx+1);
    end
    yp{j} = yr{j}(sttIdx:endIdx);
    yo{j} = yo{j}(sttIdx:endIdx);
end

%% estimate channel
idx_ce = find(~isinf(lags_ce));
switch algo
    case {'ls'}
        if nbest
            for j = 1:nMo
                hpTemp = [cell2mat(Xp(j,idx_ce)), ones(lens_ce(j),1)] \ (yp{j} - yo{j});
                for ii = 1:length(idx_ce)
                    hp{j,idx_ce(ii)} = hpTemp((ii-1)*hTotal+(1:hTotal));
                end
                nb{j} = hpTemp(end);
            end
        else
            for j = 1:nMo
                hpTemp = [cell2mat(Xp(j,idx_ce))] \ (yp{j} - yo{j} - nb{j});
                for ii = 1:length(idx_ce)
                    hp{j,idx_ce(ii)} = hpTemp((ii-1)*hTotal+(1:hTotal));
                end
            end
        end
    case {'af'}
        % initialization
        for j = 1:nMo
            for i = 1:nTx
                if isinf(lags_ce(i)), continue; end
                if length(hp{j,i}) == hTotal, continue; end
                error('hp should have been initialized');
            end
        end
        firstRun = true;
        
        % iteration
        stepsize0 = 1e-2;

        stepsize = stepsize0;
        threshold = 1e-4;        
        hps2 = hp; bps2 = nb;
        
        afIn = struct('afloss', afloss, 'nbest', nbest, ...
            'yps', {yp}, 'yos', {yo}, 'Xps', {Xp}, ...
            'hps', {hp}, 'bps', {nb}, ...
            'weight_pos', weight_pos, 'weight_posy', weight_posy, ...
            'weight_simTx', weight_simTx, 'weight_simMo', weight_simMo, ...
            'weight_smth', weight_smth, 'weight_cntr', weight_cntr, ...
            'lags', lags_ce, 'noisemodel', noisemodel, 'sameMo', sameMo);
        rvalaf = get_stats_af(afIn);
        while 1            
            % update
            for j = 1:nMo
                for i = 1:nTx
                    if isinf(lags_ce(i)), continue; end
                    hps2{j,i} = hp{j,i} - stepsize * rvalaf.grad_h{j,i};
                end
                if nbest
                    bps2{j} = nb{j} - stepsize * rvalaf.grad_b{j};
                else
                    bps2{j} = nb{j};
                end
            end
            
            % which one is better?
            afIn2 = struct('afloss', afloss, 'nbest', nbest, ...
                'yps', {yp}, 'yos', {yo}, 'Xps', {Xp}, ...
                'hps', {hps2}, 'bps', {bps2}, ...
                'weight_pos', weight_pos, 'weight_posy', weight_posy, ...
                'weight_simTx', weight_simTx, 'weight_simMo', weight_simMo, ...
                'weight_smth', weight_smth, 'weight_cntr', weight_cntr, ...
                'lags', lags_ce, 'noisemodel', noisemodel, 'sameMo', sameMo);
            rvalaf2 = get_stats_af(afIn2);
            if rvalaf2.loss > rvalaf.loss
                stepsize = stepsize / 2;
                continue;
            end 
            
            % time to stop?
            if time_to_stop( ...
                    struct('hps', {hp}, 'bps', {nb}), ...
                    struct('hps', {hps2}, 'bps', {bps2}), ...
                    threshold)
                if firstRun
                    firstRun = false;
                    stepsize = stepsize0;

                    nb(:) = {0};
                    afIn = struct('afloss', afloss, 'nbest', nbest, ...
                        'yps', {yp}, 'yos', {yo}, 'Xps', {Xp}, ...
                        'hps', {hp}, 'bps', {nb}, ...
                        'weight_pos', weight_pos, 'weight_posy', weight_posy, ...
                        'weight_simTx', weight_simTx, 'weight_simMo', weight_simMo, ...
                        'weight_smth', weight_smth, 'weight_cntr', weight_cntr, ...
                        'lags', lags_ce, 'noisemodel', noisemodel, 'sameMo', sameMo);
                    rvalaf = get_stats_af(afIn);
                    continue;
                else
                    break;
                end
            end
            
            % next iteration
            hp = hps2;
            nb = bps2;
            rvalaf = rvalaf2;
        end
    otherwise
        error('channel estimation algorithm not supported');
end

%% post process
nn = cell(nMo,1);
np = cell(nMo,1);
for j = 1:nMo
    % non-negative
    for i = 1:nTx
        if isinf(lags_ce(i)), continue; end
%         hp{j,i} = max(0, hp{j,i});
    end
%     nb{j} = max(0, nb{j});
    % noise
    if isempty(idx_ce)
        yh = yo{j} + nb{j};
    else
        yh = yo{j} + nb{j} ...
            + cell2mat(Xp(j,idx_ce)) * cell2mat(hp(j,idx_ce).');
    end
    idx = yh > max(yh) * 0.1;
    switch noisemodel
        case 'norm'
            nn{j} = std(yh-yp{j});
            np{j} = 0;
        case {'pois', 'pois0'}
            nn{j} = 0;
            np{j} = std((yh(idx)-yp{j}(idx))./sqrt(yh(idx)));
        case 'mix'
            matA = [1, mean(yh(idx)); mean(1./yh(idx)), 1];
            vecb = [sum((yh(idx)-yp{j}(idx)).^2); ...
                sum((yh(idx)-yp{j}(idx)).^2./yh(idx))];
            sol = matA \ vecb;
            assert(min(sol)>=0);
            nn{j} = sol(1).^0.5;
            np{j} = sol(2).^0.5;
    end
end

%% return
rval.hp = hp;
rval.nb = nb;
rval.nn = nn;
rval.np = np;

rval.hPre = hPre;
rval.hPost = hPost;

rval.stepsize = stepsize;
rval.sttIdx = sttIdx;
rval.endIdx = endIdx;

end

%% support functions
function ret = relu(input)
ret = input .* (input > 0);
end

function rval = time_to_stop(sol, sol_new, thld)
nMo = size(sol.hps,1);

for j = 1:nMo
    temp = max(abs([cell2mat(sol.hps(j,:).')-cell2mat(sol_new.hps(j,:).'); ...
        sol.bps{j}-sol_new.bps{j}]));
    if temp > thld
        rval = false;
        return;
    end
end
rval = true;

end

function rval = get_stats_af(params)
    switch params.afloss
        case "mean"
            rval = get_stats_af_mean(params);
        case "sum"
            rval = get_stats_af_sum(params);
    end
end

function rval = get_stats_af_mean(params)
%% inputs
yps = params.yps;
yos = params.yos;
Xps = params.Xps;
hps = params.hps;
bps = params.bps;
weight_pos = params.weight_pos;
weight_posy = params.weight_posy;
weight_simTx = params.weight_simTx;
weight_simMo = params.weight_simMo;
weight_smth = params.weight_smth;
weight_cntr = params.weight_cntr;
lags = params.lags;
noisemodel = params.noisemodel;
try sameMo = params.sameMo; catch, sameMo = false; end

%% related variables
nMo = size(Xps,1);
nTx = size(Xps,2);
hTotal = length(hps{1,find(~isinf(lags),1)});

loss = 0;
grad_h = cell(nMo,nTx); grad_h(:) = {0};
grad_b = cell(nMo,1); grad_b(:) = {0};

%% useful values
yhs = cell(nMo,1);
idxs = cell(nMo,1);
ypsErr = cell(nMo,1);
hpsTx = cell(1,nTx); hpsTx(:) = {0};
hpsMo = cell(nMo,1); hpsMo(:) = {0};
for j = 1:nMo
    yhs{j} = yos{j} + bps{j} + ...
        cell2mat(Xps(j,~isinf(lags))) * cell2mat(hps(j,~isinf(lags)).');
    idxs{j} = yhs{j} > max(yhs{j}) * 0.1;
    ypsErr{j} = yhs{j} - yps{j};
    for i = 1:nTx
        if length(hps{j,i}) ~= hTotal, continue; end
        if nMo > 1
            if sameMo
                hpsTx{1,i} = hpsTx{1,i} + hps{j,i} / nMo;
            else
                hpsTx{1,i} = hpsTx{1,i} + (hps{j,i} / norm(hps{j,i})) / nMo;
            end
        end
        if nTx > 1
            if sameMo
                hpsMo{j,1} = hpsMo{j,1} + hps{j,i} / nTx;
            else
                hpsMo{j,1} = hpsMo{j,1} + (hps{j,i} / norm(hps{j,i})) / nTx;
            end
        end
    end
end

%% loss
for j = 1:nMo
    % LS term
    switch noisemodel
        case {'norm', 'pois'}
            loss = loss + mean(ypsErr{j}.^2) ...
                + weight_posy * mean(relu(ypsErr{j}).^2);
        case 'pois0'
            loss = loss + mean(ypsErr{j}(idxs{j}).^2 ./ yhs{j}(idxs{j})) ...
                + weight_posy * mean(relu(ypsErr{j}(idxs{j})).^2 ./ yhs{j}(idxs{j}));
    end
    % non-negative term
    loss = loss + weight_pos * sum(relu([-cell2mat(hps(j,:).');-bps{j}]).^2) / hTotal;
    % CIR similarity term
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        if nMo > 1
            if sameMo
                loss = loss + weight_simTx * mean((hps{j,i}-hpsTx{1,i}).^2);
            else
                loss = loss + weight_simTx * mean((hps{j,i}-hpsTx{1,i}*norm(hps{j,i})).^2);
            end
        end
        if nTx > 1
            if sameMo
                loss = loss + weight_simMo * mean((hps{j,i}-hpsMo{j,1}).^2);
            else
                loss = loss + weight_simMo * mean((hps{j,i}-hpsMo{j,1}*norm(hps{j,i})).^2);
            end
        end
    end
    % other terms
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        temp = diff([0; hps{j,i}; 0]);
        [~, tempi] = max(hps{j,i});
        % smoothness term
        loss = loss + weight_smth * sum(relu(-temp(1:end-1).*temp(2:end))) / hTotal;
        % center term
        loss = loss + weight_cntr ...
            * sum((hps{j,i} .* ((1:hTotal).'-tempi)).^2) / hTotal.^2;
    end
end

%% gradient
for j = 1:nMo
    %% background
    % LS term
    switch noisemodel
        case {'norm', 'pois'}
            grad_b{j} = grad_b{j} + (mean(ypsErr{j})) ...
                + weight_posy * (mean(relu(ypsErr{j})));
        case 'pois0'
            e2yh = ypsErr{j}(idxs{j}) ./ yhs{j}(idxs{j});
            grad_b{j} = grad_b{j} + (mean(e2yh.*(e2yh+2))) ...
                + weight_posy * (mean(relu(e2yh).*(e2yh+2)));
    end
    % non-negative term
    grad_b{j} = grad_b{j} + (-weight_pos * relu(-bps{j}) / hTotal);

    %% CIR
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        % LS term
        switch noisemodel
            case {'norm', 'pois'}
                grad_h{j,i} = grad_h{j,i} + (Xps{j,i}' * ypsErr{j} / hTotal) ...
                    + weight_posy * (Xps{j,i}' * relu(ypsErr{j}) / hTotal);
            case 'pois0'
                e2yh = ypsErr{j}(idxs{j}) ./ yhs{j}(idxs{j});
                grad_h{j,i} = grad_h{j,i} + ...
                    (Xps{j,i}(idxs{j},:)' * (e2yh.*(e2yh+2)) / hTotal) ...
                    + weight_posy * (Xps{j,i}(idxs{j},:)' * (relu(e2yh).*(e2yh+2)) / hTotal);
        end
        % non-negative term
        grad_h{j,i} = grad_h{j,i} + (-weight_pos * relu(-hps{j,i}) / hTotal);
        % CIR similarity term
        if nMo > 1
            if sameMo
                grad_h{j,i} = grad_h{j,i} + (weight_simTx * (hps{j,i}-hpsTx{1,i}) / hTotal);
            else
                grad_h{j,i} = grad_h{j,i} + weight_simTx ...
                    * ((1-1/nMo)*eye(hTotal) + hps{j,i}/norm(hps{j,i})*(hps{j,i}/norm(hps{j,i})/nMo-hpsTx{1,i}).') ...
                    * (hps{j,i} - norm(hps{j,i})*hpsTx{1,i}) / hTotal;
            end
        end
        if nTx > 1
            if sameMo
                grad_h{j,i} = grad_h{j,i} + (weight_simMo * (hps{j,i}-hpsMo{j,1}) / hTotal);
            else
                grad_h{j,i} = grad_h{j,i} + weight_simMo ...
                    * ((1-1/nTx)*eye(hTotal) + hps{j,i}/norm(hps{j,i})*(hps{j,i}/norm(hps{j,i})/nTx-hpsMo{j,1}).') ...
                    * (hps{j,i} - norm(hps{j,i})*hpsMo{j,1}) / hTotal;            
            end
        end
        % smoothness term
        temp = diff([0; hps{j,i}; 0]);
        temp2 = -temp(1:end-1).*temp(2:end);
        [~, tempi] = max(hps{j,i});
        grad_h{j,i} = grad_h{j,i} + weight_smth ...
            * ((temp2>0) .* (temp(1:end-1)-temp(2:end))...
              + ([0; temp2(1:end-1)]>0) .* [0; -temp(1:end-2)] ...
              + ([temp2(2:end); 0]>0)   .* [temp(3:end); 0]) / hTotal;
        % center term
        grad_h{j,i} = grad_h{j,i} + weight_cntr ...
            * hps{j,i} .* ((1:hTotal).'-tempi).^2 / hTotal.^2;
    end
end

%% return
rval.loss = loss;
rval.grad_h = grad_h;
rval.grad_b = grad_b;

end

function rval = get_stats_af_sum(params)
%% inputs
yps = params.yps;
yos = params.yos;
Xps = params.Xps;
hps = params.hps;
bps = params.bps;
weight_pos = params.weight_pos;
weight_posy = params.weight_posy;
weight_simTx = params.weight_simTx;
weight_simMo = params.weight_simMo;
weight_smth = params.weight_smth;
weight_cntr = params.weight_cntr;
lags = params.lags;
noisemodel = params.noisemodel;
try sameMo = params.sameMo; catch, sameMo = false; end

%% related variables
nMo = size(Xps,1);
nTx = size(Xps,2);
hTotal = length(hps{1,find(~isinf(lags),1)});

loss = 0;
grad_h = cell(nMo,nTx); grad_h(:) = {0};
grad_b = cell(nMo,1); grad_b(:) = {0};

%% useful values
yhs = cell(nMo,1);
idxs = cell(nMo,1);
ypsErr = cell(nMo,1);
hpsTx = cell(1,nTx); hpsTx(:) = {0};
hpsMo = cell(nMo,1); hpsMo(:) = {0};
for j = 1:nMo
    yhs{j} = yos{j} + bps{j} + ...
        cell2mat(Xps(j,~isinf(lags))) * cell2mat(hps(j,~isinf(lags)).');
    idxs{j} = yhs{j} > max(yhs{j}) * 0.1;
    ypsErr{j} = yhs{j} - yps{j};
    for i = 1:nTx
        if length(hps{j,i}) ~= hTotal, continue; end
        if nMo > 1
            if sameMo
                hpsTx{1,i} = hpsTx{1,i} + hps{j,i} / nMo;
            else
                hpsTx{1,i} = hpsTx{1,i} + (hps{j,i} / norm(hps{j,i})) / nMo;
            end
        end
        if nTx > 1
            if sameMo
                hpsMo{j,1} = hpsMo{j,1} + hps{j,i} / nTx;
            else
                hpsMo{j,1} = hpsMo{j,1} + (hps{j,i} / norm(hps{j,i})) / nTx;
            end
        end
    end
end

%% loss
for j = 1:nMo
    % LS term
    switch noisemodel
        case {'norm', 'pois'}
            loss = loss + sum(ypsErr{j}.^2) ...
                + weight_posy * sum(relu(ypsErr{j}).^2);
        case 'pois0'
            loss = loss + sum(ypsErr{j}(idxs{j}).^2 ./ yhs{j}(idxs{j})) ...
                + weight_posy * sum(relu(ypsErr{j}(idxs{j})).^2 ./ yhs{j}(idxs{j}));
    end
    % non-negative term
    loss = loss + weight_pos * sum(relu([-cell2mat(hps(j,:).');-bps{j}]).^2);
    % CIR similarity term
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        if nMo > 1
            if sameMo
                loss = loss + weight_simTx * sum((hps{j,i}-hpsTx{1,i}).^2);
            else
                loss = loss + weight_simTx * sum((hps{j,i}-hpsTx{1,i}*norm(hps{j,i})).^2);
            end
        end
        if nTx > 1
            if sameMo
                loss = loss + weight_simMo * sum((hps{j,i}-hpsMo{j,1}).^2);
            else
                loss = loss + weight_simMo * sum((hps{j,i}-hpsMo{j,1}*norm(hps{j,i})).^2);
            end
        end
    end
    % other terms
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        temp = diff([0; hps{j,i}; 0]);
        [~, tempi] = max(hps{j,i});
        % smoothness term
        loss = loss + weight_smth * sum(relu(-temp(1:end-1).*temp(2:end)));
        % center term
        loss = loss + weight_cntr ...
            * sum((hps{j,i} .* ((1:hTotal).'-tempi)).^2) / hTotal;
    end
end

%% gradient
for j = 1:nMo
    %% background
    % LS term
    switch noisemodel
        case {'norm', 'pois'}
            grad_b{j} = grad_b{j} + (sum(ypsErr{j})) ...
                + weight_posy * (sum(relu(ypsErr{j})));
        case 'pois0'
            e2yh = ypsErr{j}(idxs{j}) ./ yhs{j}(idxs{j});
            grad_b{j} = grad_b{j} + (sum(e2yh.*(e2yh+2))) ...
                + weight_posy * (sum(relu(e2yh).*(e2yh+2)));
    end
    % non-negative term
    grad_b{j} = grad_b{j} + (-weight_pos * relu(-bps{j}));

    %% CIR
    for i = 1:nTx
        if isinf(lags(i)), continue; end
        % LS term
        switch noisemodel
            case {'norm', 'pois'}
                grad_h{j,i} = grad_h{j,i} + (Xps{j,i}' * ypsErr{j}) ...
                    + weight_posy * (Xps{j,i}' * relu(ypsErr{j}));
            case 'pois0'
                e2yh = ypsErr{j}(idxs{j}) ./ yhs{j}(idxs{j});
                grad_h{j,i} = grad_h{j,i} + ...
                    (Xps{j,i}(idxs{j},:)' * (e2yh.*(e2yh+2))) ...
                    + weight_posy * (Xps{j,i}(idxs{j},:)' * (relu(e2yh).*(e2yh+2)));
        end
        % non-negative term
        grad_h{j,i} = grad_h{j,i} + (-weight_pos * relu(-hps{j,i}));
        % CIR similarity term
        if nMo > 1
            if sameMo
                grad_h{j,i} = grad_h{j,i} + (weight_simTx * (hps{j,i}-hpsTx{1,i}));
            else
                grad_h{j,i} = grad_h{j,i} + weight_simTx ...
                    * ((1-1/nMo)*eye(hTotal) + hps{j,i}/norm(hps{j,i})*(hps{j,i}/norm(hps{j,i})/nMo-hpsTx{1,i}).') ...
                    * (hps{j,i} - norm(hps{j,i})*hpsTx{1,i});
            end
        end
        if nTx > 1
            if sameMo
                grad_h{j,i} = grad_h{j,i} + (weight_simMo * (hps{j,i}-hpsMo{j,1}));
            else
                grad_h{j,i} = grad_h{j,i} + weight_simMo ...
                    * ((1-1/nTx)*eye(hTotal) + hps{j,i}/norm(hps{j,i})*(hps{j,i}/norm(hps{j,i})/nTx-hpsMo{j,1}).') ...
                    * (hps{j,i} - norm(hps{j,i})*hpsMo{j,1});            
            end
        end
        % smoothness term
        temp = diff([0; hps{j,i}; 0]);
        temp2 = -temp(1:end-1).*temp(2:end);
        [~, tempi] = max(hps{j,i});
        grad_h{j,i} = grad_h{j,i} + weight_smth ...
            * ((temp2>0) .* (temp(1:end-1)-temp(2:end))...
              + ([0; temp2(1:end-1)]>0) .* [0; -temp(1:end-2)] ...
              + ([temp2(2:end); 0]>0)   .* [temp(3:end); 0]);
        % center term
        grad_h{j,i} = grad_h{j,i} + weight_cntr ...
            * hps{j,i} .* ((1:hTotal).'-tempi).^2 / hTotal;
    end
end

%% return
rval.loss = loss;
rval.grad_h = grad_h;
rval.grad_b = grad_b;

end
