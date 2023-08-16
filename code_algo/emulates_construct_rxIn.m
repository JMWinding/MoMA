function rxIn = emulates_construct_rxIn(params)
%%
isfigureCIR = false;
% isfigureCIR = true;
isfigureRx = false;
% isfigureRx = true;

%%
datanote = params.datanote;
T = params.T;
pumpstr = params.pumpstr;
Lp = params.Lp;
code = params.code;
notes = params.notes;

try type = params.type; catch, type = "ec"; end
try T2 = params.T2; catch, T2 = T; end % oversampling
assert(round(T/T2)*T2 == T, "chip interval must be divided by oversampling");
oversampling = T/T2;
try Lp2 = params.Lp2; catch, Lp2 = Lp; end % treat data bits as preamble
try hPre = params.hPre; catch, hPre = ceil(1e3/T2); end
try hPost = params.hPost; catch, hPost = ceil(3e3/T2); end
try hlen = params.hlen; catch, hlen = ceil(10e3/T2); end

%%
nTx = numel(strfind(pumpstr,"-"))+1;
nMo = size(notes,1);
try rxOffset = params.rxOffset; catch, rxOffset = zeros(nMo,1); end
if nMo > 1 && length(datanote) == 1
    datanote = repmat(datanote, [nMo,1]);
end
assert(nMo == length(datanote));

if type == "ph"
    Ka = 1.7e-5;
    try 
        pa = params.pa; pb = params.pb;
    catch 
        pa = -0.0268; pb = 15.5211;
    end
end

%% compute ground truth CIR
xChannel = struct( ...
    'xCIR', {cell(nMo,nTx)}, ...
    'noiseb', {cell(nMo,1)}, ...
    'noisen', {cell(nMo,1)}, ...
    'noisep', {cell(nMo,1)});

% CIR long/short range
for lenCIR = ["long", "short"]
switch lenCIR
    case "long"
        txOffset2 = zeros(1,nTx);
    case "short"
        hlen = hPre+hPost+1;
end

if isfigureCIR && lenCIR=="long"
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,nTx);
end
if isfigureRx && lenCIR=="long"
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,1);
end

for j = 1:nMo

foldername = datanote(j,1)+"/"+string(T)+"ms_"+pumpstr+"_"+string(Lp)+"_"+code;
[tx2, rx2] = read_txrx(foldername, notes(j,1), code);
yRx2 = rx2(2:T2:end,2);

if type == "ph"
    y0 = mean(yRx2(1:size(tx2.xChip,1)*oversampling));
    yRx2 = 10.^(-pa*yRx2-pb);
    yRx2 = (yRx2-10.^(-pa*y0-pb)) .* (1+yRx2/Ka);
end

H = ones(length(yRx2), 1+hlen*nTx);
for i = 1:nTx
    H(:,(1:hlen)+(i-1)*hlen+1) = get_decode_matrix( ...
        GenerateOversampleChips(tx2.chips(:,i),oversampling), ...
        length(yRx2), hlen, ...
        (tx2.moOffset(i)+tx2.txOffset(i))*oversampling+txOffset2(i));
end
h = H \ yRx2;

if isfigureCIR
    figure(f1);
    for i = 1:nTx
        nexttile(i+(lenCIR=="short")*nTx); box on; hold on; grid on;
        plot(h((1:hlen)+(i-1)*hlen+1));
    end
end
if isfigureRx
    figure(f2);
    nexttile(1+(lenCIR=="short")); box on; hold on; grid on;
    plot(yRx2); plot(H*h);
end

for i = 1:nTx
    xChannel.xCIR(j,i) = {h((1:hlen)+(i-1)*hlen+1)};
end

xChannel.noiseb(j,1) = {h(1)};
xChannel.noisen(j,1) = {0};
xChannel.noisep(j,1) = {mean(abs(yRx2-H*h)./sqrt(H*h))};
end

if isequal(lenCIR,"long")
    % adjust txOffset
    for i = 1:nTx
        temp = zeros(nMo,1);
        for j = 1:nMo
            [~, temp(j)] = max(xChannel.xCIR{j,i});
        end
        txOffset2(i) = round(median(temp))-1 - hPre;
    end
end

end

%% construct rxIn
rxIn = struct( ...
'yRx', {cell(nMo,1)}, ...
'xChannel', {xChannel}, ...
'xChip', {cell(nMo,nTx)}, ...
'xBit', {cell(nMo,nTx)}, ...
'pChip', {cell(nMo,nTx)}, ...
'txOffset', {cell(1,nTx)}, ...
'moOffset', {cell(nMo,nTx)}, ...
'Lp', Lp2, ...
'code', code, ...
'hPre', hPre, ...
'hPost', hPost);

for j = 1:nMo
foldername = datanote(j,1)+"/"+string(T)+"ms_"+pumpstr+"_"+string(Lp)+"_"+code;
[tx,rx] = read_txrx(foldername, notes(j,1), code);

if j == 1
    nChip = size(tx.xChip,1);
else
    assert(nChip == size(tx.xChip,1));
end
assert(Lp == tx.Lp);

yRx = rx(2:T2:end,2);
yRxp = zeros((nChip*rxIn.Lp+rxOffset(j))*oversampling,1);
for k = 1:length(yRxp)
    yRxp(k) = yRx(ceil(rand()*nChip*oversampling)+1);
end

if type == "ph"    
    y0 = mean(yRx(1:nChip*oversampling));
    yRxp = 10.^(-pa*yRxp-pb);
    yRx = 10.^(-pa*yRx-pb);
    yRxp = (yRxp-10.^(-pa*y0-pb)) .* (1+yRxp/Ka);
    yRx = (yRx-10.^(-pa*y0-pb)) .* (1+yRx/Ka);
end

rxIn.yRx(j,1) = {[yRxp; yRx]};
for i = 1:nTx
    rxIn.xChip(j,i) = {GenerateOversampleChips(2*tx.xChip(:,i)-1,oversampling)};
    switch code
        case "plain0"
            rxIn.xBit(j,i) = {2*tx.bits((Lp2-Lp)*nChip/2+1:end,i)-1};
        otherwise
            rxIn.xBit(j,i) = {2*tx.bits((Lp2-Lp)+1:end,i)-1};
    end
    rxIn.pChip(j,i) = {GenerateOversampleChips(2*tx.chips(1:Lp2*nChip,i)-1,oversampling)};
    rxIn.moOffset(j,i) = {(nChip*rxIn.Lp+tx.txOffset(1,i)+tx.moOffset(1,i)+rxOffset(j))*oversampling};
end
end
for i = 1:nTx
rxIn.txOffset(1,i) = {min(cell2mat(rxIn.moOffset(:,i)))};
for j = 1:nMo
    rxIn.moOffset(j,i) = {rxIn.moOffset{j,i} - rxIn.txOffset{1,i}};
end
rxIn.txOffset{1,i} = rxIn.txOffset{1,i} + txOffset2(i);
end

end