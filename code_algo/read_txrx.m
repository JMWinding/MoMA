function [tx, rx] = readtxrx(folder, note, code)
%% read tx
if contains(code, "goldman") || contains(code, "plain")
    settingFile = folder+"/../goldman.txt";
elseif contains(code, "gold")
    settingFile = folder+"/../gold.txt";
elseif contains(code, "ooc")
    settingFile = folder+"/../ooc.txt";
else
    error("code not defined");
end
txFile = folder+"/tx/"+note+".TXT";
tx = readtx(txFile, settingFile);

%% read rx
rxFile = folder+"/rx/"+note+".TXT";
rx = readrx(rxFile);

end

function tx = readtx(txFile, settingFile)
if contains(settingFile, "goldman")
    code = "goldman";
else
    % applies to ooc as well
    code = "gold";
end

%% CDMA code
codeId = fopen(settingFile, 'r');
l = fgetl(codeId);
A = sscanf(l, '%f %f');
nCode = A(1);
nChip = A(2);
codes = zeros(nChip,nCode);
for i = 1:nCode
    l = fgetl(codeId);
    A = sscanf(l, [repmat('%f ', [1 nChip-1]) '%f']);
    codes(:,i) = A.';
end
fclose(codeId);

%% tx
txId = fopen(txFile, 'r');
if txId == -1
    disp(txFile);
end

% first line
l = fgetl(txId);
try
    A = sscanf(l, '%f %f %f %f %f');
    Lp = A(5);
catch
    A = sscanf(l, '%f %f %f %f');
    Lp = 4;
end
T = A(1) * 1e-3;
nBit = A(2);
nTx = A(3);
assert(nChip == A(4));

% tx info
txPumps = zeros(nTx,4);
switch code
    case "gold"
        xChip = zeros(nChip,nTx);
    otherwise
        xChip = zeros(nChip*2,nTx);
end
txOffset = zeros(1,nTx);
moOffset = zeros(1,nTx);
for i = 1:nTx
    l = fgetl(txId);
    A = sscanf(l, '(%f, %f, %f, %f), (code %f)');
    txPumps(i,:) = A(1:4).';
    txPumps(i,2) = A(2) * T;
    txOffset(1,i) = A(2);
    moOffset(1,i) = A(3);
    switch code
        case "gold"
            xChip(:,i)  =  codes(:,A(5)+1);
        otherwise
            xChip(1:2:end,i) = codes(:,A(5)+1);
            xChip(2:2:end,i) = 1 - codes(:,A(5)+1);
    end
end
fclose(txId);

% tx bits
txId = fopen(txFile, 'r');
txBits = zeros(nBit,nTx);
nBits = ones(1,nTx);
txChips = zeros((nBit+Lp)*nChip*2,nTx);
nChips = ones(1,nTx);

while true
    l = fgetl(txId);
    if contains(l,'START')
        break;
    end
end

startTime = zeros(1,nTx);
while true
    l = fgetl(txId);
    if l == -1
        break;
    end
    if contains(l,'bit')
        A = sscanf(l, '%f, %f, bit %f');
        txBits(nBits(A(2)+1),A(2)+1) = A(3);
        nBits(A(2)+1) = nBits(A(2)+1) + 1;
    elseif contains(l,'prem')
        A = sscanf(l, '%f, %f, prem %f');
        if nBits(A(2)+1) == 1
            startTime(A(2)+1) = A(1);
        end
    else
        A = sscanf(l, '%f, %f, %f');
        txChips(nChips(A(2)+1),A(2)+1) = A(3);
        nChips(A(2)+1) = nChips(A(2)+1) + 1;
    end
end
fclose(txId);

% txOffset
txOffset(2:end) = round((startTime(2:end)-startTime(1)) / (T*1e3)) ...
    + txOffset(1)+moOffset(1) - moOffset(2:end);

%%
tx.Lp = Lp;
tx.xChip = xChip;
tx.txOffset = txOffset;
tx.moOffset = moOffset;
tx.bits = txBits;
tx.chips = txChips;
end

function rxSignal = readrx(rxFile)
res = 1e-3;

rxId = fopen(rxFile, 'r');
A = textscan(rxId, '%f\t%f');
fclose(rxId);

A = cell2mat(A);
A(:,1) = (A(:,1) - A(1,1)) * 1e-3;
x = A(1,1):res:A(end,1);
rxSignal = [x; interp1(A(:,1), A(:,2), x)].';
end