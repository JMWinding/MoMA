% DEBUG CODE 1
% for ijk = 1:size(notes_temp,1)
%     disp(ijk);
%     disp(squeeze(ber_temp(ijk,:,:)));
%     notes = arrayfun(@(a)(string(sprintf("%02d",a))),notes_temp(ijk,:).');
%     main_emulate_txrx
% end

%%
addpath("code_algo");
warning("off");

ismc = true;
isrepeat = false;
noisemodel = "pois";
isfigure = false;

%%
datanote = "dataset/data1";
T = 125;
T2 = T;
Lp = 16;
Lp2 = Lp;
pumpstr = "3-4";
code = "goldman";
if ~exist("notes", "var")
notes = arrayfun(@(a)(string(sprintf("%02d",a))),randi([1 40],[2 1]));
end
disp(notes);

%% get ground truth CIR
constructIn = struct(...
    "datanote", datanote, ...
    "T", T, ...
    "pumpstr", pumpstr, ...
    "Lp", Lp, ...
    "code", code, ...
    "notes", notes, ...
    "T2", T2, ...
    "Lp2", Lp2);
try constructIn.hPre = hPre; catch, constructIn.hPre = ceil(1250/T2); end
try constructIn.hPost = hPost; catch, constructIn.hPost = ceil(1750/T2); end
try constructIn.hlen = hlen; catch ; end
rxIn = emulates_construct_rxIn(constructIn);

rxIn.weights_ce = GetCEWeights(221);
rxIn.sameMo = length(datanote)==1;

%% rxOut
rxIn.algoPD = "gt";
rxIn.algoCE = "gt";
rxIn.mode = "dc";
rxOut = decode_mmo_coherent_MMoNTxSW11loop(rxIn);
disp(cell2mat(rxOut.BER));
disp("----");

% rxIn.algoPD = "gt";
% rxIn.algoCE = "af0";
% rxIn.mode = "dc";
% rxOut = decode_mmo_coherent_MMoNTxSW11loop(rxIn);
% disp(cell2mat(rxOut.BER));
% disp("----");

% rxIn.algoPD = "gt1";
% rxIn.algoCE = "af0";
% rxIn.mode = "dc";
% rxOut = decode_mmo_coherent_MMoNTxSW11loop(rxIn);
% disp(rxOut.PDOff);
% disp(cell2mat(rxOut.BER));
% disp("----");

% rxIn.algoPD = "sc";
% rxIn.algoCE = "af0";
% rxIn.debug_pd = false;
% rxIn.mode = "dc";
% rxOut = decode_mmo_coherent_MMoNTxSW11loop(rxIn);
% disp(cell2mat(rxOut.BER));
% disp("----");
disp("------------------------------");

%%
if false
    %%
    figure("units","normalized","outerposition",[0 0 1 1]);
    for j = 1:size(rxIn.xBit,1)
    for i = 1:size(rxIn.xBit,2)
        subplot(size(rxIn.xBit,1), size(rxIn.xBit,2), (j-1)*size(rxIn.xBit,1)+i);
        hold on; box on; grid on;
        stem(rxIn.xBit{j,i});
        stem(rxOut.dBit{j,i});
        yyaxis right;
        plot(rxIn.xChannel.xCIR{j,i}, "-ok");
        plot(rxOut.chan.hp{j,i}, "-xk");
    end
    end
end