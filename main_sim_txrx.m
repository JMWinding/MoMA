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
T = 125;
T2 = T;
Lp = 16;
Lp2 = Lp;
code = "goldman";
nDegree = 3;

%% get ground truth CIR
constructIn = struct(...
    "T", T, ...
    "Lp", Lp, ...
    "code", code, ...
    "nDegree", nDegree, ...
    "T2", T2, ...
    "Lp2", Lp2);
try constructIn.hPre = hPre; catch, constructIn.hPre = ceil(1250/T2); end
try constructIn.hPost = hPost; catch, constructIn.hPost = ceil(1750/T2); end
try constructIn.hlen = hlen; catch ; end
rxIn = sim_mmo_tx(constructIn);

rxIn.weights_ce = GetCEWeights(221);
rxIn.sameMo = true;

%% rxOut
% rxIn.algoPD = "gt";
% rxIn.algoCE = "gt";
% rxIn.mode = "dc";
% rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
% disp(cell2mat(rxOut.BER));
% disp("----");

% rxIn.algoPD = "gt";
% rxIn.algoCE = "af0";
% rxIn.mode = "dc";
% rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
% disp(cell2mat(rxOut.BER));
% disp("----");

rxIn.algoPD = "gt1";
rxIn.algoCE = "af0";
rxIn.mode = "dc";
rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
disp(rxOut.PDOff);
disp(cell2mat(rxOut.BER));
disp("----");

% rxIn.algoPD = "sc";
% rxIn.algoCE = "af0";
% rxIn.debug_pd = false;
% rxIn.mode = "dc";
% rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
% disp(cell2mat(rxOut.BER));
% disp("----");
disp("------------------------------");

%%
if false
    %%
    figure("units","normalized","outerposition",[0 0 1 1]);
    for j = 1:size(txOut.xBit,1)
    for i = 1:nTx
        subplot(size(txOut.xBit,1), nTx, (j-1)*size(txOut.xBit,1)+i);
        hold on; box on; grid on;
        stem(txOut.xBit{j,i});
        stem(rxOut.dBit{j,i});
    %     yyaxis right;
    %     plot(txOut.xChannel.xCIR{j,i}, "k");
    %     plot(rxOut.chan.hp{j,i}, "k");
    end
    end
end