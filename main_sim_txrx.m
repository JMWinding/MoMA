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
    "nTx", 1, ...
    "code", code, ...
    "nDegree", nDegree, ...
    "T2", T2, ...
    "Lp2", Lp2);
rxIn = sim_mmo_tx(constructIn);

rxIn.weights_ce = GetCEWeights(221);
rxIn.sameMo = true;

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