% for notes = arrayfun(@(a)sprintf("%02d",a), 0:20)
%     disp(notes);
%     main_emulate_txrx;
% end

addpath("code_algo");
warning('off');

ismc = true;
isrepeat = false;
noisemodel = 'pois';
isfigure = false;

%%
datanote = "dataset/data1";
T = 125;
Lp = 16;
Lp2 = Lp;
pumpstr = "2-3-4-5";
code = "goldman";
notes = "03";

% type = "ph";
% hPre = ceil(3e3/T);
% hPost = ceil(7e3/T);
% hlen = ceil(15e3/T);

%% get ground truth CIR
rxIn = emulates_construct_rxIn(struct(...
    'datanote', datanote, ...
    'T', T, ...
    'pumpstr', pumpstr, ...
    'Lp', Lp, ...
    'code', code, ...
    'notes', notes, ...
    'T2', T, ...
    'Lp2', Lp2));
try rxIn.type = type; catch ; end
try rxIn.hPre = hPre; catch ; end
try rxIn.hPost = hPost; catch ; end
try rxIn.hlen = hlen; catch; end

%% rxOut
rxIn.noisemodel = 'pois';

% rxIn.algoPD = 'gt';
% rxIn.algoCE = 'gt';
% rxIn.mode = 'dc';
% rxOut = decode_mmo_coherent_MMoNTxSW11(rxIn);
% disp(rxOut.BER);
% disp('----');

rxIn.algoPD = 'gt';
rxIn.algoCE = 'af0';
rxIn.mode = 'dc';
rxIn.sameMo = true;
rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
disp(rxOut.BER);
rxIn.txOffset = cellfun(@(a)(a+randi([0 2])),rxIn.txOffset,"un",0);
rxIn.weights_ce  = struct("pos", 0.1, "posy", 0, "simTx", 0.1, "simMo", 0, "smth", 0, "cntr", 0.1);
rxIn.sameMo = false;
rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
disp(rxOut.BER);
rxIn.sameMo = true;
rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
disp(rxOut.BER);
% disp('----');

% rxIn.debug_pd = false;
% rxIn.algoPD = 'sc';
% rxIn.algoCE = 'af0';
% rxIn.mode = 'dc';
% rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
% disp(rxOut.BER);
% disp('----');
disp('------------------------------');

%%
if false
    %%
    figure('units','normalized','outerposition',[0 0 1 1]);
    for j = 1:size(txOut.xBit,1)
    for i = 1:nTx
        subplot(size(txOut.xBit,1), nTx, (j-1)*size(txOut.xBit,1)+i);
        hold on; box on; grid on;
        stem(txOut.xBit{j,i});
        stem(rxOut.dBit{j,i});
    %     yyaxis right;
    %     plot(txOut.xChannel.xCIR{j,i}, 'k');
    %     plot(rxOut.chan.hp{j,i}, 'k');
    end
    end
end