% DEBUG CODE 1
% for ijk = 1:size(notes_temp,1)
%     disp(ijk);
%     disp(squeeze(ber_temp(ijk,:,:)));
%     notes = arrayfun(@(a)(string(sprintf("%02d",a))),notes_temp(ijk,:).');
%     main_emulate_txrx
% end

%%
addpath("code_algo");
warning('off');

ismc = true;
isrepeat = false;
noisemodel = 'pois';
isfigure = false;

%%
datanote = "dataset/data1";
T = 125;
Lp = 16; Lp2 = Lp;
pumpstr = "2-3-4-5";
code = "goldman";
notes = arrayfun(@(a)(string(sprintf("%02d",a))),randi([1 40],[2 1]));
disp(notes);

% T = 437;
% Lp = 2; Lp2 = Lp;
% pumpstr = "2";
% code = "plain0";
% notes = "01";

%% get ground truth CIR
rxIn = emulates_construct_rxIn(struct(...
    'datanote', datanote, ...
    'T', T, ...
    'pumpstr', pumpstr, ...
    'Lp', Lp, ...
    'code', code, ...
    'notes', notes, ...
    'T2', T, ...
    'Lp2', Lp2, ...
    'hPre', ceil(1250/T), ...
    'hPost', ceil(1750/T)));
try rxIn.hlen = hlen; catch; end

%% rxOut
rxIn.noisemodel = 'pois';

% rxIn.algoPD = 'gt';
% rxIn.algoCE = 'gt';
% rxIn.mode = 'dc';
% rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
% disp(cell2mat(rxOut.BER));
% disp('----');

% rxIn.algoPD = 'gt';
% rxIn.algoCE = 'af0';
% rxIn.mode = 'dc';
% rxIn.sameMo = true;
% rxIn.weight_ce = GetCEWeights(301);
% for i = 1:length(rxIn.txOffset)
%     rxIn.txOffset{i} = rxIn.txOffset{i} + pdoff_temp(ijk,1,i);
% end
% % rxIn.debug_ce = true;
% rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
% disp(cell2mat(rxOut.BER));
% disp('----');

rxIn.debug_pd = false;
rxIn.algoPD = 'sc';
rxIn.algoCE = 'af0';
rxIn.mode = 'dc';
rxIn.sameMo = true;
rxOut = decode_mmo_coherent_MMoNTxSW11ce(rxIn);
disp(cell2mat(rxOut.BER));
disp('----');
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