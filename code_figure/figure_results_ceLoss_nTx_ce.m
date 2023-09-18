linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ceRange = [1,5,12]+300;

%%
mol = "salt";
topo = "line";
switch topo
    case "line"
        txNameRange = ["2","3-4","2-3-4","2-3-4-5"];
        switch mol
            case "salt"
                foldernote = "1";
            case "soda"
                foldernote = "3";
        end
    case "fork"
        txNameRange = ["2-3","3-4","2-3-4","2-3-5","2-3-4-5"];
        switch mol
            case "salt"
                foldernote = "4";
            case "soda"
                foldernote = "5";
        end
end

%%
f = figure("Position", [100 100 600 300]);
box on; hold on; grid on;
algover = "11";
matfolder = "mat_"+matversion+"/mat"+foldernote+"_"+algover;

alpha = 0.1;
totalNet = 0;
datarate = 2/1.75 * 100/116; % per Tx

legendName = strings(length(ceRange),1);
berAll = nan(length(ceRange), length(txNameRange));

%% MMCDMA
for ceIdx = 1:length(ceRange)
    ceName = "ce"+string(ceRange(ceIdx));
    % L0 LS
    % L1 non-negativity
    % L2 head-tail loss
    % L3 similarity (not applicable for 1 molecule)
    % L4 smoothness
    switch mod(ceRange(ceIdx),100)
        case 1
            legendName(ceIdx) = "L0+L1+L2";
        case 5
            legendName(ceIdx) = "L0+L1";
        case 12
            legendName(ceIdx) = "L0+L2";
    end

    T = 125;
    LpName = "16";
    Lp2Name = "";
    nMo = 1;
    codeName = "goldman";
    algoName = "gt-af0";
    
    for idx = 1:length(txNameRange)
        txName = txNameRange(idx);
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
        matName = "../"+matfolder+"/"+ceName+"/"+preName+".mat";
        disp(matName);
        
        if isfile(matName)
            load(matName);
        else
            error("file not exist");
        end
        berAll(ceIdx,idx) = mean(ber_temp, "all");
%         goodrate(idx) = mean(ber_temp<=alpha,'all');
%         pertxput(domdma+2,idx) = goodrate(idx) * datarate;
    end
end

%%
b = bar(1:length(txNameRange), berAll.');

% title(topo+" "+mol);
xticks(1:length(txNameRange));
xticklabels(replace(txNameRange,["2";"3";"4";"5"],["1";"2";"3";"4"]));
xlabel("active Tx");
ylabel("mean BER");
set(f.Children(1), "FontSize", 10);

legend(legendName, "Location", "northwest");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_ceLoss_nTx_ce_"+topo+"-"+mol+algover, "fig");
% saveas(f, "jpg/figure_results_ceLoss_nTx_ce_"+topo+"-"+mol+algover, "jpg");