linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ceName = "ce"+string(cenoteFinal);

%%
f = figure("Position", [100 100 600 300]);
box on; hold on; grid on;
algover = "11";
matfolder = "mat_"+matversion+"/mat1_"+algover+"/"+ceName;

alpha = 0.1;
totalNet = 0;
datarate = 2/1.75 * 100/116; % per Tx

algoPD = "sc";
algoName1 = "sc-af0";
algoName2 = "sc-af0";

domdma = true;
xRange = 1:4;
if domdma
    pertxput = nan(3,length(xRange));
else
    pertxput = nan(2,length(xRange));
end

%% MDMA
while domdma
    nTxRange = [1,2];
    T = 437;
    txNameRange = repmat("2",size(nTxRange));
    LpName = "2";
    Lp2Name = "";
    nMo = 1;
    osName = "";
    codeName = "plain0";
    
    goodrate = nan(size(nTxRange));
    for idx = 1:length(nTxRange)
        txName = txNameRange(idx);
        nTx = numel(strfind(txName,"-"))+1;
        if nTx > 1
            pertxput(1,idx) =  pertxput(1,1);
            continue;
        end
        algoName = algoName1;
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
        matName = "../"+matfolder+osName+"/"+preName+".mat";
        disp(matName);
        
        if isfile(matName)
            load(matName);
        else
            error("file not exist");
        end
        goodrate(idx) = mean(ber_temp<=alpha,'all');
        pertxput(1,idx) = goodrate(idx) * datarate;
    end

    break;
end

%% MCDMA
while true
    nTxRange = xRange;
    T = 125;
    txNameRange = ["","2","","2-3","","2-3-4"];
    LpName = "16";
    Lp2Name = "";
    nMo = 1;
    osName = "";
    codeName = "gold";
    
    goodrate = nan(size(nTxRange));
    for idx = [(1:floor(length(nTxRange)/2))*2,(1:ceil(length(nTxRange)/2))*2-1]
        if idx == 1
            pertxput(domdma+1,1) = pertxput(domdma+1,2);
            continue;
        elseif mod(idx,2) == 1
            nTx = nTxRange(idx);
            nTx1 = floor(nTx/2);
            nTx2 = ceil(nTx/2);
            pertxput(domdma+1,idx) = (pertxput(domdma+1,2*nTx1)*nTx1+pertxput(domdma+1,2*nTx2)*nTx2) ...
                / (nTx1+nTx2);
            continue;
        end
        
        txName = txNameRange(idx);
        nTx = numel(strfind(txName,"-"))+1;
        if nTx==1
            algoName = algoName1;
        else
            if algoPD == "gt"
                algoName = algoName1;
            else
                algoName = algoName2;
            end
        end
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
        matName = "../"+matfolder+osName+"/"+preName+".mat";
        disp(matName);
        
        if isfile(matName)
            load(matName);
        else
            error("file not exist");
        end
        goodrate(idx) = mean(ber_temp<=alpha,'all');
        pertxput(domdma+1,idx) = goodrate(idx) * datarate;
    end

    break;
end

%% MMCDMA
while true
    nTxRange = xRange;
    T = 125;
    txNameRange = ["2","3-4","2-3-4","2-3-4-5","2-3-4-5-6","2-3-4-5-6-7"];
    LpName = "16";
    Lp2Name = "";
    nMo = 2;
    osName = "";
    codeName = "goldman";
    
    goodrate = nan(size(nTxRange));
    for idx = 1:length(nTxRange)
        txName = txNameRange(idx);
        nTx = numel(strfind(txName,"-"))+1;
        if nTx==1
            algoName = algoName1;
        else
            if algoPD == "gt"
                algoName = algoName1;
            else
                algoName = algoName2;
            end
        end
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
        matName = "../"+matfolder+osName+"/"+preName+".mat";
        disp(matName);
        
        if isfile(matName)
            load(matName);
        else
            error("file not exist");
        end
        goodrate(idx) = mean(ber_temp<=alpha,'all');
        pertxput(domdma+2,idx) = goodrate(idx) * datarate;
    end

    break;
end

%%
b = bar(xRange, pertxput .* xRange.^totalNet);

xticks(xRange);
xlabel("Number of colliding TX");
if totalNet
    ylabel("Network throughput (bps)");
else
    ylabel("Per TX throughput (bps)");
end
set(f.Children(1), "FontSize", 10);

if domdma
    l = legend(["MDMA", "MDMA+CDMA", "MoMA"]);
else
    l = legend(["MDMA+CDMA", "MoMA"]);
end
if totalNet
    l.Location = "northwest";
else
    l.Location = "southwest";
end
set(f.Children(1), "FontSize", 10);
% xlim([gca().XLim(1),gca().XLim(2)-2]);

restyle(2);
if totalNet
%     saveas(f, "fig/figure_results_mainbar_nTx_thpt_total"+algover, "fig");
%     saveas(f, "jpg/figure_results_mainbar_nTx_thpt_total"+algover, "jpg");
else
%     saveas(f, "fig/figure_results_mainbar_nTx_thpt_perTx"+algover, "fig");
%     saveas(f, "jpg/figure_results_mainbar_nTx_thpt_perTx"+algover, "jpg");
end