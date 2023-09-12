linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ceName = "ce301";

%%
f = figure("Position", [100 100 600 300]);
box on; hold on; grid on;
algoName = "sc-af0";
algover = "11";

%% n molecule(s)
for moName = ["1","2"] %,"3"]
while false
    TRange = [100,125,150,175,200,250];
    txName = "2-3-4-5"; nTx = numel(strfind(txName,"-"))+1;
    LpName = "4";
    Lp2Name = "_Lp16";
    osName = "";
    codeName = "goldman"; codelength = 14;
    
    hit = nan(size(TRange));
    for TIdx = 1:length(TRange)
        T = TRange(TIdx);
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+moName+"_"+algoName;        
        matName = "../mat_"+matversion+"/mat3/"+ceName+osName+"/"+preName+".mat";

        if isfile(matName)
            disp(matName);
            load(matName);
        else
            error("file not exist");
        end
        hit(TIdx) = mean(sum(pdoff_temp>ceil(ceil(1e3/T)/2),[2,3])==0);
    end
    plot(nTx./(TRange*codelength/1e3), hit, '-x', ...
        'LineWidth', linewidth_default, 'MarkerSize', markersize_default);

    break;
end

while true
    TRange = [250,200,175,125];
    txName = "2-3-4-5"; nTx = numel(strfind(txName,"-"))+1;
    LpName = "16";
    Lp2Name = "";
    osName = "";
    codeName = "goldman"; codelength = 14;
    
    hit = nan(size(TRange));
    for TIdx = 1:length(TRange)
        T = TRange(TIdx);
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+moName+"_"+algoName;        
        matName = "../mat_"+matversion+"/mat1_11/ce301"+osName+"/"+preName+".mat";

        if isfile(matName)
            disp(matName);
            load(matName);
        else
            error("file not exist");
        end
        hit(TIdx) = mean(sum(pdoff_temp>ceil(ceil(1e3/T)/2),[2,3])==0);
    end
    plot(nTx./(TRange*codelength/1e3), hit, '-x', ...
        'LineWidth', linewidth_default, 'MarkerSize', markersize_default);

    break;
end
end

%%
xlabel("Datarate per molecule (bps)");
ylabel("Correct detection rate");
set(f.Children(1), "FontSize", 10);

legend(["1 molecule", "2 molecules"], ..., "3 molecules"], ...
    "Location", "southeast");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_moreMo_rate_pd"+algover, "fig");
% saveas(f, "jpg/figure_results_moreMo_rate_pd"+algover, "jpg");