linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
f = figure("Position", [100 100 600 300]);
box on; hold on; grid on;
algover = "11";
matfolder = "mat_"+matversion+"/mat1_"+algover;

%% MMCDMA
nTxRange = [2,3,4];
TRange = [125,125,125];
txNameRange = ["3-4","2-3-4","2-3-4-5"];
LpName = "16";
Lp2Name = "";
nMo = 1; moName = string(nMo);
osName = "";
codeName = "goldman";

berplot = nan(2,length(TRange));
for idx = 1:length(TRange)
    T = TRange(idx);
    txName = txNameRange(idx);
    algoName = "sc-af0";
    
    preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
        +"_"+codeName+Lp2Name+"_"+moName+"_"+algoName;
    
    matname = "../"+matfolder+"/os"+osName+"/"+preName+".mat";
    if isfile(matname)
        load(matname);
    else
        error("file not exist");
    end

    idxgood = pdoff_temp<=ceil(ceil(1e3/T)/2);
    idxallgood = sum(~idxgood,[2,3])==0;
    ber_temp(~repmat(idxgood,[1,nMo,1])) = nan;
    berplot(1,idx) = median(ber_temp(idxallgood,:),'all','omitnan');
    berplot(2,idx) = median(ber_temp(~idxallgood,:),'all','omitnan');
end

%%
bar(nTxRange, berplot);
xticks(nTxRange);
xlabel("number of colliding TX");
ylabel("median BER");
set(f.Children(1), "FontSize", 10);

l = legend(["detected", "missing"], "Location", "southeast");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_missing"+algover, "fig");
% saveas(f, "jpg/figure_results_missing"+algover, "jpg");