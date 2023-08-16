linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
f = figure("Position", [100 100 600 300]);
box on; hold on; grid on;
algover = "11";
matfolder = "mat_"+matversion+"/mat1_"+algover;
algoName = "sc-af0";

alpha = 0.1;
totalNet = 1;
datarate = 1/1.75;

%% n molecule(s)
for moName = ["1"]
    T = 125;
    txName = "2-3-4-5"; nTx = numel(strfind(txName,"-"))+1;
    LpRange = [4,8,16,32];
    Lp2Name = "";
    osName = "";
    codeName = "goldman"; codelength = 14;
    
    goodrate = nan(size(LpRange));
    for idx = 1:length(LpRange)
        Lp = LpRange(idx);
        LpName = string(Lp);
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+moName+"_"+algoName;
        
        matname = "../"+matfolder+"/os"+osName+"/"+preName+".mat";
        if isfile(matname)
            load(matname);
        else
            error("file not exist");
        end
        goodrate(idx) = mean(ber_temp<=alpha,'all');
    end
    plot(LpRange, goodrate .* datarate .* (100/(100+Lp)) .* nTx.^totalNet, '-x', ...
        'LineWidth', linewidth_default, 'MarkerSize', markersize_default);
end

%%
xlabel("preamble length (symbols)");
ylabel("network throughput (bps)");
% set(f.Children(1), "YScale", "log");
set(f.Children(1), "FontSize", 10);

% legend(["1 molecule", "2 molecules", "3 molecules"], ...
%     "Location", "southeast");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_moreMo_Lp_thpt"+algover, "fig");
% saveas(f, "jpg/figure_results_moreMo_Lp_thpt"+algover, "jpg");