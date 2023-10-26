linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ceName = "ce"+string(cenoteFinal);

%%
T = 125;
LpName = "16";
Lp2Name = "";
osName = "";
codeName = "goldman"; codelength = 14;
algoName = "gt1-af0";
algover = "11";
topo = "line";

%%
switch topo
    case "line"
        fnotes = ["1", "3"];
        txNameRange = ["2","3-4","2-3-4","2-3-4-5"];
    case "fork"
        fnotes = ["4", "5"];
        txNameRange = ["2-3","3-4","2-3-4","2-3-5","2-3-4-5"];
end

%% n molecule(s)
legendText = ["salt-1","salt-2","salt-mix","soda-1","soda-2","soda-mix"];
noteCombRange = {[fnotes(1),"1"]; [fnotes(1),"2"]; ...
    [fnotes(1)+"-"+fnotes(2),"2",1]; ...
    [fnotes(2),"1"]; [fnotes(2),"2"]; ...
    [fnotes(1)+"-"+fnotes(2),"2",2]};

% legendText = ["salt-1","salt-mix","soda-1","soda-mix"];
% noteCombRange = {[fnotes(1),"1"]; ...
%     [fnotes(1)+"-"+fnotes(2),"2",1]; ...
%     [fnotes(2),"1"]; ...
%     [fnotes(1)+"-"+fnotes(2),"2",2]};

berplot = nan(length(noteCombRange),length(txNameRange));
for idx1 = 1:length(noteCombRange)
    noteComb = noteCombRange{idx1};
    matfolder = "mat_"+matversion+"/mat"+noteComb(1)+"_"+algover;
    moName = noteComb(2);
    
    for idx2 = 1:length(txNameRange)
        txName = txNameRange(idx2);
        
        preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
            +"_"+codeName+Lp2Name+"_"+moName+"_"+algoName;
        matName = "../"+matfolder+"/"+ceName+osName+"/"+preName+".mat";
        disp(matName);

        if isfile(matName)
            load(matName);
        else
            error("file not exist");
        end

        if length(noteComb) == 2
            berplot(idx1,idx2) = mean(ber_temp,"all");
        else
            berplot(idx1,idx2) = mean(ber_temp(:,str2double(noteComb(3)),:),"all");
        end
    end
end

%%
f = figure("Position", [100 100 600 300]);
box on; hold on; grid on;
bar(berplot.');

% title(topo);
xlabel("active Tx");
xticks(1:length(txNameRange));
xticklabels(replace(txNameRange,["2";"3";"4";"5"],["1";"2";"3";"4"]));
ylabel("mean BER");
set(f.Children(1), "FontSize", 10);

legend(legendText, "Location","northwest");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_moreMo_Tx_ce_"+topo+algover, "fig");
% saveas(f, "jpg/figure_results_moreMo_Tx_ce_"+topo+algover, "jpg");