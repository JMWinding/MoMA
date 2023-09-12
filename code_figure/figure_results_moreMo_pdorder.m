linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ceName = "ce301";

%%
T = 125;
txName = "2-3-4-5"; nTx = numel(strfind(txName,"-"))+1;
LpName = "16";
Lp2Name = ""; % "", "_Lp16"
osName = ""; % "", "_2x"
codeName = "goldman"; codelength = 14;
algoName = "sc-af0";
algover = "11";
topo = "fork";

%%
switch topo
    case "line"
        fnotes = ["1", "5"];
    case "fork"
        fnotes = ["4", "5"];
end

%%
% legendText = ["salt-1","salt-2","soda-1","soda-2","mix"];
% noteCombRange = { ...
%     ["mat"+fnotes(1),"1"]; ...
%     ["mat"+fnotes(1),"2"]; ...
%     ["mat"+fnotes(2),"1"]; ...
%     ["mat"+fnotes(2),"2"]; ...
%     ["mat"+fnotes(1)+"-"+fnotes(2),"2"]};

legendText = ["1 molecule","2 molecules"];
noteCombRange = { ...
    ["mat"+fnotes(1),"1"]; ...
    ["mat"+fnotes(1),"2"]};

hitplot = zeros(length(noteCombRange),nTx);
for idx1 = 1:length(noteCombRange)
    noteComb = noteCombRange{idx1};
    noteName = noteComb(1);
    moName = noteComb(2);

    preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
        +"_"+codeName+Lp2Name+"_"+moName+"_"+algoName;        
    matName = "../mat_"+matversion+"/"+noteName+"_11/"+ceName+osName+"/"+preName+".mat";

        if isfile(matName)
            disp(matName);
            load(matName);
        else
            error("file not exist");
        end
    
    for i = 1:nTx
        hitplot(idx1,i) = mean(abs(pdoff_temp(:,:,i))<=ceil(ceil(1e3/T)/2));
    end
end

%%
f = figure('Position', [100 100 600 300]);
box on; hold on; grid on;
bar(hitplot.');

xticks(1:nTx);
xlabel('Tx arriving order');
ylabel('correct detection rate');
% title("1/"+num2str(T*14/1e3)+" bps");
set(f.Children(1), "FontSize", 10);

legend(legendText, "Location","southwest");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_moreMo_pdorder_"+string(T)+"ms"+algover, "fig");
% saveas(f, "jpg/figure_results_moreMo_pdorder_"+string(T)+"ms"+algover, "jpg");