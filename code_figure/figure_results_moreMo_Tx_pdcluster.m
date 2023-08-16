linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ftype = "hist";
topo = "lineMo";
switch topo
    case "line"
        foldernoteRange = ["1","5","1-3"];
        ceNameRange = repmat("os",size(foldernoteRange));
        txNameRange = repmat("2-3-4-5",size(ceNameRange));
        nMoRange = [2,2,2];
    case "lineMo"
        foldernoteRange = repmat("1",[1 6]);
        ceNameRange = repmat("ce21",size(foldernoteRange));
        txNameRange = ["2","3-4","2-3-4","2-3-4-5","2-3-4-5-6","2-3-4-5-6-7"];
        nMoRange = repmat(1,[1 6]);
    case "lineFork"
        foldernoteRange = repmat("4",[1 5]);
        ceNameRange = repmat("os",size(foldernoteRange));
        txNameRange = ["2-3","3-4","2-3-4","2-3-5","2-3-4-5"];
        nMoRange = repmat(2,[1 5]);
    case "lineMix"
        foldernoteRange = repmat("1-3",[1 4]);
        ceNameRange = repmat("os",size(foldernoteRange));
        txNameRange = ["2","3-4","2-3-4","2-3-4-5"];
        nMoRange = repmat(2,[1 4]);
end

%%
algover = "11";

alpha = 0.1;
totalNet = 0;
datarate = 2/1.75 * 100/116; % per Tx

titleName = strings(length(ceNameRange),1);

%%
f = figure("Position", [100 100 1200 300]);
tiledlayout(size(ceNameRange,1),size(ceNameRange,2), ...
    "TileSpacing","compact", "Padding","compact");

for idx = 1:length(ceNameRange)
    matfolder = "mat_"+matversion+"/mat"+foldernoteRange(idx)+"PDdebug_"+algover;
    ceName = ceNameRange(idx);    
    txName = txNameRange(idx);
    nMo = nMoRange(idx);

    switch foldernoteRange(idx)
        case "1"
            titleName(idx) = strcat(repmat('salt-',[1,nMo-1]))+"salt";
        case "5"
            titleName(idx) = strcat(repmat('soda-',[1,nMo-1]))+"soda";
        case "1-3"
            titleName(idx) = "salt-soda";
    end

    T = 125;
    LpName = "16";
    Lp2Name = "";
    codeName = "goldman";
    algoName = "sc-af0";
    
    preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
        +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
    
    if isfile("../"+matfolder+"/"+ceName+"/"+preName+".mat")
        load("../"+matfolder+"/"+ceName+"/"+preName+".mat");
    else
        error("file not exist");
    end

    posTcorr = cell(size(pddebug_temp,1));
    posTratio = cell(size(pddebug_temp,1));
    for ii = 1:size(pddebug_temp,1)
        indices = find(pddebug_temp{ii}.labels(:,1)==1);
        posTcorr{ii} = pddebug_temp{ii}.corr(indices,end,:);
        posTratio{ii} = pddebug_temp{ii}.ratio(indices,:,:);
        posTratio{ii}(posTratio{ii}>1) = 1./posTratio{ii}(posTratio{ii}>1);
        posTratio{ii} = mean(posTratio{ii},2);
    end
    posFcorr = cell(size(pddebug_temp,1),1);
    posFratio = cell(size(pddebug_temp,1),1);
    for ii = 1:size(pddebug_temp,1)
        indices = find(pddebug_temp{ii}.labels(:,1)==0);
        posFcorr{ii} = pddebug_temp{ii}.corr(indices,end,:);
        posFratio{ii} = pddebug_temp{ii}.ratio(indices,:,:);
        posFratio{ii}(posFratio{ii}>1) = 1./posFratio{ii}(posFratio{ii}>1);
        posFratio{ii} = mean(posFratio{ii},2);
    end

    posTcorr = min(squeeze(cell2mat(posTcorr)),[],2,"omitnan");
    posTratio = min(squeeze(cell2mat(posTratio)),[],2,"omitnan");
    posFcorr = min(squeeze(cell2mat(posFcorr)),[],2,"omitnan");
    posFratio = min(squeeze(cell2mat(posFratio)),[],2,"omitnan");

    svmX = [[posTcorr;posFcorr],[posTratio;posFratio]];
    svmy = [repmat("true pos",[length(posTcorr),1]);repmat("false pos",[length(posFcorr),1])];

    nexttile;
    switch ftype
        case "gscatter"
            gscatter(svmX(:,1),svmX(:,2),svmy,[],"ox");
            legend off; xlim([0 1]); ylim([0 1]);
            xlabel("correlation");
            if idx == 1
                ylabel("power ratio");
                legend(["true pos","false pos"],"location","west");
            end
        case "hist"
            histRange = [0 1]; histInterval = 0.05;
            histEdges = histRange(1):histInterval:histRange(2);
            hhT = hist3([posTcorr,posTratio],"Edges",{histEdges histEdges});
            hhF = hist3([posFcorr,posFratio],"Edges",{histEdges histEdges});
            hh = (hhT-hhF)./(hhT+hhF);
%             hh = hhT-hhF;
            hh(isnan(hh)) = 0;
            hobj = pcolor(histEdges,histEdges,-hh.');
            hobj.EdgeColor = "None";
            axis on; axis equal; axis tight; axis xy;
            set(gca,"YDir","normal");
            clim([-1 1]);
%             clim([-max(abs(hh),[],"all"),max(abs(hh),[],"all")]);
            colormap("jet");
            xlabel("correlation");
            if idx == 1
                ylabel("power ratio");
            end
    end
    title(titleName(idx));

%         if ~isempty(posTcorr) && ~isempty(posFcorr)
%             svmModel = fitcsvm(svmX,svmy);
%             svmLinex = linspace(0,1,101);
%             svmLiney = -(svmModel.Beta(1)/svmModel.Beta(2)*svmLinex)-svmModel.Bias/svmModel.Beta(2);
%             plot(svmLinex,svmLiney,"k-");
%         end

end

% l = legend(["true pos", "false pos"]);
% l.Layout.Tile = "east";

%%

restyle(4);
% saveas(f, "fig/figure_results_moreMo_nTx_pdcluster_"+topo+algover, "fig");
% saveas(f, "jpg/figure_results_moreMo_nTx_pdcluster_"+topo+algover, "jpg");