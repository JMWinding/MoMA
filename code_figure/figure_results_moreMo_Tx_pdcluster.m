linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ftype = "hist";
foldernoteRange = ["1","1","1"];
ceNameRange = repmat("ce"+string(cenoteFinal),size(foldernoteRange));
txNameRange = ["2-3-4-5","2-3-4-5","2-3"];
nMoRange = [1,2,1];
codeNameRange = ["goldman","goldman","gold"];

%%
algover = "11";

alpha = 0.1;
totalNet = 0;
datarate = 2/1.75 * 100/116; % per Tx

thrdcorrRange = 0.5:0.01:0.9;
thrdratioRange = 0:0.01:0.8;

titleName = strings(length(ceNameRange),1);

%%
f = figure("Position", [100 100 1200 300], "Name", ceNameRange(1));
tiledlayout(size(ceNameRange,1),size(ceNameRange,2), ...
    "TileSpacing","compact", "Padding","compact");

for idx = 1:length(ceNameRange)
    matfolder = "mat_"+matversion+"/mat"+foldernoteRange(idx)+"PDdebug_"+algover;
    ceName = ceNameRange(idx);    
    txName = txNameRange(idx);
    nMo = nMoRange(idx);
    codeName = codeNameRange(idx);

    switch foldernoteRange(idx)
        case "1"
            titleName(idx) = strcat(repmat('salt-',[1,nMo-1]))+"salt";
        case "3"
            titleName(idx) = strcat(repmat('soda-',[1,nMo-1]))+"soda";
        case "1-3"
            titleName(idx) = "salt-soda";
    end

    T = 125;
    LpName = "16";
    Lp2Name = "";
    algoName = "sc-af0";
    
    preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
        +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
    matName = "../"+matfolder+"/"+ceName+"/"+preName+".mat";
    disp(matName);
    
    if isfile(matName)
        load(matName);
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

    %%
    falseNeg = zeros(length(thrdcorrRange),length(thrdratioRange));
    falsePos = zeros(length(thrdcorrRange),length(thrdratioRange));
    for thrdcorrIdx = 1:length(thrdcorrRange)
        thrdcorr = thrdcorrRange(thrdcorrIdx);
        for thrdratioIdx = 1:length(thrdratioRange)
            thrdratio = thrdratioRange(thrdratioIdx);
            falseNeg(thrdcorrIdx,thrdratioIdx) = 1-mean(posTcorr>=thrdcorr & posTratio>=thrdratio);
            falsePos(thrdcorrIdx,thrdratioIdx) = mean(posFcorr>=thrdcorr & posFratio>=thrdratio);
        end
    end
    falseAll = falseNeg + falsePos;
    [thrdcorrIdx,thrdratioIdx] = find(falseAll == min(falseAll,[],"all"), 1);
    title([thrdcorrRange(thrdcorrIdx), thrdratioRange(thrdratioIdx), min(falseAll,[],"all")]);

    %%
    if nMo == 1
        falseAll1 = falseAll;
    elseif nMo == 2
        falseAll2 = falseAll;
        falseAll = falseAll1 + falseAll2;

        [thrdcorrIdx,thrdratioIdx] = find(falseAll == min(falseAll,[],"all"), 1);
        disp([thrdcorrRange(thrdcorrIdx), thrdratioRange(thrdratioIdx), ...
            falseAll1(thrdcorrIdx,thrdratioIdx), falseAll2(thrdcorrIdx,thrdratioIdx)]);

        [thrdcorrIdx,thrdratioIdx] = find(falseAll == min(falseAll,[],"all"), 1, "last");
        disp([thrdcorrRange(thrdcorrIdx), thrdratioRange(thrdratioIdx), ...
            falseAll1(thrdcorrIdx,thrdratioIdx), falseAll2(thrdcorrIdx,thrdratioIdx)]);

        thrdcorrIdx = find(abs(thrdcorrRange-0.8)<1e-4);
        thrdratioIdx = find(abs(thrdratioRange-0.2)<1e-4);
        disp([thrdcorrRange(thrdcorrIdx), thrdratioRange(thrdratioIdx), ...
            falseAll1(thrdcorrIdx,thrdratioIdx), falseAll2(thrdcorrIdx,thrdratioIdx)]);
    end

end

% l = legend(["true pos", "false pos"]);
% l.Layout.Tile = "east";

%%

restyle(4);
% saveas(f, "fig/figure_results_moreMo_nTx_pdcluster_"+topo+algover, "fig");
% saveas(f, "jpg/figure_results_moreMo_nTx_pdcluster_"+topo+algover, "jpg");