linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ftype = "hist";
mol = "salt";
topo = "line";
switch topo
    case "line"
        ceNameRange = ["ce301","ce305","ce312"];
        txNameRange = repmat("2-3-4-5",size(ceNameRange));
        nMoRange = ones(size(ceNameRange));
        switch mol
            case "salt"
                foldernote = "1";
            case "soda"
                foldernote = "5";
        end
    case "fork"
        ceNameRange = ["ce301","ce305","ce312"];
        txNameRange = repmat("2-3-4-5",size(ceNameRange));
        nMoRange = ones(size(ceNameRange));
        switch mol
            case "salt"
                foldernote = "4";
            case "soda"
                foldernote = "5";
        end
end

%%
algover = "11";
matfolder = "mat_"+matversion+"/mat"+foldernote+"PDdebug_"+algover;

alpha = 0.1;
totalNet = 0;
datarate = 2/1.75 * 100/116; % per Tx

titleName = strings(length(ceNameRange),1);

%%
f = figure("Position", [100 100 1200 300]);
tiledlayout(size(ceNameRange,1),size(ceNameRange,2), ...
    "TileSpacing","compact", "Padding","compact");

for idx = 1:length(ceNameRange)
    ceName = ceNameRange(idx);    
    txName = txNameRange(idx);
    nMo = nMoRange(idx);
    
    % L0 LS
    % L1 non-negativity
    % L2 head-tail loss
    % L3 similarity (not applicable for 1 molecule)
    % L4 smoothness
    switch str2num(erase(ceName,"ce"))
        case 0
            titleName(idx) = "L0+L2+L4";
        case 1
            error("");
            titleName(idx) = "L0+L1+L2+L4";
        case 2
            titleName(idx) = "L0+L1+L2";
        case 3
            titleName(idx) = "L0+L1+L4";
        case 4
            titleName(idx) = "L0+L1";
        case 5
            error("");
            titleName(idx) = "L0+L1+L4";
        case 6
            titleName(idx) = "L0+L4";
        case 7
            titleName(idx) = "L0+L1+L2+L4";
        case 8
            titleName(idx) = "L0+L2";
        case 9
            titleName(idx) = "L0+L1+L2(-L3)";
        case 10
            titleName(idx) = "L0+L1(-L3)";
        case 11
            titleName(idx) = "L0+L2(-L3)";
    end

    T = 125;
    LpName = "16";
    Lp2Name = "";
    codeName = "goldman";
    algoName = "sc-af0";
    
    preName = "emulates_"+num2str(T)+"ms_"+txName+"_"+LpName ...
        +"_"+codeName+Lp2Name+"_"+string(nMo)+"_"+algoName;
    matName = "../"+matfolder+"/"+ceName+"/"+preName+".mat";
    
    if isfile(matName)
        disp(matName);
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
% saveas(f, "fig/figure_results_ceLoss_nTx_pdcluster_"+topo+"-"+mol+algover, "fig");
% saveas(f, "jpg/figure_results_ceLoss_nTx_pdcluster_"+topo+"-"+mol+algover, "jpg");