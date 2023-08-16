%%
expdur = hours(24);
starttime = datetime;
predur = hours(0);
if isempty(gcp('nocreate'))
    parpool(feature('numcores'));
end

%%
cenote0 = 101;
for cenote = [1,5,12]+100
for pdnote = cenote
    weights_ce = GetCEWeights(cenote);
    weights_pd = GetCEWeights(pdnote);

for datanote0c = {"1","3",["1";"3"],"2",["1";"2"], ...
        "4","5",["4";"5"]}
    datanote0 = datanote0c{1};
    datanote = "dataset/data"+datanote0;
    if length(datanote0) == 1
        sameMo = true;
    else
        sameMo = false;
    end

for nMo = length(datanote0):2
for code = ["goldman","gold","plain0"]
for T = [125,150,175,200,250,437]
for pumpstr = ["2","2-3","3-4","2-3-4","2-3-5","2-3-4-5","2-7"]
    nTx = numel(strfind(pumpstr,"-"))+1;

for Lp = [4,8,16,32]

foldername = datanote+"/"+string(T)+"ms_"+pumpstr+"_"+string(Lp)+"_"+code;
if sum(~isfolder(foldername)), continue; end

if sameMo
    ntxt = repmat(length(dir(foldername+"/tx/*.TXT")),[nMo,1]);
else
    ntxt = zeros(nMo,1);
    for jj = 1:nMo
        dtxt = dir(foldername(jj)+"/tx/*.TXT");
        ntxt(jj) = length(dtxt);
    end
end

for T2 = T
for Lp2 = Lp

%
if false
    debug_pd = false;
    mode_pd = false;
    algoPD = "gt";
    algoCE = "gt";
    loop_emulates_txrx_all
end

%
if (cenote==cenote0 && nMo==2 && isequal(datanote0,"1") && code=="goldman" ...
        && pumpstr=="2" && T==125 && Lp==16) ... % figure 6
   || (cenote==cenote0 && nMo==1 && isequal(datanote0,"1") && code=="gold" ...
        && pumpstr=="2" && T==125 && Lp==16) ... % figure 6
   || (cenote==cenote0 && nMo==1 && isequal(datanote0,"1") && code=="plain0" ...
        && pumpstr=="2" && T==437 && Lp==2) ... % figure 6
   || (ismember(cenote,[1,5,12]+100) && nMo==1 && isequal(datanote0,"1") ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 11
   || (cenote==cenote0 && ismember(nMo,[1,2]) && (numel(datanote0)==1&&ismember(datanote0,["1","3","4","5"])) ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 12
   || (cenote==cenote0 && nMo==2 && (numel(datanote0)==2&&ismember(datanote0.',[["1";"3"],["4";"5"]].',"rows")) ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 12
   || (cenote==cenote0 && nMo==1 && (numel(datanote0)==1&&ismember(datanote0,["1","2"])) ...
        && pumpstr=="2-7" && T==125 && Lp==16) ... % figure 13
   || (cenote==cenote0 && nMo==2 && isequal(datanote0,["1";"2"]) ...
        && pumpstr=="2-7" && T==125 && Lp==16) % figure 13
    debug_pd = false;
    mode_pd = false;
    algoPD = "gt";
    algoCE = "af0";
    loop_emulates_txrx_all
end

%
if false
    debug_pd = true;
    mode_pd = true;
    algoPD = "sc";
    algoCE = "af0";
    loop_emulates_txrx_all
end

%
if (cenote==cenote0 && ismember(nMo,[1,2]) && isequal(datanote0,"1") && code=="goldman" ...
        && pumpstr=="2-3-4-5" && ismember(T,[125,150,175,200,250]) && Lp==16) ... % figure 14,15
    debug_pd = false;
    mode_pd = true;
    algoPD = "sc";
    algoCE = "af0";
    loop_emulates_txrx_all
end

%
if (cenote==cenote0 && nMo==2 && isequal(datanote0,"1") && code=="goldman" ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5","2-3-4-5-6","2-3-4-5-6-7"]) && T==125 && Lp==16) ... % figure 6,9
   || (cenote==cenote0 && nMo==1 && isequal(datanote0,"1") && code=="gold" ...
        && ismember(pumpstr,["2","2-3","2-3-4"]) && T==125 && Lp==16) ... % figure 6
   || (cenote==cenote0 && nMo==1 && isequal(datanote0,"1") && code=="plain0" ...
        && pumpstr=="2" && T==437 && Lp==2) ... % figure 6
   || (cenote==cenote0 && nMo==2 && isequal(datanote0,"1") && code=="goldman" ...
        && pumpstr=="2-3-4-5" && T==125 && ismember(Lp,[4,8,16,32])) ... % figure 8
    debug_pd = false;
    mode_pd = false;
    algoPD = "sc";
    algoCE = "af0";
    loop_emulates_txrx_all
end

end
end
end
end
end
end
end
end
end
end
