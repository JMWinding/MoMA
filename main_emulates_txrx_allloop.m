%%
% expdur = hours(24);
starttime = datetime;
predur = hours(0);
if isempty(gcp('nocreate'))
    parpool(feature('numcores'));
end
addpath("code_algo");

%%
for cenote0 = 221
    cenoteOffset = cenote0-1;

for cenote = [1,5,12]+cenoteOffset
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
for code = ["goldman","gold","plain0","goldman0","ooc","ooc0"]
for T = [125,150,175,200,250,437]
for pumpstr = ["2","2-3","3-4","2-3-4","2-3-5","2-3-4-5","2-7"]
    nTx = numel(strfind(pumpstr,"-"))+1;

for Lp = [2,4,8,16,32]

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
if (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && ismember(code,"ooc0") ...
        && pumpstr=="2-3-4-5" && T==125 && Lp==16) % figure 10
    debug_pd = false;
    mode_pd = false;
    algoPD = "gt";
    algoCE = "gt";
    loop_emulates_txrx_noncoherent
end

%
if (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && ismember(code,["goldman","goldman0","ooc","ooc0"]) ...
        && pumpstr=="2-3-4-5" && T==125 && Lp==16) % figure 10
    debug_pd = false;
    mode_pd = false;
    algoPD = "gt";
    algoCE = "gt";
    loop_emulates_txrx_all
end

%
if (ismember(cenote,cenote0) && nMo==2 && isequal(datanote0,"1") && code=="goldman" ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 6 debug
   || (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && code=="gold" ...
        && ismember(pumpstr,["2","2-3"]) && T==125 && Lp==16) ... % figure 6 debug
   || (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && code=="plain0" ...
        && pumpstr=="2" && T==437 && Lp==2) ... % figure 6 debug
   || (ismember(mod(cenote,cenoteOffset),[1,5,12]) && nMo==1 && isequal(datanote0,"1") ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 11
   || (ismember(cenote,cenote0) && ismember(nMo,[1,2]) && (numel(datanote0)==1&&ismember(datanote0,["1","3","4","5"])) && code=="goldman" ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 12
   || (ismember(cenote,cenote0) && nMo==2 && (numel(datanote0)==2&&ismember(datanote0.',[["1";"3"],["4";"5"]].',"rows")) && code=="goldman" ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 12
   || (ismember(cenote,cenote0) && nMo==1 && (numel(datanote0)==1&&ismember(datanote0,["1","2"])) && code=="goldman" ...
        && pumpstr=="2-7" && T==125 && Lp==16) ... % figure 13
   || (ismember(cenote,cenote0) && nMo==2 && isequal(datanote0,["1";"2"]) && code=="goldman" ...
        && pumpstr=="2-7" && T==125 && Lp==16) % figure 13
    debug_pd = false;
    mode_pd = false;
    algoPD = "gt";
    algoCE = "af0";
    loop_emulates_txrx_all
end

%
if (ismember(cenote,cenote0) && ismember(nMo,[1,2]) && isequal(datanote0,"1") && code=="goldman" ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % debug
   || (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && code=="gold" ...
        && ismember(pumpstr,["2","2-3"]) && T==125 && Lp==16) % debug
    debug_pd = true;
    mode_pd = true;
    algoPD = "sc";
    algoCE = "af0";
    loop_emulates_txrx_all
end

%
if (ismember(cenote,cenote0) && ismember(nMo,[1,2]) && isequal(datanote0,"1") && code=="goldman" ...
        && pumpstr=="2-3-4-5" && ismember(T,[125,150,175,200,250]) && Lp==16) % figure 14,15
    debug_pd = false;
    mode_pd = true;
    algoPD = "sc";
    algoCE = "af0";
    loop_emulates_txrx_all
end

%
if (ismember(cenote,cenote0) && nMo==2 && isequal(datanote0,"1") && code=="goldman" ...
        && ismember(pumpstr,["2","3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) ... % figure 6
   || (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && code=="gold" ...
        && ismember(pumpstr,["2","2-3"]) && T==125 && Lp==16) ... % figure 6
   || (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && code=="plain0" ...
        && pumpstr=="2" && T==437 && Lp==2) ... % figure 6
   || (ismember(cenote,cenote0) && ismember(nMo,[1,2]) && isequal(datanote0,"1") && code=="goldman" ...
        && pumpstr=="2-3-4-5" && T==125 && ismember(Lp,[4,8,16,32])) ... % figure 8
   || (ismember(cenote,cenote0) && nMo==1 && isequal(datanote0,"1") && code=="goldman" ...
        && ismember(pumpstr,["3-4","2-3-4","2-3-4-5"]) && T==125 && Lp==16) % figure 9
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
end
