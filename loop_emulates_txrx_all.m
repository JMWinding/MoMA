warning('off');
subruns = 1e2;
maxIdx = 5;
if ~exist("debug_pd","var"), debug_pd = false; end
if ~exist("mode_pd","var"), mode_pd = false; end

%%
for subIdx = 1:maxIdx

% foldername
prefixname = "mat_temp/mat"+datanote0(1);
for j = 2:length(datanote0)
    prefixname = prefixname+"-"+datanote0(j);
end
if mode_pd
    prefixname = prefixname+"PD";
end
if debug_pd
    prefixname = prefixname+"debug";
end
prefixname = prefixname+"_11";
prefixname = prefixname+"/ce"+string(cenote);
if cenote ~= pdnote
    prefixname = prefixname+"-pd"+string(pdnote);
end
if T2<T
    prefixname = prefixname+"_"+string(T/T2)+"x";
end
if ~isfolder(prefixname), mkdir(prefixname); end

% filename
prefixname = prefixname+"/emulates_"+num2str(T)+"ms_"+pumpstr+"_"+num2str(Lp)+"_"+code;
if Lp2 ~= Lp
    prefixname = prefixname+"_Lp"+num2str(Lp2);
end
prefixname = prefixname+"_"+num2str(nMo)+"_"+algoPD+"-"+algoCE;
matname = prefixname+".mat";
clockname = prefixname+"_clock.mat";
submatname = prefixname+"("+string(subIdx)+").mat";
subclockname = prefixname+"_clock("+string(subIdx)+").mat";

matname0 = strrep(matname,"_temp","_author");

dothisloop = true;
if isfile(matname) %|| isfile(matname0)
    fprintf('Existing %s\n', matname);
    dothisloop = false;
elseif isfile(submatname)
    fprintf('Existing %s\n', submatname);
    dothisloop = false;
elseif isfile(subclockname)
    todo = false;
    while 1
        try
            load(subclockname);
            if endtime > datetime
                fprintf('Someone else is working on %s ......\n', matname);
                todo = false;
                break;
            else
                delete(subclockname);
                todo = true;
                break;
            end
        catch
            if ~isfile(subclockname)
                todo = false;
                break;
            end
            pause(rand()*5+5);
        end
    end
    if ~todo
        dothisloop = false;
    end
end

if dothisloop
    starttic = datetime;
    endtime = starttime + expdur;
    remdur = endtime - datetime;
    if remdur < predur
        fprintf('Remaining time may not finish next data.');
        return;
    end
    save(subclockname, 'endtime');
    
    %%
    disp(datetime('now'));
    fprintf('Working on %s ......\n', submatname);
    
    tic;
    
    if prod(ntxt) <= subIdx * subruns
        totalruns = max(prod(ntxt) - (subIdx-1)*subruns, 0);
    else
        totalruns = subruns;
    end
    disp(['working on nMo=' num2str(nMo) ', Lp=' num2str(Lp) ...
        ', runs=' num2str(totalruns) ', subIdx=' num2str(subIdx)]);
    
    notes_temp = nan(totalruns,nMo);
    if ~mode_pd, ber_temp = nan(totalruns,nMo,nTx); end
    if ~debug_pd, pdoff_temp = nan(totalruns,1,nTx); end
    if debug_pd, pddebug_temp = cell(totalruns,1); end
    
%     for kk = 1:totalruns
    parfor kk = 1:totalruns
        addpath("code_algo");
        warning('off');
        
        if prod(ntxt) <= maxIdx * subruns
            note_idx = kk+(subIdx-1)*subruns-1;
            note_temp = zeros(nMo,1);
            for jj = nMo:-1:1
                note_temp(jj) = mod(note_idx,ntxt(jj));
                note_idx = floor(note_idx/ntxt(jj));
            end
        else
            note_temp  = zeros(nMo,1);
            for jj = 1:nMo
                note_temp(jj) = floor(rand()*ntxt(jj));
            end
        end
        notes = arrayfun(@(a)sprintf("%02d",a), note_temp);
    
        constructIn = struct(...
            'datanote', datanote, ...
            'T', T, ...
            'pumpstr', pumpstr, ...
            'Lp', Lp, ...
            'code', code, ...
            'notes', notes, ...
            'T2', T2, ...
            'Lp2', Lp2);
        try constructIn.hPre = hPre; catch ; end
        try constructIn.hPost = hPost; catch ; end
        try constructIn.hlen = hlen; catch ; end
        rxIn = emulates_construct_rxIn(constructIn);
        
        rxIn.algoPD = algoPD;
        rxIn.algoCE = algoCE;
        rxIn.hPre = ceil(1250/T2);
        rxIn.hPost = ceil(1750/T2);
        rxIn.debug_pd = debug_pd;
        rxIn.sameMo = sameMo;
        rxIn.weights_ce = weights_ce;
        rxIn.weights_pd = weights_pd;
        if mode_pd, rxIn.mode = "pd"; end
        if debug_pd, rxIn.thrd_pd = struct("corr",1,"ratio",1); end
        if algoPD == "gt" && algoCE == "af0"
            pdoff = randi([0 2], [1 nTx]);
            for i = 1:length(rxIn.txOffset)
                rxIn.txOffset{i} = rxIn.txOffset{i} + pdoff(i);
            end
        end
        
        rxOut = decode_mmo_coherent_MMoNTxSW11loop(rxIn);
    
        [~,idx] = sort(cell2mat(rxIn.txOffset),'ascend');
        if ~mode_pd, ber_temp(kk,:,:) = cell2mat(rxOut.BER(:,idx)); end
        if isequal(algoPD,"gt")
            pdoff_temp(kk,1,:) = pdoff;
        elseif ~debug_pd
            pdoff_temp(kk,1,:) = rxOut.PDOff(idx);
        end
        if debug_pd, pddebug_temp{kk} = rxOut.debug_pd; end
        notes_temp(kk,:) = note_temp.';
    
    end
    
    save(submatname, "notes_temp");
    if ~mode_pd, save(submatname, "ber_temp", "-append"); end
    if ~debug_pd, save(submatname, "pdoff_temp", "-append"); end
    if debug_pd, save(submatname, "pddebug_temp", "-append"); end

    toc;
    
    disp(datetime('now'));
    fprintf('Finished with %s ......\n', submatname);
    
    %%
    delete(subclockname);
    predur = datetime - starttic;
end

end

%%
docombine = true;
for subIdx = 1:maxIdx
    submatname = prefixname+"("+string(subIdx)+").mat";
    if ~isfile(submatname)
        docombine = false;
        break;
    end
end

if docombine
    notes_temp = cell(maxIdx,1);
    if ~mode_pd, ber_temp = cell(maxIdx,1); end
    if ~debug_pd, pdoff_temp = cell(maxIdx,1); end
    if debug_pd, pddebug_temp = cell(maxIdx,1); end

    for subIdx = 1:maxIdx
        submatname = prefixname+"("+string(subIdx)+").mat";
        temp = load(submatname);
        notes_temp(subIdx,1) = {temp.notes_temp};
        if ~mode_pd, ber_temp(subIdx,1) = {temp.ber_temp}; end
        if ~debug_pd, pdoff_temp(subIdx,1) = {temp.pdoff_temp}; end
        if debug_pd, pddebug_temp(subIdx,1) = {temp.pddebug_temp}; end
    end

    notes_temp = cell2mat(notes_temp);
    save(matname, "notes_temp");
    if ~mode_pd
        ber_temp = cell2mat(ber_temp);
        save(matname, "ber_temp", "-append");
    end
    if ~debug_pd
        pdoff_temp = cell2mat(pdoff_temp);
        save(matname, "pdoff_temp", "-append");
    end
    if debug_pd
        pddebug_temp = cat(1,pddebug_temp{:});
        save(matname, "pddebug_temp", "-append");
    end

    for subIdx = 1:maxIdx
        submatname = prefixname+"("+string(subIdx)+").mat";
        delete(submatname);
    end
end