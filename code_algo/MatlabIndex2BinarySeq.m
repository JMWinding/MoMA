function seq = MatlabIndex2BinarySeq(i, nstates)
seq = zeros(nstates,1);
if nstates == 1
    if i-1 >= 1
        seq(1) = 1;
    end
else
    for k = 1:nstates
        if i-1 >= 2.^(nstates-1:-1:nstates-(k-1))*seq(1:k-1) + 2^(nstates-k)
            seq(k) = 1;
        end
    end
end
end