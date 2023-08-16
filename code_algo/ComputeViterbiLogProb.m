function rval = ComputeViterbiLogProb(y, yh, nn, np)
%%
rval = log10(normpdf(y, yh, sqrt(np*np*yh+nn*nn)));
rval = max(rval, -1e3);

end