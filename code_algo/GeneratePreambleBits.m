function mseq = GeneratePreambleBits(len)
%%
p = nextpow2(len/2+1);
pr = primpoly(p, "nodisplay");
fbconnection = de2bi(pr);

n = length(fbconnection);

rbits = [1,1,0,1,0,1,0,0,1,0];
register = flip(rbits(1:n));
% register = randi([0,1], [1,n]);

newregister = register;

mseq = zeros(len,1);
mseq(1) = register(n);
for i=2:len/2
    newregister(1)=mod(sum(fbconnection.*register),2);
    for j=2:n
        newregister(j)=register(j-1);
    end
    register=newregister;
    mseq(i)=register(n);
end

mseq(len/2+1:end) = 1 - mseq(1:len/2);
mseq = 2*mseq - 1;

end