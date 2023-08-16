function rChip = get_gold_code(nTx, degree)
nSeq = power(2,degree)-1;
nMax = get_max_num(degree);
if nTx == 0
    nTx = nMax;
end
if nTx > nMax
    error('Please select higher order polynomials');
end

%%
switch degree
    % primpoly(8,'all');
    case 3
        polynomial1 = "x^3+x^2+1";
        polynomial2 = "x^3+x+1";
    case 4
        polynomial1 = "x^4+x^3+1";
        polynomial2 = "x^4+x+1";
    case 5
        polynomial1 = "x^5+x^2+1";
        polynomial2 = "x^5+x^4+x^3+x^2+1";
    case 6
        polynomial1 = "x^6+x+1";
        polynomial2 = "x^6+x^5+x^2+x+1";
    case 7
        polynomial1 = "x^7+x^3+1";
        polynomial2 = "x^7+x^3+x^2+x+1";
    case 9
        polynomial1 = "x^9+x^4+1";
        polynomial2 = "x^9+x^6+x^4+x^3+1";
    case 10
        polynomial1 = "x^10+x^3+1";
        polynomial2 = "x^10+x^8+x^3+x^2+1";
    case 11
        polynomial1 = "x^11+x^2+1";
        polynomial2 = "x^11+x^8+x^5+x^2+1";
    otherwise
        error('not supported degree');
end

nChip = 2^degree-1;
goldseq = comm.GoldSequence("FirstPolynomial", polynomial1, ...
    "SecondPolynomial", polynomial2, ...
    "FirstInitialConditions", [zeros(1,degree-1), 1], ...
    "SecondInitialConditions", [zeros(1,degree-1), 1], ...
    "Index", 0, ...
    "SamplesPerFrame", nChip);

%%
rChip = nan(nChip, nTx);

ind = randperm(nSeq)-1;
j = 1;
for i = 1:nSeq
    release(goldseq);
    goldseq.Index = ind(i);
    seq = goldseq();
    seq = 2 * seq - 1;

    rChip(:,j) = seq;
    j = j + 1;
    if j > nTx
        break;
    end
end

rChip2 = [rChip; -rChip];
rChip2(1:2:end,:) = rChip;
rChip2(2:2:end,:) = -rChip;
rChip = rChip2;

assert(sum(isnan(rChip), 'all') == 0);

end

function nMax = get_max_num(degree)
    nMax = 2^degree-1;
end