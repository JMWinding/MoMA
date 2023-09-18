function matC = get_decode_matrix(code, len1, len2, offset, shift)
%%
if nargin < 4
    offset = 0;
end
if nargin < 5
    shift = 1;
end

%%
if isempty(len2)
    len2 = ceil((len1-offset)/shift);
end

%%
matC = zeros(len1, len2);
for i = 1:len2
    for j = 1:length(code)
        idx1 = (i-1)*shift+j+offset;
        if (idx1 > 0) && (idx1 <= len1)
            matC(idx1, i) = code(j);
        end
    end
end
end