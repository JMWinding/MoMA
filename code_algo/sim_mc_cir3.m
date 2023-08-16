function cirs = sim_mc_cir3(cirParam)
%%
try nMo = cirParam.nMo; catch nMo = 1; end
try betas = cirParam.betas; catch 
    betas = cell(nMo,1); betas(:) = {[3;0.025;1.5;1;0;0.1]};
    % W, D, d, v_parallel, v_perpendicular, tau
end
try T = cirParam.T; catch T = 1e-1; end
try Tmax = cirParam.Tmax; catch Tmax = 10; end
try mode = cirParam.mode; catch mode = 'max'; end

%%
assert(numel(betas) == nMo);

%%
t = (0:T/10:Tmax).';

cirs = cell(nMo,1);
for j = 1:nMo
    cirs{j} = zeros(size(t));
    for i = 1:size(betas{j},2)
        temp = channelmodel(betas{j}(:,i), t);
        cirs{j} = cirs{j} + temp;
    end
end

switch mode
    case 'rand'
        idx = randi([1,10]);
    case 'max'
        [~, idx] = max(cirs{1}); idx = mod(idx-1,10)+1;
    otherwise
        error('wrong sim_mc_cir mode');
end
istart = inf;
for j = 1:nMo
    cirs{j} = cirs{j}(idx:10:end);
    istart = min(istart, find(cirs{j}>1e-3,1));
end
for j = 1:nMo
    iend = find(cirs{j}>1e-3,1,'last');
    cirs{j}(1:istart-1) = [];
    cirs{j}(iend+1:end) = [];
end

end

function cir = channelmodel(b, t)
f = @(b,t) b(1) .* (4*pi*b(2).*(t-b(6))).^(-1.5) ...
    .* exp( -((b(3)-b(4)*(t-b(6))).^2+(b(5)*(t-b(6))).^2)./(4*b(2).*(t-b(6))) );
cir = f(b,t);
cir(t<=b(6)) = 0;
end