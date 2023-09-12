linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
ceName = "ce301";

%%
dataRange = ["1", "2", "1-2"];
bers = zeros(2,2);
for codeIdx = 1:length(dataRange)
    dataName = dataRange(codeIdx);
    nMo = numel(strfind(dataName,"-"))+1;

    matName = "../mat_"+matversion+"/mat"+dataName+"_11/"+ceName+"/emulates_125ms_2-7_16_goldman_"+string(nMo)+"_gt-af0.mat";

    if isfile(matName)
        disp(matName);
        load(matName);
    else
        error("result does not exist");
    end

    switch dataName
        case "1"
            bers(1,1) = mean(ber_temp, "all");
        case "2"
            bers(2,1) = mean(ber_temp, "all");
        case "1-2"
            bers(1,2) = mean(ber_temp(:,1,:), "all");
            bers(2,2) = mean(ber_temp(:,2,:), "all");
    end
end

%%
f = figure("Position", [100 100 600 300]);
hold on; box on; grid on;
b = bar(bers);
xlabel("Molecule index");
xticks(1:2);
xticklabels(["A", "B"]);
ylabel("mean BER");
legend(["w/o L3", "w/ L3"], "Location", "northwest");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_codetuple11.fig");