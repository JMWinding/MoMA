linewidth_default = 2;
markersize_default = 10;
if ~exist("matversion","var"), matversion = "author"; end

%%
codeRange = ["goldman0", "ooc0", "ooc", "goldman"];
bers = zeros(1,length(codeRange));
xtickLabels = strings(1,length(codeRange));
for codeIdx = 1:length(codeRange)
    codeName = codeRange(codeIdx);

    matname = "../mat_"+matversion+"/mat1_11/ce301/emulates_125ms_2-3-4-5_16_"+codeName+"_1_gt-gt.mat";
    if isfile(matname)
        load(matname);
    else
        error("result does not exist");
    end
    bers(codeIdx) = mean(ber_temp, "all");

    switch codeName
        case "goldman0"
            xtickLabels(codeIdx) = "Gold";
        case "ooc0"
            xtickLabels(codeIdx) = "OOC";
        case "ooc"
            xtickLabels(codeIdx) = "OOCc";
        case "goldman"
            xtickLabels(codeIdx) = "MoMA";
    end
end

%%
load("../mat_Submission/mat1_noncoherent/ce301/emulates_125ms_2-3-4-5_1_ooc0_1_noncoherent.mat");
bers = [mean(ber_temp,"all"), bers];
xtickLabels = ["[64]", xtickLabels];

%%
f = figure("Position", [100 100 600 300]);
hold on; box on; grid on;
b = bar(bers, "BarWidth",0.5);
xlim([1-b.BarWidth-0.1 length(bers)+b.BarWidth+0.1]);
xticks(1:length(bers));
xticklabels(xtickLabels);
ylabel("mean BER");
set(gca, "YScale", "log");
set(f.Children(1), "FontSize", 10);

restyle(2);
% saveas(f, "fig/figure_results_codecomplement11.fig");