[figures, pathname, ] = uigetfile('*.fig', 'MultiSelect','on');
if ~iscell(figures)
    figures = {figures};
end

for x = 1:length(figures)
    Multi_Figs = [pathname, filesep, figures{x}];
    fig = openfig(figures{x});
    dotLoc = find(figures{x} == '.');
    name = figures{x}(1:dotLoc(1)-1);
    exportgraphics(fig, name + ".pdf");
    close(fig);
end
