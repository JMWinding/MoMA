function restyle(widthratio)
set(gcf, "Units","inches", "Position",[1 1 2.32*widthratio 2.32]);

hLegend = findobj(gcf, 'Type', 'Legend');
for hLegendObj = hLegend
    hLegendObj.FontSize = 8;
    hLegendObj.FontName = "Linux Libertine";
end

hAxes = findobj(gcf, 'Type', 'Axes');
for hAxesObj = hAxes.'
    hAxesObj.FontSize = 11;
    hAxesObj.FontName = "Linux Libertine";    
    for k = 1:length(hAxesObj.Children)
        if ~isprop(hAxesObj.Children(k),"DisplayName")
            continue;
        end
        switch hAxesObj.Children(k).DisplayName
            case {"salt-1", "MDMA", "L0+L1+L2"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#00426E"); % darkest blue
            case {"salt-2", "MDMA+CDMA", "L0+L1", ...
                    "w/o L3", "all detected", "1 molecule"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#0071BD"); % dark blue
            case {"salt-mix", "MoMA", "L0+L2", ...
                    "w/ L3", "missing", "2 molecules"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#75C8FF"); % light blue
            case {}
                SetBarOrPlotColor(hAxesObj.Children(k),"#D1EFFF"); % lightest blue
            case {"soda-1"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#3D550C"); % darkest green
            case {"soda-2"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#59981A"); % dark green
            case {"soda-mix"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#81B622"); % light green
            case {"mix"}
                SetBarOrPlotColor(hAxesObj.Children(k),"#FEDE00"); % yellow
            case {"true pos", "false pos"}
            otherwise
                SetBarOrPlotColor(hAxesObj.Children(k),"#0071BD"); % dark blue
        end
    end
end

%%
function SetBarOrPlotColor(object, color)
switch get(object,"Type")
    case "line"
        object.Color = color;
        object.Marker = "d";
        object.MarkerSize = 6;
    case "bar"
        object.FaceColor = color;
end
end
end