function dummy = plot_design_mesh(rhoX, rhoY, rhoR, rhoV)

h = scatter(rhoX, rhoY, rhoR, rhoV, 'filled'); % Create a scatter plot and return a handle to the 'hggroup' object

%Obtain the axes size (in axpos) in Points

currentunits = get(gca,'Units');

set(gca, 'Units', 'Points');

axpos = get(gca,'Position');

set(gca, 'Units', currentunits);

markerRadius = rhoR./diff(ylim).*axpos(4); % Calculate Marker width in points

set(h, 'SizeData', 4 .* markerRadius .^2);

end



%================================================
