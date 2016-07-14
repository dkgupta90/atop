% To visualize all the designs
clc;
no_cycles = 1;
no_iter = 150;
nely  = 4;
fname = 'design-';
f = figure('units','normalized','position',[0 0 1 1]);
for i = 1:1:no_cycles
    for j = 63:1:no_iter
        j
        fid = fopen(['output_design/', fname, num2str(i), '_', num2str(j), '.dat']);
        data = textscan(fid, '%f%f%f');
        rhoV = cell2mat(data(1));
        rhoX = cell2mat(data(2));
        rhoY = cell2mat(data(3));
        plot_design_mesh(rhoX, rhoY, 0.01, 1- rhoV);
        hold on;
        grid on;
        ax = gca;
        ax.XTick = 0:1/nely:3;
        ax.YTick = 0:1/nely:1;
        pause (1)
        set(gca, 'XTickLabel',[], 'YTickLabel' ,[]);
        ax.XColor = [0.5 0.5 0.5];
        ax.YColor = [0.5 0.5 0.5]
        ax.GridAlpha = 1;
        end
end