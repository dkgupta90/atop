% To visualize all the designs
clc;
no_cycles = 1;
no_iter = 72;
nely  = 20;
fname = 'design-';
f = figure('units','normalized','position',[0 0 1 1]);
for i = 1:1:no_cycles
    for j = 1:1:no_iter
        j
        fid = fopen(['output_design/', fname, num2str(i), '_', num2str(j), '.dat']);
        data = textscan(fid, '%f%f%f%f');
        rhoV = cell2mat(data(1));
        rhoR = cell2mat(data(2)) ./5;
        rhoX = cell2mat(data(3));
        rhoY = cell2mat(data(4));
        plot_design_mesh(rhoX, rhoY, (rhoR./nely), rhoV);
        pause (1)
    end
end