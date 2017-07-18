% This function plots the deformation of the optimized design

data = dlmread('output_design/deform_data-1_21.dat');
    
X = data(:, 1);
Y = data(:, 2);
dx = data(:, 3);
dy = data(:, 4);
rho = data(:, 5);

del = 0.0;
X = reshape(X + del.*dx, 201, 201);
Y = reshape(Y + del.*dy, 201, 201);
rho = reshape(rho, 201, 201);
colormap('gray')
% scatter(X + del.*dx, Y + del.*dy, rho);

hold on;
%imagesc(1 - rho')
contourf(X, Y, 1 - rho, 'edgecolor','none');
set(gca,'visible','off')

data = dlmread('output_design/outline_data-1_21.dat');
    
Xo = data(:, 1);
Yo = data(:, 2);
dxo = data(:, 3);
dyo = data(:, 4);

rho_o = 0.5 .* ones(length(Xo), 1);
scatter(Xo + del.*dxo, Yo + del.*dyo, rho_o);
