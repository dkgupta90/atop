% This program uses k-means clustering approach to generate a design
% distribution of N points

clc;
clear all;
close all;

n_max = 100;

opts = statset('MaxIter', 1000);
fid = fopen('designField.dat', 'w');
for n = 1:1:n_max
    close all;
    sample_count = 1000 * n;  %samples for k-means clustering
    dim = 2;    %spatial dimensions
    xmin = -1; xmax = 1; ymin = -1; ymax = 1;   %bounds for x- and y-dimensions

    samples = rand(sample_count, dim);
    samples = samples*2 - 1;

    [idx, C] = kmeans(samples, n, 'Options', opts);
%     plot(C(:, 1), C(:, 2), '*');
%     hold on;
%     axis([-1, 1, -1, 1]);
%     hold off;
    C
    C = C';
    C = C(:);
    C = [n C'];
    fprintf(fid, '%f\t', C);
    fprintf(fid, '\n');
    %pause(1)
end

fclose(fid);