%% Generate the pics.
load data/data3.mat;
figure(1);

m = min(P_int,[], 'all');
M = max(P_int,[], 'all');

for idx = 1:8:256
    imagesc(reshape(P_int(idx, :), [32, 34]));colorbar;
    xlabel(sprintf('%d/256', idx));
    set(gca,'CLim',[m,M])
    saveas(gcf, sprintf('data/figs/save_%d.jpg', idx));
end

%% Make GIF.
for idx = 1:8:256
    A = imread(sprintf('data/figs/save_%d.jpg', idx));
    [im2, map2] = rgb2ind(A,256);
    if idx == 1
        imwrite(im2, map2, 'data/move.gif', 'DelayTime', 0.001, 'LoopCount', Inf);
    else
        imwrite(im2, map2, 'data/move.gif', 'WriteMode', 'append', 'DelayTime', 0.001);
    end
end

