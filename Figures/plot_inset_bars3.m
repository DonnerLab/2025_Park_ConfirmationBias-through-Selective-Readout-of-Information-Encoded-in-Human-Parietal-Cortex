function [x y z] = plot_inset_bars3(x,y,z, cfg)
[a,b] = size(x);
if b>1
d1 = cfg.tw(1); d2 = cfg.tw(2);
x = mean(x(:, d1:d2),2); y = mean(y(:,d1:d2),2); %z = mean(z(:,d1:d2),2);
end
    
bh = bar([1 2 3], mean([x y z], 1));
bh.FaceColor = 'flat';
bh.CData = cfg.cols;
bh.LineStyle = 'none';
bh.BarWidth = 0.8;
hold on; plot([1 1], [nanmean(x)+std(x,[],1)/sqrt(size(x, 1)) nanmean(x)-std(x,[],1)/sqrt(size(x, 1))], 'k', 'LineWidth', 1)
plot([2 2], [nanmean(y)+std(y,[],1)/sqrt(size(y, 1)) nanmean(y)-std(y,[],1)/sqrt(size(y, 1))], 'k', 'LineWidth', 1),  hold on;
plot([3 3], [nanmean(z)+std(z,[],1)/sqrt(size(z, 1)) nanmean(z)-std(z,[],1)/sqrt(size(z, 1))], 'r', 'LineWidth', 1)
hold on;
text(0.8, 0.8*cfg.ylm(2), sprintf('p=%1.4f',cfg.pval), 'FontSize', 7); hold on;
text(2.1, 0.7*cfg.ylm(1), sprintf('*p=%1.4f',cfg.pvalc), 'FontSize', 5, 'Color', 'r');
box off
ylim(cfg.ylm); xlim([0.4 3.6]);
ax = gca;
set(ax, 'xtick', [1 2 3], 'xticklabels',{'C','I','C-I'},'tickdir', 'out', 'layer', 'top')
ax.Position = [ax.Position(1)-0.05 ax.Position(2)-0.015 ax.Position(3)*0.95 ax.Position(4)];
end