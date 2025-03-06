% relevant data: gFTII_%s_parcels_pm%d_%s_np%d_msh.mat, gFTMI_%s_parcels_pm%d_%s_np%d_%s.mat

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};
cic = {'con','incon', 'all'};

load('subjects.mat'); load('goodsubs.mat')

cols1 = {};
cols1{1}{1} = cat(1, repmat([228,26,28], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{1}{2} = repmat([255 184 10], 6, 1)./255;
cols1{2}{1} = cat(1, repmat([55,126,210], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{2}{2} = repmat([255 184 10], 6, 1)./255;

cols = [10 184 0; 255 184 10]./255;


%% plot
BM_params % this also produces d1 d2, 0.1 s ~ 0.5 s, sf = 160 Hz
npl = 2;
tnpl = 12/npl;
udpt = '_udpool';
fgmu = 4;

IIcic = cell(2,1); pd = 4;
for c = 1:2
    load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', conditions{c}, fgmu, udpt, npl), 'gFTII')
    data = gFTII(:, 1);
    % Trial-baseline correction
    data = cellfun(@(x) (x - nanmean(x(:, 1:d0, 1, :), 2)), data, 'UniformOutput', false);
    data = cellfun(@(x) permute(nanmean(x(:, d1:d2, :, :), 2), [1 3 4 2]), data, 'UniformOutput', false);    
    for p = 1:180
        IIcic{c}(p, :, :) = permute(data{p}(:, pd, goodsubs(:, pd, c)), [1 3 2]);
    end
end
IIcic = cellfun(@(x) nanmean(x, 3), IIcic, 'UniformOutput', false);
IIcic_diff = cellfun(@(x) x(:, 1)-x(:, 2), IIcic, 'UniformOutput', false);

%% Do stats with the rotated null distribution
% get null from Ayelet's code
Niter = 10000;
IInulli = shuffled_mat_NN(Niter);

IInull = IIcic_diff{2}(IInulli);

nperm = size(IInull, 2);
prmcr = nan(nperm, 2);
for ii = 1:nperm
    [prmcr(ii, 1), prmcr(ii, 2)] = corr(IIcic_diff{1}, IInull(:, ii), 'tail', 'right');
end
[rho, pval] = corr(IIcic_diff{1}, IIcic_diff{2}, 'tail', 'right');
cpval = sum(prmcr(:, 1)>rho)/nperm;
fprintf('Corrected pvalue: %1.8f\n', cpval)

%% for I(S;R)
target_type = 'sp'; udpt = '_udpool';
Icic = cell(2,1); pd = 4; 
for c = 1:2
    load(sprintf('gFTMI_%s_parcels_pm%d_%s_np%d_%s.mat', conditions{c}, fgmu, udpt, npl, target_type), 'gFTMI')
    data = gFTMI(:, 1);
    % Trial-baseline correction
    data = cellfun(@(x) (x - nanmean(x(:, 1:d0, 1, :), 2)), data, 'UniformOutput', false);
    data = cellfun(@(x) permute(nanmean(x(:, d1:d2, :, :), 2), [1 3 4 2]), data, 'UniformOutput', false);    
    for p = 1:180
        Icic{c}(p, :, :) = permute(data{p}(:, pd, goodsubs(:, pd, c)), [1 3 2]);
    end
end
Icic = cellfun(@(x) nanmean(x, 3), Icic, 'UniformOutput', false);
Icic_diff = cellfun(@(x) x(:, 1)-x(:, 2), Icic, 'UniformOutput', false);

%% Do stats with the rotated null distribution
% get null from Ayelet's code
Niter = 10000;
Inulli = shuffled_mat_NN(Niter);

Inull = Icic_diff{2}(Inulli);

nperm = size(Inull, 2);
prmcr = nan(nperm, 2);
for ii = 1:nperm
    [prmcr(ii, 1), prmcr(ii, 2)] = corr(Icic_diff{1}, Inull(:, ii), 'tail', 'right');
end
[rho, pval] = corr(Icic_diff{1}, Icic_diff{2}, 'tail', 'right');
cpval = sum(prmcr(:, 1)>rho)/nperm;
fprintf('Corrected pvalue: %1.8f\n', cpval)

%% PLOT I(S;R)
target_type = 'sp';
hw = 5;
vw = 4;
fgmu = 4;
pd = 4;
if pd == 4
    dtxt = 'onlys7s8';
elseif pd == 1
    dtxt = 'onlys1s2';
end

vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

sh = scatter(Icic_diff{1}, Icic_diff{2}, 36, 'k');
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;
hold all;
[b] = regress(Icic_diff{2}, [Icic_diff{1} ones(180, 1)]);
xpt = 10^-3.*[-2 -1 0.5 0 0.5 1.5 2];
plot(xpt, b(2)+b(1).*xpt, '-k', 'LineWidth', 1)
if cpval >= 10^-4
    text(10^-3*(2.5), 10^-3*(2), sprintf('Rho: %1.2f\np: %1.4f\n(corrected)', rho,cpval), 'FontSize', 6)
elseif cpval < 10^-4
    text(10^-3*(2.5), 10^-3*(2), sprintf('Rho: %1.2f\np < 10^{-4}\n(corrected)', rho), 'FontSize', 6)
end
axis('square'); 
ax = gca;
set(ax, 'tickdir', 'out', 'xtick', 10^-3.*[-2 -1 0 1 2], 'ytick', 10^-3.*[-2 -1 0 1 2], 'FontSize', 8, 'FontName', 'Helvetica')
ax.Position = [ax.Position(1)+0.05 ax.Position(2)+0.12 0.7.*ax.Position(3:4)];
yh = findobj(ax.YAxis);
yh.SecondaryLabel.HorizontalAlignment = 'left';
yh.SecondaryLabel.Position = 10^-4.*[-15 30 -1];
xh = findobj(ax.XAxis);
xh.SecondaryLabel.HorizontalAlignment = 'left';
xh.SecondaryLabel.Position = 10^-4.*[25 -10 -1];

% xlim(10^-3.*[-0.6 1.2]); ylim(10^-3.*[-0.6 1.2]);
xlim(10^-3.*[-2 2]); ylim(10^-3.*[-2.5 3]);
xlabel(sprintf('Cue I(S;R)\nConsistent - Inconsistent'), 'FontSize', 7)
ylabel(sprintf('Choice I(S;R)\nConsistent - Inconsistent'), 'FontSize', 7)
[lh,icons] = legend('Parcels'); 
lh.FontSize = 8;
lh.Position = [0.58 0.83 0.3 0.11];
lh.Box = 'off';

ic = findobj(icons, 'type', 'patch');
ic.MarkerEdgeColor = 'w';
ic.MarkerSize = 8;
ic.Vertices = [0.4 0.5];


sname = sprintf('Idiff_CC_corr_spintest_%s_%s', target_type, dtxt);
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')


%% for I(R;E)
target_type = 'est'; udpt = '_udpool';
Icic = cell(2,1); pd = 4;
for c = 1:2
    load(sprintf('gFTMI_%s_parcels_pm%d_%s_np%d_%s.mat', conditions{c}, fgmu, udpt, npl, target_type), 'gFTMI')
    data = gFTMI(:, 1);
    % Trial-baseline correction
    data = cellfun(@(x) (x - nanmean(x(:, 1:d0, 1, :), 2)), data, 'UniformOutput', false);
    data = cellfun(@(x) permute(nanmean(x(:, d1:d2, :, :), 2), [1 3 4 2]), data, 'UniformOutput', false);    
    for p = 1:180
        Icic{c}(p, :, :) = permute(data{p}(:, pd, goodsubs(:, pd, c)), [1 3 2]);
    end
end
Icic = cellfun(@(x) nanmean(x, 3), Icic, 'UniformOutput', false);
Icic_diff = cellfun(@(x) x(:, 1)-x(:, 2), Icic, 'UniformOutput', false);

%% Do stats with the rotated null distribution
% get null from Ayelet's code
Niter = 10000;
Inulli = shuffled_mat_NN(Niter);

Inull = Icic_diff{2}(Inulli);

nperm = size(Inull, 2);
prmcr = nan(nperm, 2);
for ii = 1:nperm
    [prmcr(ii, 1), prmcr(ii, 2)] = corr(Icic_diff{1}, Inull(:, ii), 'tail', 'right');
end
[rho, pval] = corr(Icic_diff{1}, Icic_diff{2}, 'tail', 'right');
cpval = sum(prmcr(:, 1)>rho)/nperm;
fprintf('Corrected pvalue: %1.8f\n', cpval)

%% PLOT I(R;E)
target_type = 'est';
hw = 5;
vw = 4;
fgmu = 4;
pd = 4;
if pd == 4
    dtxt = 'onlys7s8';
elseif pd == 1
    dtxt = 'onlys1s2';
end
vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

sh = scatter(Icic_diff{1}, Icic_diff{2}, 36, 'k');
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;
hold all;
[b] = regress(Icic_diff{2}, [Icic_diff{1} ones(180, 1)]);
xpt = 10^-3.*[-2 -1 0.5 0 0.5 1.5 2 3];
plot(xpt, b(2)+b(1).*xpt, '-k', 'LineWidth', 1)
if cpval >= 10^-4
    text(10^-3*(3), 10^-3*(2), sprintf('Rho: %1.2f\np: %1.4f\n(corrected)', rho,cpval), 'FontSize', 6)
elseif cpval < 10^-4
    text(10^-3*(3), 10^-3*(2), sprintf('Rho: %1.2f\np < 10^{-4}\n(corrected)', rho), 'FontSize', 6)
end
axis('square'); 
ax = gca;
set(ax, 'tickdir', 'out', 'xtick', 10^-3.*[-2 -1 0 1 2], 'ytick', 10^-3.*[-2 -1 0 1 2 3], 'FontSize', 8, 'FontName', 'Helvetica')
ax.Position = [ax.Position(1)+0.05 ax.Position(2)+0.12 0.7.*ax.Position(3:4)];
yh = findobj(ax.YAxis);
yh.SecondaryLabel.HorizontalAlignment = 'left';
yh.SecondaryLabel.Position = 10^-4.*[-20 35 -1];
xh = findobj(ax.XAxis);
xh.SecondaryLabel.HorizontalAlignment = 'left';
xh.SecondaryLabel.Position = 10^-4.*[25 -10 -1];

% xlim(10^-3.*[-2 2.5]); ylim(10^-3.*[-2 3.5]);
xlabel(sprintf('Cue I(R;E)\nConsistent - Inconsistent'), 'FontSize', 7)
ylabel(sprintf('Choice I(R;E)\nConsistent - Inconsistent'), 'FontSize', 7)
[lh,icons] = legend('Parcels'); 
lh.FontSize = 8;
lh.Position = [0.58 0.83 0.3 0.11];
lh.Box = 'off';

ic = findobj(icons, 'type', 'patch');
ic.MarkerEdgeColor = 'w';
ic.MarkerSize = 8;
ic.Vertices = [0.4 0.5];


sname = sprintf('Idiff_CC_corr_spintest_%s_%s', target_type, dtxt);
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')

%% PLOT II(S;R;E)
hw = 5;
vw = 4;
fgmu = 4;
pd = 4;
if pd == 4
    dtxt = 'onlys7s8';
elseif pd == 1
    dtxt = 'onlys1s2';
end
vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

sh = scatter(IIcic_diff{1}, IIcic_diff{2}, 36, 'k');
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;
hold all;
[b] = regress(IIcic_diff{2}, [IIcic_diff{1} ones(180, 1)]);
plot(10^-3.*[-0.6 0 0.5 1.2], b(2)+b(1).*10^-3.*[-0.6 0 0.5 1.2], '-k', 'LineWidth', 1)
if cpval >= 10^-4
    text(10^-3*(1.5), 10^-3*(0.8), sprintf('Rho: %1.2f\np: %1.4f\n(corrected)', rho,cpval), 'FontSize', 6)
elseif cpval < 10^-4
    text(10^-3*(1.5), 10^-3*(0.8), sprintf('Rho: %1.2f\np < 10^{-4}\n(corrected)', rho), 'FontSize', 6)
end
axis('square'); 
ax = gca;
set(ax, 'tickdir', 'out', 'xtick', 10^-3.*[-0.5 0 0.5 1], 'ytick', 10^-3.*[-0.5 0 0.5 1], 'FontSize', 8, 'FontName', 'Helvetica')
ax.Position = [ax.Position(1)+0.05 ax.Position(2)+0.12 0.7.*ax.Position(3:4)];
yh = findobj(ax.YAxis);
yh.SecondaryLabel.HorizontalAlignment = 'left';
yh.SecondaryLabel.Position = [-7.0000e-04 0.0012 -1];
xh = findobj(ax.XAxis);
xh.SecondaryLabel.HorizontalAlignment = 'right';
xh.SecondaryLabel.Position = [0.00195 -0.00040 -1];

% xlim(10^-3.*[-0.6 1.2]); ylim(10^-3.*[-0.6 1.2]);
xlabel(sprintf('Cue II(S;R;E)\nConsistent - Inconsistent'), 'FontSize', 7)
ylabel(sprintf('Choice II(S;R;E)\nConsistent - Inconsistent'), 'FontSize', 7)
[lh,icons] = legend('Parcels'); 
lh.FontSize = 8;
lh.Position = [0.58 0.83 0.3 0.11];
lh.Box = 'off';

ic = findobj(icons, 'type', 'patch');
ic.MarkerEdgeColor = 'w';
ic.MarkerSize = 8;
ic.Vertices = [0.4 0.5];


sname = sprintf('IIdiff_CC_corr_spintest_%s', dtxt);
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')

%% II(S;R;E): Do this for 22 ROIs
BM_params % this also produces d1 d2, 0.1 s ~ 0.5 s, sf = 160 Hz
npl = 2;
tnpl = 12/npl;
udpt = '_udpool';
fgmu = 4;

IIcic = cell(2,1); pd = 4;
if pd == 4
    dtxt = 'onlys7s8';
elseif pd == 1
    dtxt = 'onlys1s2';
end
for c = 1:2
    load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', conditions{c}, fgmu, udpt, npl), 'gFTII')
    data = gFTII(:, 1);
    % Trial-baseline correction
    data = cellfun(@(x) (x - nanmean(x(:, 1:d0, 1, :), 2)), data, 'UniformOutput', false);
    data = cellfun(@(x) permute(nanmean(x(:, d1:d2, :, :), 2), [1 3 4 2]), data, 'UniformOutput', false);
    for p = 1:180
        IIcic{c}(p, :, :) = permute(data{p}(:, pd, goodsubs(:, pd, c)), [1 3 2]);
    end

end
IIcic = cellfun(@(x) nanmean(x, 3), IIcic, 'UniformOutput', false);
IIcic_diff = cellfun(@(x) x(:, 1)-x(:, 2), IIcic, 'UniformOutput', false);

IIcic_diff22 = cell(2, 1);
for c = 1:2
    for a = 1:22
        areas = glasser_group(a);
        aind = ismember(gFTII(:, 2), areas);
        tmp = nanmean(IIcic_diff{c}(aind, 1));
        IIcic_diff22{c}(a, 1) = tmp;
    end
end

% PLOT
hw = 5;
vw = 4;
fgmu = 4;

vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

sh = scatter(IIcic_diff22{1}, IIcic_diff22{2}, 36, 'k');
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;
hold all;
[rho, pval] = corr(IIcic_diff22{1}, IIcic_diff22{2}, 'tail', 'right');
[b] = regress(IIcic_diff22{2}, [IIcic_diff22{1} ones(22, 1)]);
plot(10^-3.*[-0.6 0 0.5 1.2], b(2)+b(1).*10^-3.*[-0.6 0 0.5 1.2], '-k', 'LineWidth', 1)
if pval >= 10^-4
    text(10^-3*(1.5), 10^-3*(0.8), sprintf('Rho: %1.2f\np: %1.4f\n', rho,pval), 'FontSize', 6)
elseif pval < 10^-4
    text(10^-3*(1.5), 10^-3*(0.8), sprintf('Rho: %1.2f\np < 10^{-4}', rho), 'FontSize', 6)
end
axis('square'); 
ax = gca;
set(ax, 'tickdir', 'out', 'xtick', 10^-3.*[-0.5 0 0.5 1], 'ytick', 10^-3.*[-0.5 0 0.5 1], 'FontSize', 8, 'FontName', 'Helvetica', ...
    'XTickLabelRotation', 0)
ax.Position = [ax.Position(1)+0.05 ax.Position(2)+0.12 0.7.*ax.Position(3:4)];
yh = findobj(ax.YAxis);
yh.SecondaryLabel.HorizontalAlignment = 'left';
yh.SecondaryLabel.Position = [-7.0000e-04 0.0012 -1];
xh = findobj(ax.XAxis);
xh.SecondaryLabel.HorizontalAlignment = 'right';
xh.SecondaryLabel.Position = [0.00195 -0.00040 -1];

xlim(10^-3.*[-0.6 1.2]); ylim(10^-3.*[-0.6 1.2]);
xlabel(sprintf('Cue II(S;R;E)\nConsistent - Inconsistent'), 'FontSize', 7)
ylabel(sprintf('Choice II(S;R;E)\nConsistent - Inconsistent'), 'FontSize', 7)
[lh,icons] = legend('ROI'); 
lh.FontSize = 8;
lh.Position = [0.58 0.83 0.3 0.11];
lh.Box = 'off';

ic = findobj(icons, 'type', 'patch');
ic.MarkerEdgeColor = 'w';
ic.MarkerSize = 8;
ic.Vertices = [0.4 0.5];


sname = sprintf('IIdiff22_CC_corr_uncorrected_%s', dtxt);
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')


%% II(S;R): Do this for 22 ROIs
BM_params % this also produces d1 d2, 0.1 s ~ 0.5 s, sf = 160 Hz
target_type = 'sp';
npl = 2;
tnpl = 12/npl;
udpt = '_udpool';
fgmu = 4;

Icic = cell(2,1); pd = 4;
if pd == 4
    dtxt = 'onlys7s8';
elseif pd == 1
    dtxt = 'onlys1s2';
end

for c = 1:2
    load(sprintf('gFTMI_%s_parcels_pm%d_%s_np%d_%s.mat', conditions{c}, fgmu, udpt, npl, target_type), 'gFTMI')
    data = gFTMI(:, 1);
    % Trial-baseline correction
    data = cellfun(@(x) (x - nanmean(x(:, 1:d0, 1, :), 2)), data, 'UniformOutput', false);
    data = cellfun(@(x) permute(nanmean(x(:, d1:d2, :, :), 2), [1 3 4 2]), data, 'UniformOutput', false);    
    for p = 1:180
        Icic{c}(p, :, :) = permute(data{p}(:, pd, goodsubs(:, pd, c)), [1 3 2]);
    end
end
Icic = cellfun(@(x) nanmean(x, 3), Icic, 'UniformOutput', false);
Icic_diff = cellfun(@(x) x(:, 1)-x(:, 2), Icic, 'UniformOutput', false);

Icic_diff22 = cell(2, 1);
for c = 1:2
    for a = 1:22
        areas = glasser_group(a);
        aind = ismember(gFTMI(:, 2), areas);
        tmp = nanmean(Icic_diff{c}(aind, 1));
        Icic_diff22{c}(a, 1) = tmp;
    end
end

% PLOT
hw = 5;
vw = 4;
fgmu = 4;

vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

sh = scatter(Icic_diff22{1}, Icic_diff22{2}, 36, 'k');
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;
hold all;
[rho, pval] = corr(Icic_diff22{1}, Icic_diff22{2}, 'tail', 'right');
[b] = regress(Icic_diff22{2}, [Icic_diff22{1} ones(22, 1)]);
plot(10^-3.*[-0.6 0 0.5 1.2], b(2)+b(1).*10^-3.*[-0.6 0 0.5 1.2], '-k', 'LineWidth', 1)
if pval >= 10^-4
    text(10^-3*(1.5), 10^-3*(0.2), sprintf('Rho: %1.2f\np: %1.4f\n', rho,pval), 'FontSize', 6)
elseif pval < 10^-4
    text(10^-3*(1.5), 10^-3*(0.2), sprintf('Rho: %1.2f\np < 10^{-4}', rho), 'FontSize', 6)
end
axis('square'); 
ax = gca;
set(ax, 'tickdir', 'out', 'xtick', 10^-3.*[-0.5 0 0.5 1], 'ytick', 10^-3.*[-1 -0.5 0 0.5], 'FontSize', 8, 'FontName', 'Helvetica', ...
    'XTickLabelRotation', 0)
ax.Position = [ax.Position(1)+0.05 ax.Position(2)+0.12 0.7.*ax.Position(3:4)];
yh = findobj(ax.YAxis);
yh.SecondaryLabel.HorizontalAlignment = 'left';
yh.SecondaryLabel.Position = 10^-4.*[-5.5 5.5 -1];
xh = findobj(ax.XAxis);
xh.SecondaryLabel.HorizontalAlignment = 'left';
xh.SecondaryLabel.Position = 10^-4.*[12 -10 -1];

xlim(10^-3.*[-0.6 1.2]); ylim(10^-3.*[-1.2 0.6]);
xlabel(sprintf('Cue I(S;R)\nConsistent - Inconsistent'), 'FontSize', 7)
ylabel(sprintf('Choice I(S;R)\nConsistent - Inconsistent'), 'FontSize', 7)
[lh,icons] = legend('ROI'); 
lh.FontSize = 8;
lh.Position = [0.58 0.83 0.3 0.11];
lh.Box = 'off';

ic = findobj(icons, 'type', 'patch');
ic.MarkerEdgeColor = 'w';
ic.MarkerSize = 8;
ic.Vertices = [0.4 0.5];


sname = sprintf('Idiff22_CC_corr_uncorrected_%s', fig_dir, dtxt);
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')