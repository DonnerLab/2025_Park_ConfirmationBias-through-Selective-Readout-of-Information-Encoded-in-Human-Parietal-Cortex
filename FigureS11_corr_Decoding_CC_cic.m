% relevant data: %s_pm%d_btS_btE_btSE_abs.mat, gFTII_%s_parcels_pm%d_%s_np%d_msh.mat

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};
cic = {'con','incon', 'all'};

subj = get_sinfo(setdiff(1:34, [12 14 19 21]), 2);

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

DCcic = cell(2,1); pd = 4;
for c = 1:2
    load(sprintf('%s_pm%d_btS_btE_btSE_abs.mat', condt{c}, fgmu), 'gCORR','ROI')
    for p=1:size(gCORR,1)
        for cs = 1:2
            DCcic{c}(p, cs, :) = permute(nanmean(gCORR{p,1}{cs, pd}(3, d1:d2, :), 2), [3 1 2]);
        end
    end
end
DCcic = cellfun(@(x) nanmean(x, 3), DCcic, 'UniformOutput', false);
DCcic_diff = cellfun(@(x) x(:, 1)-x(:, 2), DCcic, 'UniformOutput', false);

%% Do stats with the rotated null distribution
% get null from Ayelet's code
Niter = 10000;
DCnulli = shuffled_mat_NN(Niter, mdir);

DCnull = DCcic_diff{2}(DCnulli);

nperm = size(DCnull, 2);
prmcr = nan(nperm, 2);
for ii = 1:nperm
    [prmcr(ii, 1), prmcr(ii, 2)] = corr(DCcic_diff{1}, DCnull(:, ii), 'tail', 'right');
end
[rho, pval] = corr(DCcic_diff{1}, DCcic_diff{2}, 'tail', 'right');
cpval = sum(prmcr(:, 1)>rho)/nperm;
fprintf('Corrected pvalue: %1.8f\n', cpval)


%% PLOT
setxt = {'$$corr(S,\hat{S})$$','$$corr(E,\hat{E})$$','$$corr(\hat{S},E)$$'};

hw = 5;
vw = 4;
fgmu = 4;

vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

sh = scatter(DCcic_diff{1}, DCcic_diff{2}, 36, 'k');
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;
hold all;
[b] = regress(DCcic_diff{2}, [DCcic_diff{1} ones(180, 1)]);
plot(10^-2.*[-1 0 1 2 3], b(2)+b(1).*10^-2.*[-1 0 1 2 3], '-k', 'LineWidth', 1)
% text(10^-2*(2), 10^-2*(0), sprintf('Rho: %1.2f\np: %1.8f\n(corrected)', rho,cpval))
text(10^-2*(2.9), 10^-2*(0.3), sprintf('Rho: %1.2f\np < 10^{-5} \n(corrected)', rho), 'FontSize', 7)
axis('square'); 
ax = gca;
set(ax, 'tickdir', 'out', 'xtick', 10^-2.*[-1 0 1 2 3], 'ytick', 10^-2.*[-1 0 1 2], 'FontSize', 8, 'FontName', 'Helvetica')
ax.Position = [ax.Position(1)+0.13 ax.Position(2)+0.16 0.7.*ax.Position(3:4)];
ax.YAxis.Exponent = -2;
ax.XAxis.Exponent = -2;
yh = findobj(ax.YAxis);
yh.SecondaryLabel.HorizontalAlignment = 'left';
% yh.SecondaryLabel.Position = [-7.0000e-04 0.0012 -1];
xh = findobj(ax.XAxis);
xh.SecondaryLabel.HorizontalAlignment = 'right';
xh.SecondaryLabel.Position = [0.062 -0.006 -1];

xlim(10^-2.*[-2 4]); ylim(10^-2.*[-1 3]);
xlabel(sprintf('Cue %s \nConsistent - Inconsistent', setxt{3}), 'FontSize', 8, 'Interpreter', 'Latex', 'FontName', 'Helvetica')
ylabel(sprintf('Choice %s \nConsistent - Inconsistent', setxt{3}), 'FontSize', 8, 'Interpreter', 'Latex', 'FontName', 'Helvetica')
[lh,icons] = legend('Parcels'); 
lh.FontSize = 8;
lh.Position = [0.58 0.83 0.3 0.11];
lh.Box = 'off';

ic = findobj(icons, 'type', 'patch');
ic.MarkerEdgeColor = 'w';
ic.MarkerSize = 8;
ic.Vertices = [0.4 0.5];

sname = sprintf('DCdiff_CC_corr_spintest');
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')


%% do correlation between II cic and DC cic
load('goodsubs.mat'); 
BM_params % this also produces d1 d2, 0.1 s ~ 0.5 s, sf = 160 Hz
npl = 2;
tnpl = 12/npl;
udpt = '_udpool';
fgmu = 4;

DCcic = cell(2,1); pd = 4;
for c = 1:2
    load(sprintf('%s_pm%d_btS_btE_btSE_abs.mat', condt{c}, fgmu), 'gCORR','ROI')
    for p=1:size(gCORR,1)
        for cs = 1:2
            DCcic{c}(p, cs, :) = permute(nanmean(gCORR{p,1}{cs, pd}(3, d1:d2, :), 2), [3 1 2]);
        end
    end
end
DCcic = cellfun(@(x) nanmean(x, 3), DCcic, 'UniformOutput', false);
DCcic_diff = cellfun(@(x) x(:, 1)-x(:, 2), DCcic, 'UniformOutput', false);

IIcic = cell(2,1); pd = 4;
for c = 1:2
    load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', mdir, conditions{c}, fgmu, udpt, npl), 'gFTII')
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



%% plot 
hw = 5;
vw = 4;
fgmu = 4;

for c = 1:2
    
    null = IIcic_diff{c}(nulli);
    
    nperm = size(null, 2);
    prmcr = nan(nperm, 2);
    for ii = 1:nperm
        [prmcr(ii, 1), prmcr(ii, 2)] = corr(DCcic_diff{c}, null(:, ii), 'tail', 'right');
    end
    [rho, pval] = corr(DCcic_diff{c}, IIcic_diff{c}, 'tail', 'right');
    cpval = sum(prmcr(:, 1)>rho)/nperm;
    fprintf('Corrected pvalue: %1.8f\n', cpval)
    
    
    vis = figure;
    vis.Units = 'centimeter';
    vis.Position = [1 1 hw vw];
    
    
    sh = scatter(DCcic_diff{c}, IIcic_diff{2}, 36, 'k');
    sh.MarkerFaceColor = 'k';
    sh.MarkerEdgeColor = 'w';
    sh.MarkerFaceAlpha = 1;
    sh.MarkerEdgeAlpha = 1;
    hold all;
    [b] = regress(IIcic_diff{c}, [DCcic_diff{c} ones(180, 1)]);
    plot(10^-2.*[-1 0 1 2 3], b(2)+b(1).*10^-2.*[-1 0 1 2 3], '-k', 'LineWidth', 1)
    % text(10^-2*(2), 10^-2*(0), sprintf('Rho: %1.2f\np: %1.8f\n(corrected)', rho,cpval))
    if cpval < 1/2*nperm
        text(10^-2*(3.2), 10^-3*(0.16), sprintf('Rho: %1.2f\np < 10^{-5} \n(corrected)', rho), 'FontSize', 7)
    else
        text(10^-2*(3.2), 10^-3*(0.16), sprintf('Rho: %1.2f\np: %1.5f \n(corrected)', rho, cpval), 'FontSize', 7)
    end
    axis('square'); 
    ax = gca;
    set(ax, 'tickdir', 'out', 'xtick', 10^-2.*[-1 0 1 2 3], 'ytick', 10^-3.*[0 1], 'FontSize', 8, 'FontName', 'Helvetica')
    ax.Position = [ax.Position(1)+0.1 ax.Position(2)+0.18 0.7.*ax.Position(3:4)];
    ax.YAxis.Exponent = -3;
    ax.XAxis.Exponent = -2;
    yh = findobj(ax.YAxis);
    yh.SecondaryLabel.HorizontalAlignment = 'left';
    yh.SecondaryLabel.Position = [-0.01 0.0012 -1];
    xh = findobj(ax.XAxis);
    xh.SecondaryLabel.HorizontalAlignment = 'right';
    xh.SecondaryLabel.Position = [0.05 -0.0005 -1];
    
    xlim(10^-2.*[-1 3.5]); ylim(10^-3.*[-0.7 1.3]);
    xlabel(sprintf('%s %s \nConsistent - Inconsistent', condt{c}, setxt{3}), 'FontSize', 8, 'Interpreter', 'Latex', 'FontName', 'Helvetica')
    ylabel(sprintf('%s II(S;R;E)\nConsistent - Inconsistent', condt{c}), 'FontSize', 7)
    [lh,icons] = legend('Parcels');
    lh.FontSize = 8;
    lh.Position = [0.58 0.83 0.3 0.11];
    lh.Box = 'off';
    
    ic = findobj(icons, 'type', 'patch');
    ic.MarkerEdgeColor = 'w';
    ic.MarkerSize = 8;
    ic.Vertices = [0.4 0.5];
    
    sname = sprintf('%s_DCdiff_IIdiff_corr_spintest', condt{c});
    % saveas(vis, sname, 'png')
%     exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%         'BackgroundColor','none')
end