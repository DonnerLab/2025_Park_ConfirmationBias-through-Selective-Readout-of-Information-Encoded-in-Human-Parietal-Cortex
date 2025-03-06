% relevant data: gFTMI_SE_int12_cic_pm0-5_np%d%s.mat,
% gFTII_%s_parcels_pm%d_%s_np%d_msh.mat, gFTMI_SE_int12_cic_pm0-5_np%d%s.mat

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};
cic = {'con','incon', 'all'};


load('subjects.mat')


cols1 = {};
cols1{1}{1} = cat(1, repmat([228,26,28], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{1}{2} = repmat([255 184 10], 6, 1)./255;
cols1{2}{1} = cat(1, repmat([55,126,210], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{2}{2} = repmat([255 184 10], 6, 1)./255;


%% correlation II vs. I: one area
BM_params
udpt = '_udpool';
fgmu = 4;
npl = 2;
tnpl = 12/npl;

pd = 4; 

indx = 3;

load(sprintf('gFTMI_SE_int12_cic_pm0-5_np%d%s.mat', npl, udpt), 'fmise')
load('goodsubs.mat', 'goodsubs'); 


[areas, ~, areaid] = glasser_group(indx);
dSE = cell(2, 1); dII = cell(2, 1);
rho = []; pval = [];
for c = 1:2
    SEc = permute(fmise{c}(fgmu+1, 1, pd, goodsubs(:, pd, c)), [4 3 1 2]);
    SEi = permute(fmise{c}(fgmu+1, 2, pd, goodsubs(:, pd, c)), [4 3 1 2]);
    dSE{c} = SEc - SEi;
    
    load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', conditions{c}, fgmu, udpt, npl), 'gFTII')
    pin = ismember(gFTII(:, 2), areas);
    data = gFTII(pin, 1);
    IIc = {};
    IIi = {};
    for p = 1:length(areas)
        IIc{p} = nanmean(data{p}(1, d1:d2, pd, goodsubs(:, pd, c)), 2);
        IIi{p} = nanmean(data{p}(2, d1:d2, pd, goodsubs(:, pd, c)), 2);
    end
    IIc = cellfun(@(x) permute(x, [4 3 1 2]), IIc, 'UniformOutput', false);
    IIi = cellfun(@(x) permute(x, [4 3 1 2]), IIi, 'UniformOutput', false);
    IIc = nanmean(cell2mat(IIc), 2); IIi = nanmean(cell2mat(IIi), 2);
    dII{c} = IIc - IIi;
end

%% average Cue and Choice
dSE_cc = mean([dSE{1} dSE{2}], 2);
dII_cc = mean([dII{1} dII{2}], 2);

cols = [80 80 80]./255;

hw = 5;
vw = 4;

vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];
% annotation(vis, 'textbox', 'Position',[0.05 0.9 0.1 0.1],'String', sprintf('corr(\\DeltaI(S,E), \\DeltaII(S,R,E)), con - incon, dot %d', dot), 'LineStyle', 'none', 'FontSize', 10)

[x_st, y_st, x_w, y_h] = multi_axes(hw, vw, 1, 1, [1 1], [0.15 0.05 0.1 0.25], [0 0]);
ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st, y_st, x_w, y_h]);

sh = scatter(ax, dII_cc, dSE_cc, 36, cols); hold all;
sh.MarkerFaceColor = 'k';
sh.MarkerEdgeColor = 'w';
sh.MarkerFaceAlpha = 1;
sh.MarkerEdgeAlpha = 1;

[B,BINT,R,RINT,STATS] = regress(dSE_cc, [dII_cc ones(length(dII_cc), 1)]);
[rho, pval] = corr(dSE_cc, dII_cc, 'type', 'Spearman');

plot(dII_cc, B(1).*dII_cc+B(2), 'color', cols, 'LineWidth', 1)

axis('square'); 
hold on;
box off

set(ax, 'ytick', [0 0.2 0.4], 'xtick', 10^-3.*[-1 0 1 2], 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 8)

ax.XAxis.TickLength = [0.02 0];
ax.YAxis.TickLength = [0.02 0];
if indx == 17
    set(ax, 'ytick', [0 0.2 0.4], 'xtick', 10^-3.*[-1 0 1 2], 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 8)
    xh = findobj(ax.XAxis);
    xh.SecondaryLabel.HorizontalAlignment = 'right';
    xh.SecondaryLabel.Position = [0.0023 -0.019 -1];
    ylim(ax, [-0.1 0.5]); xlim(ax, [-0.001 0.0015])
    text(0.0005, 0.02, sprintf('r = %1.3f\np = %1.5f', rho, pval), 'FontSize', 7, 'Color', cols)
elseif indx == 3
    set(ax, 'ytick', [0 0.2 0.4], 'xtick', 10^-3.*[0 1 2], 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 8)
    xh = findobj(ax.XAxis);
    xh.SecondaryLabel.HorizontalAlignment = 'right';
    xh.SecondaryLabel.Position = [0.0042 -0.01 -1];
    ylim(ax, [-0.05 0.45]); xlim(ax, [-0.001 0.003])
    text(0.001, 0.05, sprintf('r = %1.3f\np = %1.5f', rho, pval), 'FontSize', 7, 'Color', cols)
end
ylabel(ax, 'I(S;E)_{con-incon}')
xlabel(ax, 'II(S;R;E)_{con-incon}')
hold on;


sname = sprintf('corr_II_ISE_cic_np%d_pm%d_dot%d_%s_avgCC', npl, fgmu, pd, areaid);
% saveas(vis, sname, 'png')
% 
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')


%% plot brain map for Cue Choice averaged

BM_params
alp = 0.05; al = 1;
fgmu = 4; npl = 2;
axw = 5;
axh = 8;

test_type = 'P'; % or 'P' (permutation)
udpt = '_udpool';
clim = [-1 1];
ctick = [-1 1];
pvalt = cell(2, 1);
tvalt = cell(2, 1);
hh = cell(2, 1);

msk = false;
lgg = logical([1 0]);
hemv = {'lateral','medial'};
tails = {'both','right'};
tl = 1; ll = 2; pd = 4;

thflag = lgg(ll); % threhold
if thflag
    sup = 'th';
    cbt = 'Rho (fdr)';
    lp = [4 4 0];
else
    sup = [];
    cbt = 'Rho';
    lp = [4 6 0];
end


load(sprintf('gFTMI_SE_int12_cic_pm0-5_np%d%s.mat', npl, udpt), 'fmise')
load('goodsubs.mat', 'goodsubs');

dSE = cell(2, 1); dII = cell(2, 1);
for c = 1:2
    SEc = permute(fmise{c}(fgmu+1, 1, pd, goodsubs(:, pd, c)), [4 3 1 2]);
    SEi = permute(fmise{c}(fgmu+1, 2, pd, goodsubs(:, pd, c)), [4 3 1 2]);
    dSE{c} = SEc - SEi;
    
    load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', conditions{c}, fgmu, udpt, npl), 'gFTII')
    
    for a = 1:22
        areas = glasser_group(a);
        pin = ismember(gFTII(:, 2), areas);
        data = gFTII(pin, 1);
        IIc = {};
        IIi = {};
        for p = 1:length(areas)
            IIc{p} = nanmean(data{p}(1, d1:d2, pd, goodsubs(:, pd, c)), 2);
            IIi{p} = nanmean(data{p}(2, d1:d2, pd, goodsubs(:, pd, c)), 2);
        end
        IIc = cellfun(@(x) permute(x, [4 3 1 2]), IIc, 'UniformOutput', false);
        IIi = cellfun(@(x) permute(x, [4 3 1 2]), IIi, 'UniformOutput', false);
        IIc = nanmean(cell2mat(IIc), 2); IIi = nanmean(cell2mat(IIi), 2);
        dII{c}(:, a) = IIc - IIi;        
        clear data
    end
end
dSE_cc = 0.5.*(dSE{1}+dSE{2});
dII_cc = 0.5.*(dII{1}+dII{2});

rho = []; pval = [];
for a = 1:22
    [rho(a,1), pval(a,1)] = corr(dII_cc(:, a), dSE_cc, 'type', 'Spearman', 'tail', 'right');
end


vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 axw axh]; % - - width height
% axes matrix
[x_st, y_st, x_w, y_h] = multi_axes(axw, axh, 1, 2, [3 4], [-0.2, 0.02, 0.005, 0.001], [0.08 0.02]);
% fdr threshold
th_rho = rho;
if thflag && ~msk
    hh = fdr(pval, alp(al));
    th_rho(~hh) = nan;
elseif thflag && msk
    tmp = ones(22,1);
    hh = fdr(pval,alp(al));
    tmp(hh) = nan;
end
% convert back to 180 parcels
th_rho180 = nan(180, 1); mask180 = nan(180,1);
for a = 1:22
    areas = glasser_group(a);
    aind = ismember(gFTII(:, 2), areas);
    th_rho180(aind, 1) = th_rho(a);
%     mask180(aind, 1) = tmp(a);
end

for hi = 1:2
    ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st, y_st(hi), x_w, y_h]); hold all;
    cfg = BM_cfg(ax, 'L', clim, vws(hi, :), tmap22, cbt, 'on', mdir);
    cfg.mask_alpha = 0.7;
    vwt = hemv{hi};
%     hm = plot_bm_opacity(cfg, th_rho180, mask180);
    hm = plot_bm(cfg, th_rho180);
    camlight(10, 10)
    hm.Label.Position = [2 hm.Label.Position(2) hm.Label.Position(3)];
    hm.FontSize = 9;
    hm.Position = [hm.Position(1)+0.05 0.07+hm.Position(2) 0.9*hm.Position(3) 0.5*hm.Position(4)];
    hm.Ticks = ctick;
end

sname = sprintf('GRP_CCavg_corrIII_BM_pm%d%s_np%d%s_%s_a%d_%s_%s_onlys7s8_right', fgmu, udpt, npl, sup, tw, 100*(1-alp(al)), tails{tl}, test_type)
print(vis, sname, '-dpng', '-r300')

%% plot the sample-wise averaged info measures TC: CC averaged
BM_params
udpt = '_udpool';
fgmu = 4;
npl = 2;
tnpl = 12/npl;


pds = 4:6;
indx = 3;


cols = [10 184 0; 255 184 10]./255;

hw = 8;
vw = 3;


psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

alp = 0.5;
mrs = 2;
lnw1 = 1;
lnw2 = 0.25;
bv1 = -0.08; bv1t = -0.0004;
bv2 = -0.12;
bt = 0.04; bt1 = 0.0002;

load(sprintf('gFTMI_SE_int12_cic_pm0-5_np%d%s.mat', npl, udpt), 'fmise')
load('goodsubs.mat', 'goodsubs');

[x_st, y_st, x_w, y_h] = multi_axes(hw, vw, 2, 1, [1 1.5], [0.15 0.015 0.02 0.25], [0.13 0.01]);
[areas, arealb, areaid] = glasser_group(indx);

csSE = cell(2, 1); csII = cell(2, 1);
for c = 1:2
    load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', conditions{c}, fgmu, udpt, npl), 'gFTII')

    pin = ismember(gFTII(:, 2), areas);
    data = gFTII(pin, 1);
    for pd = 1:tnpl
        for cs = 1:2
            SE = permute(fmise{c}(fgmu+1, cs, pd, goodsubs(:, pd, c)), [4 3 1 2]);
            sin = find(goodsubs(:, pd, c));
            tmp = nan(30, 1);
            tmp(sin) = SE;
            csSE{c}(:, pd, cs) = tmp;
            IIt = {};
            for p = 1:length(areas)
                IIt{p} = nanmean(data{p}(cs, d1:d2, pd, goodsubs(:, pd, c)), 2);
            end
            IIt = cellfun(@(x) permute(x, [4 3 1 2]), IIt, 'UniformOutput', false);
            IIt = nanmean(cell2mat(IIt), 2);
            tmp = nan(30, 1);
            tmp(sin) = IIt;
            csII{c}(:, pd, cs) = tmp;
        end
    end
end

dSE_CC = 0.5.*(csSE{1}+csSE{2});
dII_CC = 0.5.*(csII{1}+csII{2});

%% PLOT interval 2 time course: CC averaged
vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];
% plot I(S;E)
ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st(1), y_st, x_w, y_h]);
for cs = 1:2
    plot(ax, [1:3;1:3], ...
        [nanmean(dSE_CC(:, 4:6, cs), 1)+nanstd(dSE_CC(:, 4:6, cs), [], 1)./sqrt(sum(goodsubs(:, 4:6, c),1)); nanmean(dSE_CC(:, 4:6, cs), 1)-nanstd(dSE_CC(:, 4:6, cs), [], 1)./sqrt(sum(goodsubs(:, 4:6, c),1))], ...
        'color', cols(cs, :), 'LineWidth', 1);
    hold on;
    plot(ax, 1:3, nanmean(dSE_CC(:, 4:6, cs), 1), 'color', cols(cs, :), 'LineWidth', 2); hold on;
    plot(ax, 1:3, nanmean(dSE_CC(:, 4:6, cs), 1), 'marker', 'o', 'markersize', 3,'markeredgecolor', 'w', 'markerfacecolor', cols(cs, :),'LineWidth', 0.1, 'LineStyle', 'none');
    hold on;
    box off
    
    
    ax.XAxis.TickLength = [0.02 0];
    ax.YAxis.TickLength = [0.02 0];
    xlim(ax, [0.5 3.5]); ylim(ax, [-0.08 0.4])
    ylabel(ax, 'I(S;E) (bit)')
    
    pval = []; nperm = 5000;
    % from baseline
    for pd = 4:6
        tmp1e = squeeze(dSE_CC(:, pd, cs));
        mu_orig = nanmean(tmp1e);
        mu_boot = zeros(size(tmp1e, 2), nperm);
        for b=1:nperm
            % randomly assign signs
            rs = sign(rand(size(tmp1e, 1), 1)-0.5);
            mu_boot(:,b) = nanmean(rs.*tmp1e);
        end
        pval(pd-3) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
    end
    pval(pval==0) = 10^-12;
    h = fdr(pval, 0.01);
    bh = bar(find(h)-0.2*(cs-1), bv1+bt.*sign(pval(h)));
    bh.FaceColor = cols(cs, :);
    bh.EdgeColor = 'none';
    bh.BarWidth = 0.2;
    bh.BaseValue = bv1;
    hold on;
    set(ax, 'xticklabels', psmpt(4:6), 'ytick', 10^-1.*[0 2 4], 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 9, 'layer', 'top')
end

% con-incon difference
pvalcic = [];
for pd = 4:6
    tmp1e = squeeze(dSE_CC(:, pd, 1));
    tmp2e = squeeze(dSE_CC(:, pd, 2));
    tmpe = tmp1e-tmp2e;
    mu_orig = nanmean(tmpe);
    tmpce = cat(1, tmp1e, tmp2e);
    mu_boot = zeros(size(tmpe, 2), nperm);
    for b=1:nperm
        % randomly switch labels between con and incon
        rs = cat(1, ones(size(tmpe,1), 1), 2*ones(size(tmpe,1), 1));
        rs = logical(rs(randperm(length(rs)))-1);
        tmp1v = tmpce(rs, :);
        tmp2v = tmpce(~rs, :);
        tmpv = tmp1v-tmp2v;
        mu_boot(:,b) = nanmean(tmpv);
    end
    pvalcic(pd-3) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
end
pvalcic(pvalcic==0) = 10^-12;
h1 = fdr(pvalcic, 0.01);
b2 = bar(find(h1)+0.2, bv1+bt.*sign(pvalcic(h1)));
b2.FaceColor = [0 0 0];
b2.EdgeColor = 'none';
b2.BarWidth = 0.2;
b2.BaseValue = bv1;
% plot II(S;R;E)
ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st(2), y_st, x_w, y_h]);
for cs = 1:2
    plot(ax, [1:3;1:3], ...
        [nanmean(dII_CC(:, 4:6, cs), 1)+nanstd(dII_CC(:, 4:6, cs), [], 1)./sqrt(sum(goodsubs(:, 4:6, c),1)); nanmean(dII_CC(:, 4:6, cs), 1)-nanstd(dII_CC(:, 4:6, cs), [], 1)./sqrt(sum(goodsubs(:, 4:6, c),1))], ...
        'color', cols(cs, :), 'LineWidth', 1);
    hold on;
    plot(ax, 1:3, nanmean(dII_CC(:, 4:6, cs), 1), 'color', cols(cs, :), 'LineWidth', 2); hold on;
    plot(ax, 1:3, nanmean(dII_CC(:, 4:6, cs), 1), 'marker', 'o', 'markersize', 3,'markeredgecolor', 'w', 'markerfacecolor', cols(cs, :),'color', cols(cs, :), 'LineWidth', 0.1, 'LineStyle', 'none');
    hold on;
    box off
    
    
    ax.XAxis.TickLength = [0.02 0];
    ax.YAxis.TickLength = [0.02 0];
    xlim(ax, [0.5 3.5]); ylim(ax, 10^-3.*[-0.4 2])
    ax.YAxis.Exponent = -3;
    ylabel(ax, 'II(S;R;E) (bit)')
    
    pval = []; nperm = 5000;
    % from baseline
    for pd = 4:6
        tmp1e = squeeze(dII_CC(:, pd, cs));
        mu_orig = nanmean(tmp1e);
        mu_boot = zeros(size(tmp1e, 2), nperm);
        for b=1:nperm
            % randomly assign signs
            rs = sign(rand(size(tmp1e, 1), 1)-0.5);
            mu_boot(:,b) = nanmean(rs.*tmp1e);
        end
        pval(pd-3) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
    end
    pval(pval==0) = 10^-12;
    h = fdr(pval, 0.01);
    bh = bar(find(h)-0.2*(cs-1), bv1t+bt1.*sign(pval(h)));
    bh.FaceColor = cols(cs, :);
    bh.EdgeColor = 'none';
    bh.BarWidth = 0.2;
    bh.BaseValue = bv1t;
    hold on;
end
% con-incon difference
pvalcic = [];
for pd = 4:6
    tmp1e = squeeze(dII_CC(:, pd, 1));
    tmp2e = squeeze(dII_CC(:, pd, 2));
    tmpe = tmp1e-tmp2e;
    mu_orig = nanmean(tmpe);
    tmpce = cat(1, tmp1e, tmp2e);
    mu_boot = zeros(size(tmpe, 2), nperm);
    for b=1:nperm
        % randomly switch labels between con and incon
        rs = cat(1, ones(size(tmpe,1), 1), 2*ones(size(tmpe,1), 1));
        rs = logical(rs(randperm(length(rs)))-1);
        tmp1v = tmpce(rs, :);
        tmp2v = tmpce(~rs, :);
        tmpv = tmp1v-tmp2v;
        mu_boot(:,b) = nanmean(tmpv);
    end
    pvalcic(pd-3) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
end
pvalcic(pvalcic==0) = 10^-12;
h1 = fdr(pvalcic, 0.01);
b2 = bar(find(h1)+0.2, bv1t+bt1.*sign(pvalcic(h1)));
b2.FaceColor = [0 0 0];
b2.EdgeColor = 'none';
b2.BarWidth = 0.2;
b2.BaseValue = bv1t;

set(ax, 'xticklabels', psmpt(4:6), 'ytick', 10^-3.*[0 1 2], 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 9, 'layer', 'top')

sname = sprintf('CC_tc_II_ISE_cic_np%d_pm%d_%s_msh', npl, fgmu, areaid);
% saveas(vis, sname, 'png')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%     'BackgroundColor','none')


%% plot brain map MI: Cue Choice averaged
BM_params
alp = 0.05; al = 1;
fgmu = 4; npl = 2;
axw = 5;
axh = 8;
load('gminsp_un_pooled.mat', 'gminsp_udpooled')

test_type = 'P'; % or 'P' (permutation)
udpt = '_udpool';
clim = [-1 1];
ctick = [-1 1];
pvalt = cell(2, 1);
tvalt = cell(2, 1);
hh = cell(2, 1);

msk = false;
lgg = logical([1 0]);
hemv = {'lateral','medial'};
tails = {'both','right'};
tl = 1; ll = 2; pd = 4;

thflag = lgg(ll); % threhold
if thflag
    sup = 'th';
    cbt = 'Rho (fdr)';
    lp = [4 4 0];
else
    sup = [];
    cbt = 'Rho';
    lp = [4 6 0];
end
target_type = 'sp';
load(sprintf('gFTMI_SE_int12_cic_pm0-5_np%d%s.mat', npl, udpt), 'fmise')
load('goodsubs.mat', 'goodsubs');

dSE = cell(2, 1); dMI = cell(2, 1);
for c = 1:2
    SEc = permute(fmise{c}(fgmu+1, 1, pd, goodsubs(:, pd, c)), [4 3 1 2]);
    SEi = permute(fmise{c}(fgmu+1, 2, pd, goodsubs(:, pd, c)), [4 3 1 2]);
    dSE{c} = SEc - SEi;
    
    load(sprintf('gFTMI_%s_parcels_pm%d_%s_np%d_%s.mat', conditions{c}, fgmu, udpt, npl, target_type), 'gFTMI')
    
    for a = 1:22
        areas = glasser_group(a);
        pin = ismember(gFTMI(:, 2), areas);
        data = gFTMI(pin, 1);
        MIc = {};
        MIi = {};
        for p = 1:length(areas)
            MIc{p} = nanmean(data{p}(1, d1:d2, pd, goodsubs(:, pd, c)), 2);
            MIi{p} = nanmean(data{p}(2, d1:d2, pd, goodsubs(:, pd, c)), 2);
        end
        MIc = cellfun(@(x) permute(x, [4 3 1 2]), MIc, 'UniformOutput', false);
        MIi = cellfun(@(x) permute(x, [4 3 1 2]), MIi, 'UniformOutput', false);
        MIc = nanmean(cell2mat(MIc), 2); MIi = nanmean(cell2mat(MIi), 2);
        dMI{c}(:, a) = MIc - MIi;        
        clear data
    end
end
dSE_cc = 0.5.*(dSE{1}+dSE{2});
dMI_cc = 0.5.*(dMI{1}+dMI{2});

rho = []; pval = [];
for a = 1:22
    [rho(a,1), pval(a,1)] = corr(dMI_cc(:, a), dSE_cc, 'type', 'Spearman', 'tail', 'right');
end


vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 axw axh]; % - - width height
% axes matrix
[x_st, y_st, x_w, y_h] = multi_axes(axw, axh, 1, 2, [3 4], [-0.2, 0.02, 0.005, 0.001], [0.08 0.02]);
% fdr threshold
th_rho = rho;
if thflag && ~msk
    hh = fdr(pval, alp(al));
    th_rho(~hh) = nan;
elseif thflag && msk
    tmp = ones(22,1);
    hh = fdr(pval,alp(al));
    tmp(hh) = nan;
end
% convert back to 180 parcels
th_rho180 = nan(180, 1); mask180 = nan(180,1);
for a = 1:22
    areas = glasser_group(a);
    aind = ismember(gFTMI(:, 2), areas);
    th_rho180(aind, 1) = th_rho(a);
%     mask180(aind, 1) = tmp(a);
end

for hi = 1:2
    ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st, y_st(hi), x_w, y_h]); hold all;
    cfg = BM_cfg(ax, 'L', clim, vws(hi, :), tmap22, cbt, 'on', mdir);
    cfg.mask_alpha = 0.7;
    vwt = hemv{hi};
%     hm = plot_bm_opacity(cfg, th_rho180, mask180);
    hm = plot_bm(cfg, th_rho180);
    camlight(10, 10)
    hm.Label.Position = [2 hm.Label.Position(2) hm.Label.Position(3)];
    hm.FontSize = 9;
    hm.Position = [hm.Position(1)+0.05 0.07+hm.Position(2) 0.9*hm.Position(3) 0.5*hm.Position(4)];
    hm.Ticks = ctick;
end

sname = sprintf('GRP_CCavg_corrMI_BM_pm%d%s_np%d%s_%s_a%d_%s_%s_onlys7s8_right', fgmu, udpt, npl, sup, tw, 100*(1-alp(al)), tails{tl}, test_type)
% print(vis, sname, '-dpng', '-r300')