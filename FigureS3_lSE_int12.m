%% Parameters
pc_local = 1; % 1 is cluster
if pc_local == 2
    mdir = '/Users/hpark/Servers/mountpoint1';
    fig_dir = '/Users/hpark/Research_wAlan/_manuscript/figures';
elseif pc_local == 1
    mdir = '/home/hamepark';
    fig_dir = '/home/hamepark/P03/fig';
end
subj = get_sinfo(setdiff(1:34, [12 14 19 21]), 2);

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};

dots = {'s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','s_{11}','s_{12}'};
pdots = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};


%% Get Data
fgmu = 0;

nldd = {}; n = 1;
npl = 1;

for f = 1:length(fgmu)
    fmise = cell(2, 1);
    for c = 1:2
        for s = 1:length(subj)
            sname = sprintf('%s/P03/MI/%s/ISE_%s_%s_pm%d_np%d_CC.mat', mdir, subj{s}, subj{s}, conditions{c}, fgmu(f), npl)
            try
                load(sname, 'MIv', 'MInull')
                fmise{c}(:, :, s) = MIv - nanmean(MInull, 3); % cs x pd: averaged across shuffles/subsampling
            catch
                nldd{n, 1} = sname; n = n+1;
                fmise{c}(:, :, s) = nan;
            end
        end
    end
    save(sprintf('%s/P03/MI/gFTMI_SE_int12_pm%d_np%d_CCmatched.mat', mdir, fgmu(f), npl), 'fmise')
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOT FULL TRIAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alp = 0.5;
mrs = 3;
lnw1 = 2;
lnw2 = 0.25;
bv1 = -0.08;
bv2 = -0.12;
bt = 0.04;

% cols = [10 184 0; 255 184 10; 0 0 0]./255; % con, incon
cols = [130 130 130; 0 0 0]./255; % Cue/Choice

fw = 15;
fh = 4;
vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 fw fh];

fgmu = 0;
npl = 1;
tnpl = 12/npl;
tnp = tnpl/2; % divide by number of intervals(=2)

fdrp = 0.01;
nperm = 5000;

load(sprintf('%s/P03/MI/gFTMI_SE_int12_pm%d_np%d_CCmatched.mat', mdir, fgmu, npl), 'fmise')

cfg = []; cfg.nperm = 5000; cfg.ttest = false; cfg.perm = 'unpaired'; cfg.tail = 'both';
cs = 3; % all samples 
for c = 1:2
    [x_st, y_st, x_w, y_h] = multi_axes(fw, fh, 2, 1, [1 2.2], [0.1 0.02 0.002 0.15],[0.05 0.02]);
    ax = axes('Units', 'centimeters', 'Position', [x_st(c), y_st, x_w, y_h]);
    for int = 1:2
        plot(1+tnp*(int-1):tnp*int, squeeze(fmise{c}(cs, 1+tnp*(int-1):tnp*int, :)), 'LineStyle', '-','color', cols(c, :), 'LineWidth', 0.1);
        hold on;
        
        plot(1+tnp*(int-1):tnp*int, squeeze(nanmean(fmise{c}(cs, 1+tnp*(int-1):tnp*int, :), 3)), 'marker', 'o', 'markersize', mrs, 'color', cols(c, :), ...
            'markerfacecolor', cols(c, :), 'markeredgecolor', 'w', 'LineWidth', lnw1)
        hold on;
    end
    % do stats
    pval1 = [];
    
    % from baseline
    cfg.perm = 'zero';
    for pd = 1:tnpl
        x = permute(fmise{c}(cs, pd, :), [3 1 2]); y = zeros(length(subj), 1);
        pval1(pd,cs) = perm_test_behav(x,y,cfg);
    end
    
    pval1(pval1==0) = 10^-30;
    h = fdr(pval1(:,cs), fdrp);
    % bh = bar(find(h(12*(cs-1)+(1:12)))-0.2*(cs-1), bv1+bt.*sign(pval1(h(12*(cs-1)+(1:12)))));
    bh = bar(find(h), bv1+bt.*sign(pval1(h,cs)));
    bh.FaceColor = cols(c, :);
    bh.EdgeColor = 'none';
    bh.BarWidth = 0.2;
    bh.BaseValue = bv1;
    hold on;
    hold on;
    plot([0.7 tnpl+0.3], [0 0], ':k', 'LineWidth', 0.5)
    hold on;
    axo = gca;
    if npl == 1
        set(axo, 'xtick', [1 6 7 12], 'xticklabels', dots([1, 6, 7, 12]), 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    else
        set(axo, 'xtick', 1:tnpl, 'xticklabels', pdots(1:tnpl), 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    end
    set(axo, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    %     ax.XAxis.LineWidth = 1.5;
    axo.XAxis.TickLength = [0.02 0];
    %     ay.YAxis.LineWidth = 1.5;
    axo.YAxis.TickLength = [0.02 0];
    xlim(axo, [0.5 tnpl+0.5]); ylim(axo, [bv2 0.7])
    if c == 1
        ylabel(axo, 'I(S;E) (bit)','FontSize', 8)
    end
    title(condt{c}, 'FontSize', 8)
    box off    
end
sname = sprintf('%s/FigureS3_ISE_int12_pm%d_np%d_CC_all', fig_dir, fgmu, npl);
saveas(vis, sname, 'png')
exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
    'BackgroundColor','none')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLE MEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
