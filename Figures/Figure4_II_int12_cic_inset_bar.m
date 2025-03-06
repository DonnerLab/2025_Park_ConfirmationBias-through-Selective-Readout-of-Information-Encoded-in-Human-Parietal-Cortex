% relevant data: gFTII_%s_parcels_pm%d_%s_np%d_msh.mat, gFTII_%s_parcels_pm%d_%s_np%d_msh_PC1.mat

conditions = {'all', 'attn', 'choice'};
cic = {'con','incon'};
condt = {'Cue','Choice'};


load('subjects.mat'); load('goodsubs.mat')

cols = {};
cols{1} = [10 184 0]./255;
cols{2} = [255 184 10]./255;


%% plot

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

fgmu = 4;

pl = true;
if pl
    udpt = '_udpool';
elseif ~pl
    udpt = []; 
end

% % if PC1
% udpt = 'udpool';

np = 2;
npl = 2;

% for permutation test
addpath('fieldtrip-20171001')
ft_defaults

BM_params

ARG.pcrit = 0.05;  
ARG.Nperm = 5000;
ARG.T_Range =[xtime(1) xtime(end)]; 
ARG.tw = [d1 d2];

cfg = [];
cfg.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg.minsize = 1; % minimal cluster size
cfg.pval = 0.05; % threshold to select signifciant clusters
cfg.df = length(subj)-1; %degrees of freedom. Only needed for effect size.
cfg.smooth = true;
cfg.ttest = false;
cfg.perm = 'unpaired';

ylm = 10^-3.*[-0.1 1.5];

% if PC1
% ylm = 10^-3.*[-0.015 2.5];

indx = 3;
[areas,~,areaid] = glasser_group(indx);

pds = [3 5]; % this is the pooled sample index: 1: s1s2, 2: s3s4, 3: s5s6, 4: s7s8, 5: s9s10, 6: s11s12

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 5 5];

pvalcic = [];
for c = 2:3
    minsp = gminsp_udpooled{fgmu+2, c-1};
    load(sprintf('%s/P03/II/gFTII_%s_parcels_pm%d_%s_np%d_msh.mat', mdir, conditions{c}, fgmu, udpt, npl), 'gFTII')
    % if PC1, load gFTIIpc1 and set
%     gFTII = gFTIIpc1;
    
    inda = ismember(gFTII(:, 2), areas);
    data = gFTII(inda, 1);
    npd = size(gFTII{1,1}, 3);
    
    for pd = 1:length(pds)
        subplot(2, length(pds), pd+length(pds)*(c-2)); hold all;
        dat = cell(2, 1);

        xline(0, ':', 'HandleVisibility', 'off'); yline(0, '-', 'HandleVisibility', 'off')
        for cs = 1:2
            tmp = cellfun(@(x) permute((x(cs, :, pds(pd), goodsubs(:, pd, c))), [1 2 4 3]), data, 'UniformOutput', false);
            tmp = cell2mat(tmp);
            tmp = permute(nanmean(tmp, 1), [3 2 1]); dat{cs} = tmp;
        end
        
        %%%%%%%%%%% test 0.1-0.5 window %%%%%%%%%%%
        [~, pvalcic(c-1,pd)] = perm_test(dat{1},dat{2},ARG,cfg);
        text(1.3, 1.3*10^-3, sprintf('*p=%1.4f',pvalcic(c-1,pd)));
        cdat = cellfun(@(x) nanmean(x(:, ARG.tw(1):ARG.tw(2)), 2), dat, 'UniformOutput', false);
        bh = bar([1 2], mean([cdat{1} cdat{2}], 1));
        bh.FaceColor = 'flat';
        bh.CData = [cols{1};cols{2}];
        bh.LineStyle = 'none';
        bh.BarWidth = 0.8;
        hold on; plot([1 1], [nanmean(cdat{1})+std(cdat{1},[],1)/sqrt(size(cdat{1}, 1)) nanmean(cdat{1})-std(cdat{1},[],1)/sqrt(size(cdat{1}, 1))], 'k', 'LineWidth', 1)
        plot([2 2], [nanmean(cdat{2})+std(cdat{2},[],1)/sqrt(size(cdat{2}, 1)) nanmean(cdat{2})-std(cdat{2},[],1)/sqrt(size(cdat{2}, 1))], 'k', 'LineWidth', 1)

        box off
        ylim(ylm); xlim([0.4 2.6]); if pd == 1, ylabel(condt{c-1}), end
        ax = gca;
        set(ax, 'xtick', [1 2], 'xticklabels',{'C','I'},'ytick', 10^-3.*[0 1 2], 'tickdir', 'out', 'layer', 'top')
        ax.Position = [ax.Position(1)+0.01 ax.Position(2)-0.01 ax.Position(3)*0.9 ax.Position(4)*0.8];
    end
end
set(findall(vis, '-property','FontName'),'FontName','Helvetica')
set(findall(vis, '-property','FontSize'),'FontSize',7)
sname = sprintf('%s/GRP_II_%s_pm%d%s_np%d_insetBar_msh_s5s6', fig_dir, areaid, fgmu, udpt, npl);
%     saveas(vis, sname, 'png')
exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
           'BackgroundColor','none')

