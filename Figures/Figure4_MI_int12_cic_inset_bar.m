% relevant data: gFTMI_%s_parcels_pm%d_%s_np%d.mat, gFTMI_%s_parcels_pm%d_%s_np%d_PC1.mat


conditions = {'all', 'attn', 'choice'};
cic = {'con','incon'};
condt = {'Cue','Choice'};

load('subjects.mat');

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
np = 2;

% for permutation test
addpath('fieldtrip-20171001')
ft_defaults

ARG.pcrit = 0.05;  
ARG.Nperm = 2000;
% ARG.T_Range =[-0.1 1]; % time range to display

cfg = [];
cfg.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg.minsize = 1; % minimal cluster size
cfg.pval = 0.05; % threshold to select signifciant clusters
cfg.df = length(subj)-1; %degrees of freedom. Only needed for effect size.


BM_params

ylm = 10^-3.*[-0.25 5];

indx = 3;

load('gminsp_un_pooled.mat')

pds = [1 4]; % this is the pooled sample index: 1: s1s2, 2: s3s4, 3: s5s6, 4: s7s8, 5: s9s10, 6: s11s12

for tar = 1:2
    if tar == 1
        target_type = 'sp';
        infot = 'I(S;R)';
    elseif tar == 2
        target_type = 'est';
        infot = 'I(R;E)';
    end
    vis = figure;
    vis.Units = 'centimeters';
    vis.Position = [2 2 5 5];

    pvalcic = [];
    for c = 2:3
        minsp = gminsp_udpooled{fgmu+2, c-1};
        load(sprintf('gFTMI_%s_parcels_pm%d_%s_np%d_%s_PC1.mat', conditions{c}, fgmu, udpt, npl, target_type), 'gFTMI')
        [areas,~,areaid] = glasser_group(indx);
        
        inda = ismember(gFTMI(:, 2), areas);
        data = gFTMI(inda, 1);
        npd = size(gFTMI{1,1}, 3);
        
        for pd = 1:length(pds)
            subplot(2, length(pds), pd+length(pds)*(c-2)); hold all;
            dat = cell(2, 1);
            minsp = min(gminsp_udpooled{fgmu+2, c-1}(:, :, pds(pd)), [], 2);
            goodsub = minsp >= 27;
            xline(0, ':', 'HandleVisibility', 'off'); yline(0, '-', 'HandleVisibility', 'off')
            for cs = 1:2
                tmp = cellfun(@(x) permute((x(cs, :,pds(pd), goodsub)), [1 2 4 3]), data, 'UniformOutput', false);
                tmp = cell2mat(tmp);
                tmp = permute(nanmean(tmp, 1), [3 2 1]); dat{cs} = tmp;
            end
            
            %%%%%%%%%%% test 0.1-0.5 window %%%%%%%%%%%
            nperm = 5000;
            cdat = cellfun(@(x) nanmean(x(:, d1:d2), 2), dat, 'UniformOutput', false);
            tmp2 = cdat{1}-cdat{2};
            mu_orig = mean(tmp2);
            tmpc = cat(1, cdat{1}, cdat{2});
            mu_boot = zeros(size(tmp2,2),nperm);
            for b=1:nperm
                % randomly switch labels between baseline and data
                rs = cat(1, ones(size(tmp2,1), 1), 2*ones(size(tmp2,1), 1));
                rs = logical(rs(randperm(length(rs)))-1);
                tmp1v = tmpc(rs, :);
                tmp2v = tmpc(~rs, :);
                tmpv = tmp1v-tmp2v;
                mu_boot(:,b) = mean(tmpv);
            end
            pvalcic(c-1, pd) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
            text(1.3, 3*10^-3, sprintf('*p=%1.4f',pvalcic(c-1, pd)));
            
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
            set(ax, 'xtick', [1 2], 'xticklabels',{'C','I'},'ytick', 10^-3.*[0 2 4], 'tickdir', 'out', 'layer', 'top')
            ax.Position = [ax.Position(1)-0.05 ax.Position(2)-0.015 ax.Position(3)*1.02 ax.Position(4)*0.8];
        end
    end
    sname = sprintf('GRP_MI_%s%d_pm%d%s_np%d_inset_se_%s_PC1', areas{1}, sum(indx), fgmu, udpt, np, target_type);
%     saveas(vis, sname, 'png')
    hgexport(vis, sname)
end

%% check used_samples
figure

for pd = 1:6
    subplot(1, 6, pd)
    tmp1 = used_samples(1, :, pd);
    tmp1 = cell2mat(tmp1);
    tmp1 = tmp1(:);
    h1 = histogram(tmp1, -40:2:40); 
    h1.FaceColor = 'g';
    h1.FaceAlpha = 0.3;
    hold on;
    tmp2 = used_samples(2, :, pd);
    tmp2 = cell2mat(tmp2);
    tmp2 = tmp2(:);
    h1 = histogram(tmp2, -40:2:40); 
    h1.FaceColor = 'r';
    h1.FaceAlpha = 0.3;
    hold on;
       
end
    