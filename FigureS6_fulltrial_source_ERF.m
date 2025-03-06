% relevant data: ERF_source_%s_CC_fulltrial_abs.mat

%% Parameters
conditions = {'attn', 'choice'};

load('subjects.mat')

cols = [120 120 120; 0 0 0]./255; 
% plot together for each dot
titles = {'s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','s_{11}','s_{12}'};

sf = 400;
etimep = -0.5:1/sf:4.5; % full length of ERF

d0s = {0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 2.70, 2.85, 3.00, 3.15, 3.30, 3.45};
d0 = cellfun(@(x) find(etimep>x-0.1/sf & etimep<x+0.1/sf), d0s, 'UniformOutput', false);
d0 = cell2mat(d0);


% number of samples to pool
npl = 1;
tnpl = round(12/npl);
sess_id = {[1 2 3 4], [1 2], [3 4]}; % 1,2 is Cue(attn), 3,4 is Choice sessions


%% Plot
% for permutation test
addpath('fieldtrip-20171001')
ft_defaults

ARG.pcrit = 0.01;  
ARG.Nperm = 2000;

cfg = [];
cfg.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg.minsize = 2; % minimal cluster size
cfg.pval = 0.01; % threshold to select signifciant clusters
cfg.df = length(subj)-1; %degrees of freedom. Only needed for effect size.


%% plot full trial
load('gRT_simple.mat');
rto = 0.95;
load(sprintf('ERF_source_%s_CC_fulltrial_abs.mat', areaid),'cdata')
condt = {'Cue','Choice'};

hw = 15;
vw = 7;
vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 hw vw];
% annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('N = %d, Condition: %s, %s', length(subj), conditions{c}, msre), 'LineStyle', 'none', 'fontsize', 16)

ylm = 10^-12.*[-0.01 0.05];

donset = [0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 2.70, 2.85, 3.00, 3.15, 3.30, 3.45];
xdonset = repmat(donset', 1, 2);
ydonset = repmat(ylm, 12, 1);

for c = 1:2
    subplot(2, 1, c)
    if c == 1
        % plot attention cue
        acue_x = [1.95 1.95+0.5 1.95+0.5 1.95];
        acue_y = [ylm(1) ylm(1) ylm(2) ylm(2)];
        ph = patch(acue_x, acue_y, [0 0 0], 'Handlevisibility', 'off');
        ph.EdgeColor = 'none';
        ph.FaceAlpha = 0.15;
        hold on;
    end
    yline(0,':k', 'LineWidth', 0.1, 'Handlevisibility', 'off')
    hold on;
     % plot RT
    rts = rto + mean(gRT{c+1}(setdiff(1:34, [12 14 19 21]), 1), 1);
    plot([rts rts], ylm, 'g', 'LineWidth', 2, 'Handlevisibility', 'off'), hold on; 
    plot([etimep(d0);etimep(d0)], ylm, ':r', 'LineWidth', 0.1, 'Handlevisibility', 'off')
    data1 = cdata{c};
    tmp1 = permute(data1, [2 1]); % subjects x time
    % baseline for permutation test
    baseline = mean(tmp1(:, 1:d0(1)), 2);
    
    plot(etimep, nanmean(tmp1-baseline, 1), 'color', cols(c, :), 'LineStyle', '-', 'LineWidth', 0.1);
    
    %%%%%%%%%% do permutation test between baseline and data %%%%%%%%%%
    tmp2 = repmat(baseline, 1, size(tmp1, 2));
    tmp = tmp1-tmp2;
    tval = sqrt(length(subj))*(mean(tmp, 1)./std(tmp, [], 1));
    tmpc = cat(1, tmp1, tmp2);
    tval_boot = zeros(size(tmp,2),ARG.Nperm);
    for b=1:ARG.Nperm
        % randomly switch labels between baseline and data
        rs = cat(1, ones(length(subj), 1), 2*ones(length(subj), 1));
        rs = logical(rs(randperm(length(rs)))-1);
        tmp1v = tmpc(rs, :);
        tmp2v = tmpc(~rs, :);
        tmpv = tmp1v-tmp2v;
        tval_boot(:,b) = sqrt(length(subj))*(mean(tmpv, 1)./std(tmpv, [], 1));
    end
    [PosClus,NegClus] = eegck_clusterstats(cfg,tval',tval_boot);
    [mask,index] = eegck_stripclusters(PosClus,NegClus,[size(tval, 2) 1]);
    if ~isempty(index)
        plot(etimep(index), 0,'.','Color',cols(c,:),'markersize',1,'HandleVisibility','off');
    end
    hold on;
    title(condt{c})
    ylabel('fT')
    ylim(ylm)
    box off

    xlim([etimep(1) etimep(end)])
    set(gca, 'tickdir', 'out', 'xtick',etimep(d0),'xticklabels', titles)    
end



% % save figure
% fig_file = sprintf('FigureS5_ERF_timecourse_fulltrial_source_%s_abs', areaid);
% hgexport(vis, fullfile(fig_dir, fig_file))

