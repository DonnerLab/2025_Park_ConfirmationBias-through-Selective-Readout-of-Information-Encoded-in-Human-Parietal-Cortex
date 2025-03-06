% relevant files: gSPII_all_%s_estim_rtco_parcels_bins.mat,
% gSPMI_all_%s_rtco_parcels_bins, gSPMI_all_%s_est_rtco_parcels_bins

%% Parameters

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};

sess_id = {[1 2 3 4], [1 2], [3 4]};

load('subjects.mat')

cind = {};
cind{1}{1} = 1;
cind{2}{1} = 1;
cind{2}{2} = 1;

sf = 200;
etimep = -0.1:1/sf:1;

d0 = find(etimep==0);


%% plot (up down mean) (estim)
ctype = 'est'; 
c = 1;
load(sprintf('gSPMI_all_%s_%s_rtco_parcels_bins.mat', conditions{c}, ctype), 'gMI', 'bins', 'gntr')

% colors
clear cols
if c == 2
    cols{1} = [228,26,28]./255; % int 1
    cols{2} = [80, 200, 120; % consistent
        242	133	0]./255; % inconsistent
    cols{3} = [0 0 0];
elseif c == 3
    cols{1} = [55,126,210]./255; % int 1
    cols{2} = [80, 200, 120; % consistent
        242	133	0]./255; % inconsistent
    cols{3} = [0 0 0];
end


%% PLOT
indx = [1:4];
yul = [-0.0002 0.005];

addpath(sprintf('%s/fieldtrip-20171001', mdir))
ft_defaults

ARG.pcrit = 0.01;  % for behav regression
ARG.Nperm = 2000;
ARG.T_Range =[-0.1 1]; % time range to display
% ARG.T_Range =[-0.1 0.7]; % time range to display

cfg = [];
cfg.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg.minsize = 4; % minimal cluster size
cfg.pval = 0.01; % threshold to select signifciant clusters
% cfg.df = length(subj)-1; %degrees of freedom. Only needed for effect size.


vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 28 8];
% annotation(vis, 'textbox', [0.1 0.9 0.1 0.1], 'string', sprintf('II PCA (sample): N = %d, %s, II(S,R,E)', length(subj), condt{c-1}), ... 
%     'LineStyle', 'none', 'fontsize', 16)
for a = 1:length(indx)
    [areas, arealb] = glasser_group(indx(a));
    aind = ismember(gMI(:, 2), areas);
    tmp1 = gMI{aind, 1}; 
    tmp1 = cellfun(@(x) permute(x, [3 2 1]), tmp1, 'UniformOutput', false);
    tmp1 = cell2mat(tmp1');
    tmp1 = nanmean(tmp1, 3);
    subplot(1, length(indx)+1, a)
    plot(etimep, tmp1-mean(tmp1(:, 1:d0), 2))
    title(arealb)
    if a == 1
        legend(string(bins))
    end
end
subplot(1, length(indx)+1, length(indx)+1)
bar(bins, mean(gntr./bins.^2))
ylim([0 300])
xlabel('n bins')
title('nSamples/bin')

% save figure
fig_file = sprintf('GRP_%s_timecourse_MI_sample_pca_rtco_int%d_nbins', conditions{c}, int);
% sname = fullfile(fig_dir, fig_file)
% hgexport(vis, sname)
