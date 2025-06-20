% relevant data: gFTII_%s_parcels_pm%d_%s_np%d_msh.mat, gFTII_%s_parcels_pm%d_%s_np%d_msh_PC1.mat
% gFTMI_%s_parcels_pm%d_%s_np%d.mat, gFTMI_%s_parcels_pm%d_%s_np%d_PC1.mat
% can be found here: https://www.fdr.uni-hamburg.de/record/16918

conditions = {'all', 'attn', 'choice'};
cic = {'con','incon'};
condt = {'Cue','Choice'};
subj = get_sinfo(setdiff(1:34, [12 14 19 21]), 2);

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
cfg.smooth = false;
cfg.ttest = false;
cfg.perm = 'unpaired';


pds = [1 4];
indx = 3;

% get good subject index
load('gminsp_un_pooled.mat')

hw = 6;
vw = 4;

%% PLOT
% if PC1, 
% udpt = 'udpool' ;
udpt = '_udpool';

target_type = []; % '_est'; % '_sp'


% this is for PC1
if isempty(target_type)
    iflg = 'II';
    infot = 'II(S;R;E)';
    ylm = 10^-3.*[-0.7 4];
    mnsp = 81;
elseif strcmp(target_type, '_sp')
    iflg = 'MI';
    ylm = 10^-3.*[-2 10];
    infot = 'I(S;R)';
    mnsp = 27;
elseif strcmp(target_type, '_est')
    iflg = 'MI';
    ylm = 10^-3.*[-1.6 7];
    infot = 'I(R;E)';
    mnsp = 27;
end

% this is for all PC averaged
% if isempty(target_type)
%     iflg = 'II';
%     infot = 'II(S;R;E)';
%     ylm = 10^-3.*[-0.7 3];
%     mnsp = 81;
% elseif strcmp(target_type, '_sp')
%     iflg = 'MI';
%     ylm = 10^-3.*[-1.6 7];
%     infot = 'I(S;R)';
%     mnsp = 27;
% elseif strcmp(target_type, '_est')
%     iflg = 'MI';
%     ylm = 10^-3.*[-1.6 7];
%     infot = 'I(R;E)';
%     mnsp = 27;
% end

hw1 = 4; vw1 = 6;

clear gFTII gFTMI
for c = 2:3
    vis = figure;
    vis.Units = 'centimeters';
    vis.Position = [2 2 hw vw];
    % annotation(vis, 'textbox',[0.05 0.9 0.1 0.1], 'String', ...
    %     sprintf('gen means: \\pm(%d:0), data size > 81', fgmu(f)), 'LineStyle', 'none')
    [x_st, y_st, x_w, y_h] = multi_axes(hw, vw, length(pds), 1, [1 1], [0.13 0.01 0.1 0.2], [0.06 0.001]);
    
    minsp = gminsp_udpooled{fgmu+2, c-1};
    load(sprintf('gFT%s_%s_parcels_pm%d_%s_np%d%s_msh_PC1.mat', iflg, iflg, conditions{c}, fgmu, udpt, np, target_type))
    [areas,~,areaid] = glasser_group(indx);
    % if it is for PC1
    gFTII = gFTIIpc1;
    if strcmp(iflg, 'II')
        inda = ismember(gFTII(:, 2), areas);
        data = gFTII(inda, 1);
        npd = size(gFTII{1,1}, 3);
    elseif strcmp(iflg, 'MI')
        inda = ismember(gFTMI(:, 2), areas);
        data = gFTMI(inda, 1);
        npd = size(gFTMI{1,1}, 3);
    end
    
    pvalcic = [];
    for pd = 1:length(pds)
        ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st(pd), y_st, x_w, y_h]); hold all;
        dat = cell(2, 1); smdat = cell(2, 1);
        minsp = min(gminsp_udpooled{fgmu+2, c-1}(:, :, pds(pd)), [], 2);
        goodsub = minsp >= mnsp;
        xline(0, ':', 'HandleVisibility', 'off'); yline(0, '-', 'HandleVisibility', 'off')
        for cs = 1:2
            tmp = cellfun(@(x) permute((x(cs, :, pds(pd), goodsub)), [1 2 4 3]), data, 'UniformOutput', false);
            tmp = cell2mat(tmp);
            tmp = permute(nanmean(tmp, 1), [3 2 1]); dat{cs} = tmp;
            
            % smooth with gaussian
            w = gausswin(7, 0.2); % make it taper at the end
            w = repmat(w', size(tmp, 1), 1);
            % pad data
            dstmp = cat(2, repmat(mean(tmp(:,1:3),2),1,3), tmp, repmat(mean(tmp(:,end-2:end),2),1,3)); dsdat = [];
            ii = 1;
            for t = 1:length(xtime)
                dsdat(:, ii) = mean(w.*dstmp(:, t:t+6), 2);
                ii = ii+1;
            end
            smdat{cs} = dsdat;
            ph = patch([xtime fliplr(xtime)], ...
                [mean(dsdat, 1)+std(dsdat, [], 1)/sqrt(size(tmp, 1)) fliplr(mean(dsdat, 1)-std(dsdat, [], 1)/sqrt(size(tmp, 1)))], ...
                cols{cs});
            ph.LineStyle = 'none';
            ph.FaceAlpha = 0.2;
            hold on;
            plot(xtime, mean(dsdat, 1), 'color', cols{cs}, 'LineWidth', 1);
            
            % baseline for permutation test
            tmpb = cellfun(@(x) permute((x(cs, :, 1, goodsub)), [1 2 4 3]), data, 'UniformOutput', false);
            tmpb = cell2mat(tmpb);
            tmpb = permute(nanmean(tmpb, 1), [3 2 1]);
            % smoothe baseline data
            stmpb = cat(2, repmat(mean(tmpb(:,1:3),2),1,3), tmpb, repmat(mean(tmpb(:,end-2:end),2),1,3)); sbsdat = [];
            ii = 1;
            for t = 1:length(xtime)
                sbsdat(:, ii) = mean(w.*stmpb(:, t:t+6), 2);
                ii = ii+1;
            end
            baseline = mean(sbsdat(:, 1:d0), 2);
            
            %%%%%%%%%% do permutation test between baseline and data %%%%%%%%%%
            baseline = repmat(baseline, 1, size(dsdat, 2));
            index = perm_test(dsdat, baseline, ARG, cfg);
            if ~isempty(index)
                plot(xtime(index), 10^-4-0.0004-0.0002*(cs-1),'.','Color',cols{cs},'markersize',3,'HandleVisibility','off');
            end
            hold on;
        end
        clear mask index PosClus NegClus
        %%%%%%%%%% do permutation test between con/incon %%%%%%%%%%
        [index, pvalcic(pd)] = perm_test(smdat{1}, smdat{2}, ARG, cfg);
        if ~isempty(index)
            plot(xtime(index), 10^-4-0.0002,'.','Color','k','markersize',3,'HandleVisibility','off');
        end
        clear mask index PosClus NegClus
        
        hold on;
        xlim([xtime(1) xtime(end)]); ylim(ylm)
        
        xlabel('Time (s)')
        if pd == 1
            ylabel(sprintf('%s (bit)',infot))
        end
        th = title(sprintf('%s',psmpt{pds(pd)}));
        th.Position = [th.Position(1) th.Position(2)+0.0007 th.Position(3)];
        %%%%%%%%%%% test 0.1-0.5 window %%%%%%%%%%%
        if pvalcic(pd)<0.05
            ph = patch([0.1 0.5 0.5 0.1], [ylm(1) ylm(1) ylm(2) ylm(2)], 'k', 'HandleVisibility', 'off');
            ph.FaceAlpha = 0.08;
            ph.LineStyle = 'none';
            uistack(ph, 'bottom')
            text(0.1, ylm(2)*0.98, sprintf('p=%1.4f',pvalcic(pd)));
        end
        box off
        if c == 3 && pd == npd
            legend('Consistent','Inconsistent','Location','NorthEastOutside')
            legend boxoff
        end
        ax = gca;
        set(ax, 'xtick', [0 0.2 0.4 0.6], 'tickdir', 'out', 'layer', 'top')
        ax.Position = [ax.Position(1)-0.05 ax.Position(2)-0.015 ax.Position(3)*1.02 ax.Position(4)*0.8];
    end
    pvalcic(pvalcic==0) = 10^-30; h1 = fdr(pvalcic, 0.05);
    set(findall(vis, '-property','FontName'),'FontName','Helvetica')
    set(findall(vis, '-property','FontSize'),'FontSize',7)
    sname = sprintf('GRP_%s_%s_%s_pm%d%s_np%d%s_msh_PC1', iflg, areaid, conditions{c}, fgmu, udpt, np, target_type);
    saveas(vis, sname, 'png')
%     exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%                'BackgroundColor','none')
end




