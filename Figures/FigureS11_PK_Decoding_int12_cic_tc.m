% relevant data: %s_%s_pm%d_btS_btE_btSE_abs.mat

conditions = {'all', 'attn', 'choice'};
cic = {'con','incon'};
condt = {'Cue','Choice'};
subj = get_sinfo(setdiff(1:34, [12 14 19 21]), 2);


clear cols
cols{1} = [10 184 0; 255 184 10; 100 100 100]./255;
cols{2} = [10 184 0; 255 184 10; 0 0 0]./255;


cfgb = [];
cfgb.cols = [cols{1}(1:2, :); 0 0 0];

%% plot

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

fgmu = 0;

pl = true;
if pl
    udpt = '_udpool';
elseif ~pl
    udpt = []; 
end
np = 2; npl = 2;
tnpl = round(12/npl);

% for permutation test
addpath('fieldtrip-20171001')
ft_defaults

sf = 160; 

xtime = -0.1:1/sf:0.7;
bs0 = find(xtime==0);
d1 = find(xtime == 0.1);
d2 = find(xtime >= 0.5-0.1/sf & xtime <= 0.5+0.1/sf);

ARG.pcrit = 0.05;  
ARG.Nperm = 2000;
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

BM_params
    
pds = [1 4];

%% plot decoding score
[areas,arean,areaid] = glasser_group(10);

fgmu = 0;
if fgmu > 0
    anot = sprintf('\\pm%d', fgmu);
elseif fgmu == 0
    anot = 'All';
end
vis2 = figure;
vis2.Units = 'centimeters';
vis2.Position = [5 5 32 24];
annotation(vis2, 'textbox', [0.08 0.9 0.1 0.1], 'string', sprintf('Decoding performance, %s, Gen\\mu: %s', arean, anot), 'LineStyle', 'none', 'fontsize', 10)
ylm = [-0.0125 0.06; -0.01 0.04; -0.03 0.03];

setxt = {'$$corr(S,\hat{S})$$','$$corr(E,\hat{E})$$','$$corr(\hat{S},E)$$'};
tht = [0.86 0.42];

absol = false;
bscorr = false;


for c = 1:2
    igSCR = cell(3, 6);

    load(sprintf('%s_%s_pm%d_btS_btE_btSE_abs.mat', condt{c}, areaid, fgmu), 'gbtS','gbtE','gbtSE','gCORR', 'ROI')
    
    tmp = gCORR(1:6);
    for parind = 1:length(areas)
        gSCR = tmp{parind};
        igSCR = cellfun(@(x,y) cat(4,x,y), igSCR, gSCR, 'UniformOutput', false);
    end
    clear gSCR
    % average across parcels
    SCR = cellfun(@(x) nanmean(x, 4), igSCR, 'UniformOutput', false); % S,E x times x subj
    annotation(vis2, 'textbox', [0.08 tht(c) 0.1 0.1], 'string', sprintf('%s', condt{c}), 'LineStyle', 'none', 'fontsize', 10)

    for pd = 1:tnpl
        for se = 1:3 % stimulus / estimation / stimulus^hat * estimation
            subplot(6,6,6*(se-1)+18*(c-1)+pd), 
            xline(0, ':k', 'HandleVisibility', 'off'), yline(0, ':k', 'HandleVisibility', 'off'), hold on;
            for cs = 1:2 % con/incon/all
                % smooth with gaussian
                if se == 3
                    if absol
                        z1 = permute(abs(SCR{cs,pd}(se, :, :)), [3 2 1]);
                    else
                        z1 = permute(SCR{cs,pd}(se, :, :), [3 2 1]);
                    end
                else
                    z1 = permute(SCR{cs,pd}(se, :, :), [3 2 1]);
                end
                w = gausswin(7, 0.2); % make it taper at the end
                w = repmat(w', size(z1,1), 1);
                % Pad data
                tmp1 = cat(2, repmat(mean(z1(:,1:3),2),1,3), z1, repmat(mean(z1(:,end-2:end),2),1,3)); dsz1 = [];
                ii = 1;
                for t = 1:length(xtime)
                    dsz1(:, ii) = mean(w.*tmp1(:, t:t+6), 2);
                    ii = ii+1;
                end
                
                ph = patch([xtime fliplr(xtime)], ...
                [mean(dsz1, 1)+std(dsz1, [], 1)/sqrt(size(dsz1, 1)) fliplr(mean(dsz1, 1)-std(dsz1, [], 1)/sqrt(size(dsz1, 1)))], ...
                    cols{c}(cs,:));
                ph.LineStyle = 'none';
                ph.FaceAlpha = 0.2; hold on;
                plot(xtime, mean(dsz1, 1),'color',cols{c}(cs,:)), hold on;

                baseline = permute(mean(SCR{cs,1}(se, 1:bs0, :),2), [3 2 1])*ones(1, length(xtime));
                index = perm_test(permute(SCR{cs,pd}(se, :, :), [3 2 1]),baseline, ARG,cfg);
                if ~isempty(index)
                    plot(xtime(index), -((ylm(se,2)-ylm(se,1))/20)*cs,'.','Color',cols{c}(cs,:),'markersize',3,'HandleVisibility','off');
                end
            end
            hold on;
            if cs ~= 3
                % compare between con/incon
                if se == 3
                    if absol
                        aSCR = cellfun(@(x) abs(x), SCR, 'UniformOutput', false); 
                        index1 = perm_test(permute(aSCR{1,pd}(se, :, :), [3 2 1]),permute(aSCR{2,pd}(se, :, :), [3 2 1]), ARG,cfg);
                    else
                        index1 = perm_test(permute(SCR{1,pd}(se, :, :), [3 2 1]),permute(SCR{2,pd}(se, :, :), [3 2 1]), ARG,cfg);
                    end
                else
                    index1 = perm_test(permute(SCR{1,pd}(se, :, :), [3 2 1]),permute(SCR{2,pd}(se, :, :), [3 2 1]), ARG,cfg);
                end
                if ~isempty(index1)
                    plot(xtime(index1), -((ylm(se,2)-ylm(se,1))/20)*0,'.','Color','k','markersize',3,'HandleVisibility','off');
                end   
            end
            xlabel('time (s)'); ylim(ylm(se,:));
            if pd == 1,ylabel(setxt{se}, 'Interpreter', 'Latex'),end
            if c == 1 && se == 1, title(psmpt{pd}), end
            set(gca, 'tickdir', 'out')
        end 
    end
end

allAxes = findall(vis2,'type','axes');
for a = 1:length(allAxes)
    allAxes(a).Position(2) = allAxes(a).Position(2)-0.01;
end
sname = sprintf('GRP_RegSamp_%s_%s_allParcels_SCORE_int12_pm%d_udpool_np%d_cic', conditions{c}, areaid, fgmu, npl);
saveas(vis2, sname, 'png')
exportgraphics(vis2, [sname '.eps'],'ContentType','vector',...
       'BackgroundColor','none')



%% plot bar graphs of various ways of contrasting and averaging correlation values
%% PLOT for manuscript
indx = [3 10];


hw1 = 8; vw1 = 8;

ylm = [-0.01 0.03]; 

cfg1 = [];
cfg1.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg1.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg1.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg1.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg1.minsize = 1; % minimal cluster size
cfg1.pval = 0.05; % threshold to select signifciant clusters
cfg1.df = length(subj)-1; %degrees of freedom. Only needed for effect size.
cfg1.smooth = false;
cfg1.ttest = false;
cfg1.perm = 'unpaired';
cfg1.nperm = 5000;

cfgc = [];
cfgc.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfgc.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfgc.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfgc.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfgc.minsize = 1; % minimal cluster size
cfgc.pval = 0.05; % threshold to select signifciant clusters
cfgc.df = length(subj)-1; %degrees of freedom. Only needed for effect size.
cfgc.smooth = false;
cfgc.ttest = false;
cfgc.perm = 'zero';
cfgc.nperm = 5000;

np = 2;
for c = 1:2
    vis1 = figure;
    vis1.Units = 'centimeters';
    vis1.Position = [2 2 hw1 vw1];
    [x_st1, y_st1, x_w1, y_h1] = multi_axes(hw1, vw1, length(pds), 2, [1 1], [0.18 0.005 0.1 0.1], [0.05 0.1]);
    arealb = cell(2,1);
    for a = 1:length(indx)
        [areas,arean,arealb{a}] = glasser_group(indx(a));
        igSCR = cell(3, 6); igSCRneg = cell(3, 6);
        for parind = 1:length(areas)
            load(sprintf('%s_pm%d_btS_btE_btSE_%s_abs.mat', condt{c}, fgmu, areas{parind}), 'gSCR','gSCRneg')
            igSCR = cellfun(@(x,y) cat(4,x,y), igSCR, gSCR, 'UniformOutput', false);
            igSCRneg = cellfun(@(x,y) cat(3,x,y), igSCRneg, gSCRneg, 'UniformOutput', false);
        end
        clear gSCR
        % SCR{cs,pd}(se, :, :, :), se = 1:3 % stimulus / estimation / stimulus^hat * estimation 
        % se × time × cp × subj
        % average across parcels
        SCR = cellfun(@(x) nanmean(x, 4), igSCR, 'UniformOutput', false); % S,E x times x subj
        SCRneg = cellfun(@(x) nanmean(x, 3), igSCRneg, 'UniformOutput', false); % S,E x times x subj
        z3 = cellfun(@(x) permute(x(3, :, :), [3 2 1]), SCR, 'UniformOutput', false);
        for pd = 1:length(pds)
            %%%%% contrast between con and incon after averaging within time window %%%%%  
            ax3 = axes(vis1, 'Units', 'centimeters', 'Position', [x_st1(pd), y_st1(a), x_w1, y_h1]); hold all;
            con = abs(mean(z3{1, pds(pd)}(:, d1:d2), 2)); incon = abs(mean(z3{2, pds(pd)}(:, d1:d2), 2));
            pvalcic = perm_test_behav(con, incon, cfg1);
            conts = con - incon;
            pvalc = perm_test_behav(conts, incon, cfgc);
            cfgb.pval = pvalcic; cfgb.ylm = ylm(1,:);cfgb.tw = [d1 d2];cfgb.pvalc = pvalc;
            plot_inset_bars3(con, incon, conts, cfgb);
            set(gca, 'ytick', [])
            if pd == 1
                set(gca, 'ytick', [ylm(1) 0 ylm(2)])
                ylabel(sprintf('%s', arean))
            else
                set(gca, 'ytick', [])
            end
            if a == 1
                title(psmpt{pds(pd)})
            end
        end
        set(findall(vis1, '-property','FontName'),'FontName','Helvetica')
        set(findall(vis1, '-property','FontSize'),'FontSize',8)
    end
    sname2 = sprintf('GRP_Decoding_%s_%s_%s_pm%d_np%d_bar_contrast_muabs', arealb{1},arealb{2}, condt{c}, fgmu, np);
    exportgraphics(vis1, [sname2 '.eps'],'ContentType','vector',...
        'BackgroundColor','none')
end
%% PLOT
indx = [3 10];


hw1 = 18; vw1 = 8;

ylm = [-0.05 0.05]; 

cfg1 = [];
cfg1.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg1.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg1.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg1.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg1.minsize = 1; % minimal cluster size
cfg1.pval = 0.05; % threshold to select signifciant clusters
cfg1.df = length(subj)-1; %degrees of freedom. Only needed for effect size.
cfg1.smooth = false;
cfg1.ttest = false;
cfg1.perm = 'unpaired';
cfg1.nperm = 5000;

cfgc = [];
cfgc.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfgc.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfgc.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfgc.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfgc.minsize = 1; % minimal cluster size
cfgc.pval = 0.05; % threshold to select signifciant clusters
cfgc.df = length(subj)-1; %degrees of freedom. Only needed for effect size.
cfgc.smooth = false;
cfgc.ttest = false;
cfgc.perm = 'zero';
cfgc.nperm = 5000;

np = 2;
for c = 1:2
    vis1 = figure;
    vis1.Units = 'centimeters';
    vis1.Position = [2 2 hw1 vw1];
    [x_st1, y_st1, x_w1, y_h1] = multi_axes(hw1, vw1, 3*length(pds), 2, [1 1], [0.1 0.005 0.1 0.1], [0.1 0.1]);
    arealb = cell(2,1);
    for a = 1:length(indx)
        [areas,arean,arealb{a}] = glasser_group(indx(a));
        igSCR = cell(3, 6); igSCRneg = cell(3, 6);
        for parind = 1:length(areas)
            load(sprintf('%s_pm%d_btS_btE_btSE_%s_abs.mat', condt{c}, fgmu, areas{parind}), 'gSCR','gSCRneg')
            igSCR = cellfun(@(x,y) cat(4,x,y), igSCR, gSCR, 'UniformOutput', false);
            igSCRneg = cellfun(@(x,y) cat(3,x,y), igSCRneg, gSCRneg, 'UniformOutput', false);
        end
        clear gSCR
        % SCR{cs,pd}(se, :, :, :), se = 1:3 % stimulus / estimation / stimulus^hat * estimation 
        % se × time × cp × subj
        % average across parcels
        SCR = cellfun(@(x) nanmean(x, 4), igSCR, 'UniformOutput', false); % S,E x times x subj
        SCRneg = cellfun(@(x) nanmean(x, 3), igSCRneg, 'UniformOutput', false); % S,E x times x subj
        z3 = cellfun(@(x) permute(x(3, :, :), [3 2 1]), SCR, 'UniformOutput', false);
        for pd = 1:length(pds)
            %%%%% basic contrast between con/incon %%%%%            
            ax1 = axes(vis1, 'Units', 'centimeters', 'Position', [x_st1(3*(pd-1)+1), y_st1(a), x_w1, y_h1]); hold all;
            con = z3{1, pds(pd)}; incon = z3{2, pds(pd)};
            [~, pvalcic] = perm_test(con, incon, ARG,cfg);
            conts = con - incon; conts = mean(conts(:, d1:d2), 2);
            pvalc = perm_test_behav(conts, incon, cfgc);
            cfgb.pval = pvalcic; cfgb.ylm = ylm(1,:);cfgb.tw = [d1 d2]; cfgb.pvalc = pvalc;
            plot_inset_bars3(con, incon, conts, cfgb)
            if pd == 1
                set(gca, 'ytick', [ylm(1) 0 ylm(2)])
                ylabel(sprintf('%s', arean))
            else
                set(gca, 'ytick', [])
            end
            
            %%%%% contrast between con and negative incon %%%%%  
            ax2 = axes(vis1, 'Units', 'centimeters', 'Position', [x_st1(3*(pd-1)+2), y_st1(a), x_w1, y_h1]); hold all;
            con = SCRneg{1, pds(pd)}(3, :)'; incon = SCRneg{2, pds(pd)}(3, :)';
            pvalcic = perm_test_behav(con, incon, cfg1);
            conts = SCRneg{3, pds(pd)}(3, :)';
            pvalc = perm_test_behav(conts, incon, cfgc);
            cfgb.pval = pvalcic; cfgb.ylm = ylm(1,:);cfgb.tw = [d1 d2]; cfgb.pvalc = pvalc;
            plot_inset_bars3(con, incon, conts, cfgb);
            set(gca, 'ytick', [])
            
            %%%%% contrast between con and incon after averaging within time window %%%%%  
            ax3 = axes(vis1, 'Units', 'centimeters', 'Position', [x_st1(3*(pd-1)+3), y_st1(a), x_w1, y_h1]); hold all;
            con = abs(mean(z3{1, pds(pd)}(:, d1:d2), 2)); incon = abs(mean(z3{2, pds(pd)}(:, d1:d2), 2));
            pvalcic = perm_test_behav(con, incon, cfg1);
            conts = con - incon;
            pvalc = perm_test_behav(conts, incon, cfgc);
            cfgb.pval = pvalcic; cfgb.ylm = ylm(1,:);cfgb.tw = [d1 d2];cfgb.pvalc = pvalc;
            plot_inset_bars3(con, incon, conts, cfgb)
            set(gca, 'ytick', [])
        end
        set(findall(vis1, '-property','FontName'),'FontName','Helvetica')
        set(findall(vis1, '-property','FontSize'),'FontSize',8)
    end
    sname2 = sprintf('GRP_Decoding_%s_%s_%s_pm%d_np%d_bar_contrasts', arealb{1},arealb{2}, condt{c}, fgmu, np);
    exportgraphics(vis1, [sname2 '.eps'],'ContentType','vector',...
        'BackgroundColor','none')
end

