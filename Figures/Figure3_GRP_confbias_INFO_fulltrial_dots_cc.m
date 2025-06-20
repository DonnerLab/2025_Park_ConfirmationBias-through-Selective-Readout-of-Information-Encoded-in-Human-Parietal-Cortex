% relevant data: gFTII_all_%s_parcels_%s_nosh.mat, gFTII_all_%s_parcels_%s_sh.mat, 
% gFTMI_all_%s_parcels_%s_nosh.mat, gFTMI_all_%s_parcels_%s_sh.mat
% can be found in: https://www.fdr.uni-hamburg.de/record/16918

%% Parameters
conditions = {'all', 'attn', 'choice'};

load('subjects.mat')

cols1 = {};
cols1{1} = [120 120 120]./255;
cols1{2} = [0 0 0]./255;



%%  get data and plot information for areas
target_type = 'sp'; 
IIdata = cell(2,1);

indx = 3;
dotax = [1:6 8:13];
for c = 2:3
    file_name = sprintf('gFTII_all_%s_parcels_%s_nosh.mat', conditions{c}, target_type) 
    load(file_name)
    file_name = sprintf('gFTII_all_%s_parcels_%s_sh.mat', conditions{c}, target_type) 
    sh = load(file_name);
    clear tmp*
    for a = 1:length(indx)
        areas = glasser_group(indx(a));
        inda = ismember(gFTII(:, 2), areas);
        tmp = gFTII(inda, 1);
        tmpsh = sh.gFTII(inda, 1);
        data = {}; datash = {};
        for d = 1:length(dotax)
            for jj = 1:length(tmp)
                tmp2{jj,1} = tmp{jj}{dotax(d)};
                tmpsh2{jj,1} = tmpsh{jj}{dotax(d)};
            end            
            tmp3 = cell2mat(cellfun(@(x) permute(x, [3 2 1]), tmp2, 'UniformOutput', false));
            tmpsh3 = cell2mat(cellfun(@(x) permute(x, [3 2 1]), tmpsh2, 'UniformOutput', false));
            data{d} = permute(nanmean(tmp3, 1), [2 3 1]);
            datash{d} = permute(nanmean(tmpsh3, 1), [2 3 1]);            
        end        
        for d = 1:length(dotax)
            IIdata{c-1}{a}{d} = data{d}-datash{d};
        end
    end
end

%% plot II for areas
sf = 200;
etimep = -0.5:1/sf:4.5; % full length of ERF
d0s = {0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 2.70, 2.85, 3.00, 3.15, 3.30, 3.45};
% get samples of onset times
d0 = cellfun(@(x) find(etimep>x-0.5/sf & etimep<x+0.5/sf), d0s, 'UniformOutput', false); 
d0 = cell2mat(d0);

target_type = 'sp'; 

titles = {'s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','s_{11}','s_{12}'};


%% PLOT only 2nd interval
% 2nd interval
epochs = nan(6, 2);
for d = 7:12
    epochs(d-6, :) = [d0(d)-round(sf*0.1) d0(d) + round(sf*0.7)];    
end

xtime = -0.1:1/sf:0.7;
bs0 = find(xtime == 0);
vw = 2;
hw = 8; % cm
% create custom axes positions
wp = 0.55;
wpx = 0.75;
hfa = 0.08;
hla = 0.007; % distance from left
vfa = 0.18; % distance from bottom
vla = 0.015;
epochx = epochs(:, 2)-epochs(:, 1);
epochw = epochs(end, 2);
x_w = wpx*(epochx./epochw)/sum(epochx./epochw);
xspace = (1-wpx-hfa-hla)/(size(epochs, 1)-1); % spacing between axes
x_st = [];
x_st(1) = hfa; 
for d = 2:size(epochs, 1)
    tmp = cumsum(x_w(1:d-1));
    x_st(d) = hfa + tmp(end) + xspace*(d-1);
end
yspace = (1-wp-vfa-vla)/length(indx);
y_h = wp*1/length(indx);
y_st = [];
y_st(1) = vfa;
for a = 2:length(indx)
    y_st(a) = vfa + (y_h+yspace)*(a-1);
end
y_st = fliplr(y_st);

%% PLOT
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


vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 hw vw];
% annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('N = %d, Condition: %s, %s', length(subj), conditions{c}, msre), 'LineStyle', 'none', 'fontsize', 16)

ylm = [-0.001 0.004];

for a = 1:length(indx)
    [~,~, gname] = glasser_group(indx(a));    
    for d = 1:6
        ax = axes('Position', [x_st(d) y_st(a) x_w(d) y_h]);
        
        plot([etimep(epochs(d, 1)) etimep(epochs(d, 2))], [0 0], '-k', 'LineWidth', 0.1, 'Handlevisibility', 'off')
        hold on;
        plot([etimep(epochs(d, 1))+0.1 etimep(epochs(d, 1))+0.1], ylm, ':k', 'LineWidth', 0.1, 'Handlevisibility', 'off')

        % II
        for c = 2:3
            data1 = IIdata{c-1}{a}{d+6}(:, 1:round(sf*0.8));     
%             if d == 6
%                 data1 = data1(:, 1:round(sf*1.1)+1);
%             end
            tmp1 = data1;
           % smooth with gaussian
            w = gausswin(7, 0.2); % make it taper at the end
            w = repmat(w', size(tmp1, 1), 1);
            % pad data
            dstmp = cat(2, repmat(mean(tmp1(:,1:3),2),1,3), tmp1, repmat(mean(tmp1(:,end-2:end),2),1,3)); dsdat = [];
            ii = 1;
            for t = 1:length(xtime)
                dsdat(:, ii) = mean(w.*dstmp(:, t:t+6), 2);
                ii = ii+1;
            end
            
            % baseline for permutation test
            tmpb = IIdata{c-1}{a}{1};   
            % smoothe baseline data
            stmpb = cat(2, repmat(mean(tmpb(:,1:3),2),1,3), tmpb, repmat(mean(tmpb(:,end-2:end),2),1,3)); sbsdat = [];
            ii = 1;
            for t = 1:length(xtime)
                sbsdat(:, ii) = mean(w.*stmpb(:, t:t+6), 2);
                ii = ii+1;
            end
            baseline = mean(sbsdat(:, 1:bs0), 2);
            plot(etimep(epochs(d, 1):epochs(d, 2)), nanmean(dsdat-baseline, 1), 'color', cols1{c-1}, 'LineStyle', '-', 'LineWidth', 0.1);
            
            %%%%%%%%%% do permutation test between baseline and data %%%%%%%%%%
            baseline = repmat(baseline, 1, size(dsdat, 2));
            index = perm_test(dsdat, baseline, ARG, cfg);
            if ~isempty(index)
                plot(xtime(index), 10^-4-0.0004-0.0002*(cs-1),'.','Color',cols{cs},'markersize',3,'HandleVisibility','off');
            end
            hold on;
            
        end
        ylim(ylm)
        box off
        if d == 1
            hold on;
            plot([etimep(round(0.3*sf)) etimep(round(0.6*sf))], 3*10^-3.*[1 1], '-k', 'LineWidth', 5)
        end
        xlim([etimep(epochs(d,1)) etimep(epochs(d,2))])
        set(gca, 'tickdir', 'out')
        
        if d ==1
            set(gca, 'xtick', (etimep(epochs(d, 1))-0.05+0.1))
            ylabel(sprintf('%s\nII (bit)', gname))
        else
            set(gca, 'xtick',(etimep(epochs(d, 1))-0.05+0.1),'yticklabels', [])
        end
        if length(indx)~=1 && a == 1
            set(gca, 'xtick', [])
        end
    end
end

% save figure
fig_file = sprintf('II_timecourse_int2_%s_pca_est', target_type);
