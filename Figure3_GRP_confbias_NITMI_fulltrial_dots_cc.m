
%% Parameters
conditions = {'all', 'attn', 'choice'};

load('subjects.mat')

cols1 = {};
cols1{1}{1} = [228,26,28]./255;
cols1{1}{2} = [10 184 0; 255 184 10]./255;
cols1{2}{1} = [55,126,210]./255;
cols1{2}{2} = [10 184 0; 255 184 10]./255;
cols1{3}{1} = [0 0 0]./255;

%%  get data and plot MI for areas
target_type = 'sp'; 
Idata = cell(2,1);
Idata1 = cell(2,1);
Idata2 = cell(2,1);
indx = 3;
dotax = [1:6 8:13];
for c = 2:3
    file_name = sprintf('%s/P03/MI/gFTMI_all_%s_parcels_%s_nosh.mat', mdir, conditions{c}, target_type) % cue/choice number of trials the same
    load(file_name)
    file_name = sprintf('%s/P03/MI/gFTMI_all_%s_parcels_%s_sh.mat', mdir, conditions{c}, target_type) % cue/choice number of trials the same
    sh = load(file_name);
    clear tmp*
    for a = 1:length(indx)
        areas = glasser_group(indx(a));
        inda = ismember(gFTMI(:, 2), areas);
        tmp = gFTMI(inda, 1);
        tmpsh = sh.gFTMI(inda, 1);
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
            Idata{c-1}{a}{d} = data{d}-datash{d};
            Idata1{c-1}{a}{d} = data{d};
            Idata2{c-1}{a}{d} = datash{d};
        end
    end
end


%% plot MI for areas
sf = 200;
etimep = -0.5:1/sf:4.5; % full length of ERF
d0s = {0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 2.70, 2.85, 3.00, 3.15, 3.30, 3.45};
% get samples of onset times
d0 = cellfun(@(x) find(etimep>x-0.5/sf & etimep<x+0.5/sf), d0s, 'UniformOutput', false); 
d0 = cell2mat(d0);

% epochs = nan(13, 2);
% % for first dot, include 500 ms baseline period
% epochs(1, :) = [1 d0(1)+round(sf*1)]; % baseline period
% for d = 2:6
%     epochs(d, :) = [d0(d)-round(sf*0.1) d0(d) + round(sf*1)];    
% end
% % between int1 and int2
% epochs(7, :) = [d0(6) d0(7)];
% for d = 8:12
%     epochs(d, :) = [d0(d-1)-round(sf*0.1) d0(d-1) + round(sf*1)];
% end
% % for the 12th dot include post stim period
% epochs(13, :) = [d0(12)-round(sf*0.1) length(etimep)]; % after second interval
% % epochs(epochs(:, 2)>length(etimep), 2) = length(etimep);

% without break
epochs = nan(12, 2);
% for first dot, include 500 ms baseline period
for d = 1:12
    epochs(d, :) = [d0(d)-round(sf*0.1) d0(d) + round(sf*0.7)];    
end

load(sprintf('%s/P03/behav/gRT_simple.mat', mdir), 'gRT')


%% plot together for each dot

% titles = {'{s^{I}}_1','{s^{I}}_2','{s^{I}}_3','{s^{I}}_4','{s^{I}}_5','{s^{I}}_6','{<s^{I}>}_{1-6}','{s^{II}}_1','{s^{II}}_2','{s^{II}}_3','{s^{II}}_4','{s^{II}}_5','{s^{II}}_6'};
titles = {'s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','s_{11}','s_{12}'};

%% Plot
vw = 4;
hw = 40; % cm
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

% RT offset for interm and estimation;
rto = 0.95; 

% for permutation test

addpath(sprintf('%s/fieldtrip-20171001', mdir))
ft_defaults

ARG.pcrit = 0.01;  
ARG.Nperm = 2000;
% ARG.T_Range =[-0.1 1]; % time range to display

cfg = [];
cfg.critvaltype = 'par';%'par'; 'prctile' % type of threshold to apply. Usual 'par'
cfg.critval = abs(tinv(ARG.pcrit/2,length(subj)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
cfg.conn = 4; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
cfg.clusterstatistic = 'maxsum'; %'max' 'maxsize' 'maxsum'
cfg.minsize = 2; % minimal cluster size
cfg.pval = 0.01; % threshold to select signifciant clusters
cfg.df = length(subj)-1; %degrees of freedom. Only needed for effect size.

%% PLOT

vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 hw vw];
% annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('N = %d, Condition: %s, %s', length(subj), conditions{c}, msre), 'LineStyle', 'none', 'fontsize', 16)

if strcmp(target_type, 'est')
    ylm = [-0.0007 0.007];
    ep = -3;
elseif strcmp(target_type, 'sp')
    ylm = [-0.002 0.02];
    ep = -2;
end
for a = 1:length(indx)
    [~,~, gname] = glasser_group(indx(a));    
    for d = 1:length(dotax)
        ax = axes('Position', [x_st(d) y_st(a) x_w(d) y_h]);
        
        plot([etimep(epochs(d, 1)) etimep(epochs(d, 2))], [0 0], '-k', 'LineWidth', 1, 'Handlevisibility', 'off')
        hold on;
%         if d == 1
%             plot([etimep(epochs(d, 1))+0.5 etimep(epochs(d, 1))+0.5], ylm, ':k', 'LineWidth', 0.5, 'Handlevisibility', 'off')
%         elseif ismember(d, [2:6 8:13])
            plot([etimep(epochs(d, 1))+0.1 etimep(epochs(d, 1))+0.1], ylm, ':k', 'LineWidth', 1, 'Handlevisibility', 'off')
%         end
%         if d == 7
%             % plot RT
%             rts = rto + mean([mean(gRT{2}(setdiff(1:34, [12 14 19 21]), 1), 1) mean(gRT{3}(setdiff(1:34, [12 14 19 21]), 1), 1)]);
%             plot([rts rts], [0 ylm(end)], 'g', 'LineWidth', 1, 'Handlevisibility', 'off')
%             hold on;
%         end
        % MI
        for c = 2:3
%             if d == 7 && c == 2
%                 % plot attention cue
%                 acue_x = [1.95 1.95+0.5 1.95+0.5 1.95];
%                 acue_y = [ylm(1) ylm(1) ylm(2) ylm(2)];
%                 ph = patch(acue_x, acue_y, [0 0 0], 'Handlevisibility', 'off');
%                 ph.EdgeColor = 'none';
%                 ph.FaceAlpha = 0.2;
%             end

            data1 = Idata{c-1}{a}{d};   
            if d == 1
                data1 = data1(:, round(sf*0.45)+1:end-round(sf*0.3));
            else
                data1 = data1(:, 1:round(sf*0.8)+1);
            end
            tmp1 = data1;
            
            % baseline for permutation test
            baseline = Idata{c-1}{a}{1};   
%             baseline = baseline';
            baseline = mean(baseline(:, 1:d0(1)), 2);
            
%             tmp1m = mean([tmp1(:, 1:3:end-2); tmp1(:, 2:3:end-1); tmp1(:, 3:3:end)], 1);
            plot(etimep(epochs(d, 1):epochs(d, 2)), nanmean(tmp1-baseline, 1), 'color', cols1{c-1}{1}(1, :), 'LineStyle', '-', 'LineWidth', 0.1);
%             xdat = etimep(epochs(d, 1):3:epochs(d, 2));
%             xdat = xdat(1:length(tmp1m));
%             plot(xdat, tmp1m, 'color', cols1{c-1}{1}(1, :), 'LineStyle', '-', 'LineWidth', 0.5);
            ax = gca;
            ax.YAxis.Exponent = ep;
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
            detimep = etimep(epochs(d, 1):epochs(d, 2));
            if ~isempty(index)
                plot(detimep(index), 0.1*ylm(1)+0.4*ylm(1)*(c-2),'.','Color',cols1{c-1}{1}(1, :),'markersize',1,'HandleVisibility','off');
            end
            hold on;
        end
        ylim(ylm)
        box off
        
%         if a == 1
%             title(titles{d})
%         elseif a == length(indx)
%             xlabel('time (s)')
%         end
        if d == 1
            hold on;
            plot([etimep(round(0.65*sf)) etimep(round(1.05*sf))], 1.8*10^-2.*[1 1], '-k', 'LineWidth', 5)            
            ylabel(sprintf('%s\nI (bit)', gname))
        end
        xlim([etimep(epochs(d,1)) etimep(epochs(d,2))])
        set(gca, 'tickdir', 'out')
        
        if d ==1
            set(gca, 'xtick', (etimep(epochs(d, 1))-0.05+0.1))
        else
            set(gca, 'xtick',(etimep(epochs(d, 1))-0.05+0.1),'yticklabels', [])
        end
        if length(indx)~=1 && a == 1
            set(gca, 'xtick', [])
        end
    end
end


% % save figure
fig_file = sprintf('Figure2_MI_timecourse_fulltrial_%s_pca_%d', target_type, sum(indx));
hgexport(vis, fullfile(fig_dir, fig_file))



