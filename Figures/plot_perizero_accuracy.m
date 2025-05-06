%% Parameters
pc_local = 2; % 1 is cluster
if pc_local == 2
    mdir = '/Users/hpark/Servers/mountpoint1';
%     mdir = '/Users/hpark/Research_wAlan/';
    fig_dir = '/Users/hpark/Research_wAlan/_manuscript/figures';
    dat_dir = '/Users/hpark/Research_wAlan/data';
elseif pc_local == 1
    mdir = '/home/hamepark';
    fig_dir = '/home/hamepark/P03/fig';
end
subj = get_sinfo(setdiff(1:34, [12 14 19 21]), 2);

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};

conditions = {'attn', 'choice'};

sess_id = {[1 2 3 4], [1 2], [3 4]}; % 1,2 is Cue(attn), 3,4 is Choice sessions

fgmu = 4;
%% get peri-zero trials accuracy and do stats between Cue and Choice

cols = [150 150 150; 0 0 0]./255;
gaccu = cell(2, 1); saccu = cell(2, 1);
for s = 1:length(subj)
    fname = sprintf('%s/P03/behav/%s_4decode_rtco.mat', mdir, subj{s})
    beh = load(fname); % behavioral data
    for c = 1:2
        
        
        clear tind vind
        tind = ismember(beh.sess, sess_id{c+1});
        genmus = beh.gen_means(tind);
        if c == 1
            interm = beh.interm_cues(tind);
        elseif c == 2
            interm = beh.binary_choices(tind);
        end
        
        % valid trials based on interm button press
        ind1 = beh.binary_choices(tind)~=99; % missed interm button press (includes both intervals)
        % valid trials based on estimation
        ind2 = ~isnan(beh.estimations(tind)); % missed estimation
        
        if fgmu > 0
            % valid trials2: based on generative mean:
            zgmu = [-fgmu:fgmu];
            zind = ismember(genmus, zgmu);
        elseif fgmu == 0 % use all trials
            zind = ~isnan(genmus(tind));
        end
        pztrls = zind & ind1 & ind2;
        % for the correct trials
        correcttrls = interm(pztrls).*genmus(pztrls)>0;
        gaccu{c}(s, 1) = nanmean(correcttrls);
        
        fname = sprintf('%s/P03/behav/%s_4decode_rtco.mat', mdir, subj{s})
        beh = load(fname); % behavioral data
        clear tind vind
        tind = ismember(beh.sess, sess_id{c+1});
        sampmus = mean(beh.samples(tind, 1:6), 2);
        if c == 1
            interm = beh.interm_cues(tind);
        elseif c == 2
            interm = beh.binary_choices(tind);
        end
        
        % valid trials based on interm button press
        ind1 = beh.binary_choices(tind)~=99; % missed interm button press (includes both intervals)
        % valid trials based on estimation
        ind2 = ~isnan(beh.estimations(tind)); % missed estimation
        
        if fgmu > 0
            % valid trials2: based on generative mean:
            zgmu = -fgmu:fgmu;
            zind = ismember(genmus, zgmu);
        elseif fgmu == 0 % use all trials
            zind = ~isnan(genmus(tind));
        end
        pztrls = zind & ind1 & ind2;
        % for the correct trials
        correcttrls = interm(pztrls).*sampmus(pztrls)>0;
        saccu{c}(s, 1) = nanmean(correcttrls);
    end
    
end

%% plot bar graphs and compare
hw = 5;
vw = 5;

vis = figure;
vis.Units = 'centimeter';
vis.Position = [1 1 hw vw];

choice_subs = setdiff(1:30, [1 3 5]); choice_subs  = choice_subs';

% subplot(1, 2, 1)
% bh1 = bar([0.9 2.1], 100.*[nanmean(gaccu{1}(choice_subs)) nanmean(gaccu{2}(choice_subs))]);
% bh1.FaceColor = 'flat';
% bh1.CData = cols;
% bh1.BarWidth = 0.5;
hold all;
data = [gaccu{1}(choice_subs) gaccu{2}(choice_subs)];
ri = rand(length(data(:, 1)), 1); ri = 0.3*(ri-0.5);
ph1 = scatter(1+ri, 100*data(:, 1), 25,repmat([0 0 0], length(data(:, 1)), 1), 'filled');
ph2 = scatter(2+ri, 100*data(:, 2), 25, repmat(0.6.*[1 1 1], length(data(:, 1)), 1), 'filled');
ph1.MarkerFaceAlpha = 0.3;
ph2.MarkerFaceAlpha = 0.3;

plot([1-0.3 1+0.3], [mean(100*data(:, 1)) mean(100*data(:, 1))], '-k', 'LineWidth', 1)
plot([2-0.3 2+0.3], [mean(100*data(:, 2)) mean(100*data(:, 2))], '-k', 'color', [0.6 0.6 0.6], 'LineWidth', 1)
upl = 1.05*100*max(data(:));
plot([1 2], upl*[1 1], '-k', 'LineWidth', 1); plot([1 1], [upl-1 upl], '-k'); plot([2 2], [upl-1 upl], '-k')
[p,h,stats] = signrank(gaccu{1}(choice_subs), gaccu{2}(choice_subs));
if p<0.001
    text(1, upl+2, sprintf('p<0.001'), 'FontSize', 8, 'FontName', 'Helvetica')
else
    text(1, upl+2, sprintf('p=%1.3f',p), 'FontSize', 8, 'FontName', 'Helvetica')
end
% hold on; plot([0.9 0.9], 100.*[nanmean(gaccu{1}(choice_subs))+std(gaccu{1}(choice_subs),[],1)/sqrt(length(choice_subs)) nanmean(gaccu{1}(choice_subs))-std(gaccu{1}(choice_subs)/sqrt(length(choice_subs)),[],1)], 'k', 'LineWidth', 2)
% plot([2.1 2.1], 100.*[nanmean(gaccu{2}(choice_subs))+std(gaccu{2}(choice_subs),[],1)/sqrt(length(choice_subs)) nanmean(gaccu{2}(choice_subs))-std(gaccu{2}(choice_subs),[],1)/sqrt(length(choice_subs))], 'r', 'LineWidth', 2)
title('Near-0 trials')
ylim([45 78]), xlim([0 3]), box off
set(gca, 'xtick', [0.9 2.1], 'xticklabels', condt, 'tickdir', 'out', 'ytick', [50 60 70])
ylabel('Accuracy (%)')
% subplot(1, 2, 2)
% bh2 = bar([1 2], nanmean([saccu{1} saccu{2}], 1));
% bh2.FaceColor = 'flat';
% bh2.CData = cols;
% hold on; plot([1 1], [nanmean(saccu{1})+std(saccu{1},[],1) nanmean(saccu{1})-std(saccu{1},[],1)], 'k', 'LineWidth', 2)
% plot([2 2], [nanmean(saccu{2})+std(saccu{2},[],1) nanmean(saccu{2})-std(saccu{2},[],1)], 'r', 'LineWidth', 2)
% title('sample(1-6) mean')
% ylim([0 0.8]), box off
% set(gca, 'xticklabels', condt, 'tickdir', 'out')

sname = sprintf('%s/perizo_accuracy_dp.eps', fig_dir);
saveas(vis, sname)
exportgraphics(vis,sname,'ContentType','vector')

%% do stats
pval = [];
nperm = 10000;
tmpe = gaccu{1}(choice_subs)-gaccu{2}(choice_subs);
mu_orig = mean(tmpe);
tmpce = cat(1, gaccu{1}(choice_subs), gaccu{2}(choice_subs));
mu_boot = zeros(size(tmpe, 2), nperm);
for b=1:nperm
    % randomly switch labels between Cue and Choice
    rs = cat(1, ones(size(tmpe,1), 1), 2*ones(size(tmpe,1), 1));
    rs = logical(rs(randperm(length(rs)))-1);
    tmp1v = tmpce(rs, :);
    tmp2v = tmpce(~rs, :);
    tmpv = tmp1v-tmp2v;
    mu_boot(:,b) = mean(tmpv);
end
pval(1) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
% tmpe = saccu{1}-saccu{2};
% mu_orig = mean(tmpe);
% tmpce = cat(1, saccu{1}, saccu{2});
% mu_boot = zeros(size(tmpe, 2), nperm);
% boot_ind = zeros(nperm, size(tmpce, 1));
% for b=1:nperm
% %     % randomly switch labels between Cue and Choice
% %     rs = cat(1, ones(size(tmpe,1), 1), 2*ones(size(tmpe,1), 1));
% %     rs = logical(rs(randperm(length(rs)))-1);
% %     tmp1v = tmpce(rs, :);
% %     tmp2v = tmpce(~rs, :);
% %     tmpv = tmp1v-tmp2v;
% %     mu_boot(:,b) = mean(tmpv);
%     rs = randperm(size(tmpce, 1));
%     boot_ind(b, :) = rs;
% end
% boot_ind = unique(boot_ind, 'rows');
% boot_ind = boot_ind';
% for b = 1:nperm
%     tmp1v = tmpce(boot_ind(1:30, b), 1);
%     tmp2v = tmpce(boot_ind(31:60, b), 1);
%     mu_boot(:, b) = mean(tmp1v-tmp2v);
% end
% pval(2) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;

% Wilcoxon signed rank test
[p,h,stats] = signrank(gaccu{1}(choice_subs), gaccu{2}(choice_subs))
[p,h,stats] = signrank(saccu{1}(choice_subs), saccu{2}(choice_subs))

