function estimation_histograms()
clc;close all;
subjects = [1:11, 13 15:18 20 22:34];
all_subj.sj = [];
all_subj.estim = [];
all_subj.norm_estim = [];
all_subj.estim_rt = [];
all_subj.gen_means = [];
all_subj.sample_means = [];
all_subj.interm_choice = [];
all_subj.interm_cue = [];
all_subj.interm_rt = [];
all_subj.sess = [];
all_subj.interm_resp_accu = [];
all_subj.outlier = [];
for sj = subjects
    load(sprintf('behavior_trials_rtco/S%02d_4decode_rtco.mat', sj));
    all_subj.sj = [all_subj.sj; sj*ones(length(estimations), 1)];
    all_subj.estim = [all_subj.estim; estimations];
    all_subj.norm_estim = [all_subj.norm_estim; nanzscore(estimations)];
    all_subj.estim_rt = [all_subj.estim_rt; estim_rt];
    all_subj.gen_means = [all_subj.gen_means; gen_means];
    all_subj.sample_means = [all_subj.sample_means; sample_means];
    all_subj.interm_choice = [all_subj.interm_choice; binary_choices];
    all_subj.interm_cue = [all_subj.interm_cue; interm_cues];
    all_subj.interm_rt = [all_subj.interm_rt; interm_rt];
    all_subj.sess = [all_subj.sess; sess];
    all_subj.interm_resp_accu = [all_subj.interm_resp_accu; interm_resp_accu];
    q = prctile(estim_rt, [25, 75]);
    int_quart = 1.5*iqr(estim_rt);
    outliers = zeros(length(estimations), 1);
    outliers(find(estim_rt >= q(2)+int_quart | estim_rt <= q(1)-int_quart)) = 1;
    all_subj.outlier = [all_subj.outlier; outliers];
    clearvars -except all_subj subjects load_data
end
all_subj = struct2table(all_subj);

n_subj = length(unique(all_subj.sj));
subjects = unique(all_subj.sj)';
lw = 1;
choice_col = [0, 0, 0];
cue_col = [0.8, 0.8, 0.8];
% plot cue and choice histograms separately in one figure and pooled across
% conditions in the other
% only show stairs not bars
% tiledlayout(6,5,'TileSpacing','tight', 'Padding', 'compact');
clf;
for i = 1:n_subj
    subplot(6,5,i);hold on;
    bw = range(all_subj.estim(all_subj.sj == subjects(i)))/20;
    histogram(all_subj.estim(all_subj.sj == subjects(i) & isnan(all_subj.interm_cue)), 'BinWidth', bw, 'DisplayStyle','stairs', 'EdgeColor', choice_col, 'EdgeAlpha', 0.6, 'LineWidth', lw);
    histogram(all_subj.estim(all_subj.sj == subjects(i) & ~isnan(all_subj.interm_cue)), 'BinWidth', bw, 'DisplayStyle','stairs', 'EdgeColor', cue_col, 'EdgeAlpha', 0.8, 'LineWidth', lw);
    axis tight;box("off");
    xlim = get(gca, 'XLim');
    set(gca, 'XLim', [-max(abs(xlim)), max(abs(xlim))]);
    [~, pCh] = HartigansDipSignifTest(all_subj.estim(all_subj.sj == subjects(i) & isnan(all_subj.interm_cue) & ~isnan(all_subj.estim)), 1000, 0);
    [~, pC] = HartigansDipSignifTest(all_subj.estim(all_subj.sj == subjects(i) & ~isnan(all_subj.interm_cue) & ~isnan(all_subj.estim)), 1000, 0);
    if pCh < 0.001 && pC < 0.001
        title({sprintf('S %d', subjects(i)), 'p < 0.001, p < 0.001'}, "FontWeight","bold",'Fontsize',7);
    elseif pC < 0.001 && pCh >= 0.001
        title({sprintf('S %d', subjects(i)), sprintf('p = %.3f, p < 0.001', pCh)}, "FontWeight","bold",'Fontsize',7);
    elseif pCh < 0.001 && pC >= 0.001
        title({sprintf('S %d', subjects(i)), sprintf('p < 0.001, p = %.3f', pC)}, "FontWeight","bold",'Fontsize',7);
    else
        title({sprintf('S %d', subjects(i)), sprintf('p = %.3f, p =%.3f', pCh, pC)}, "FontWeight","bold",'Fontsize',7);
    end
end

suplabel('Estimations', 'x');
suplabel('Bin Count', 'y');
suplabel('Estimation histograms per subject', 't')
h=gcf;
set(h, 'PaperType', 'A4');
set(h,'PaperOrientation','portrait');
exportgraphics(gcf,'FigS1.pdf','ContentType','vector');

clf;hold on;lw=2;
bw = range(all_subj.norm_estim)/75;
histogram(all_subj.norm_estim(isnan(all_subj.interm_cue)), 'BinWidth', bw, 'DisplayStyle','stairs', 'EdgeColor', choice_col, 'EdgeAlpha', 0.6, 'LineWidth', lw);
histogram(all_subj.norm_estim(~isnan(all_subj.interm_cue)), 'BinWidth', bw, 'DisplayStyle','stairs', 'EdgeColor', cue_col, 'EdgeAlpha', 0.8, 'LineWidth', lw);
axis tight;box("off");
set(gca, 'XLim', [-4, 4], 'YLim', [0, 1600]);
legend('Choice', 'Cue', 'Location','best'); legend('Box','off');
[~, pCh] = HartigansDipSignifTest(all_subj.norm_estim(isnan(all_subj.interm_cue) & ~isnan(all_subj.norm_estim)), 1000, 0);
[~, pC] = HartigansDipSignifTest(all_subj.norm_estim(~isnan(all_subj.interm_cue) & ~isnan(all_subj.norm_estim)), 1000, 0);
if pCh < 0.001 && pC < 0.001
        title({'All subjects', 'p < 0.001, p < 0.001'}, "FontWeight","bold",'Fontsize',7);
    elseif pC < 0.001 && pCh >= 0.001
        title({'All subjects', sprintf('p = %.3f, p < 0.001', pCh)}, "FontWeight","bold",'Fontsize',7);
    elseif pCh < 0.001 && pC >= 0.001
        title({'All subjects', sprintf('p < 0.001, p = %.3f', pC)}, "FontWeight","bold",'Fontsize',7);
    else
        title({'All subjects', sprintf('p = %.3f, p =%.3f', pCh, pC)}, "FontWeight","bold",'Fontsize',7);
    end
xlabel('Estimations normalsied within subjects');
ylabel('Bin Count');

h=gcf;
set(h, 'PaperType', 'A4');
set(h,'PaperOrientation','portrait');
exportgraphics(gcf,'Fig_1g.pdf','ContentType','vector');
