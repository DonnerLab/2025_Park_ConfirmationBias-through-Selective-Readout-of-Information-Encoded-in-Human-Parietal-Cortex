
%% Parameters
load('subjects')

conditions = {'attn', 'choice'};
cic = {'con','incon'};

sess_id = {[1 2 3 4], [1 2], [3 4]};
condt = {'Cue','Choice'};

% colormap gray;
qtl_vals = [0.2, 0.6, 1.0];

cols = {};
cols{1} = [10 184 0]./255; % consistent
cols{2} = [255 184 10]./255; % inconsistent

colslr = {};
colslr{1} = [55,126,210; 45,116,200]./255; % left, light, dark
colslr{2} = [228,26,28; 218,16,18]./255; % right, light, dark

ugmu = union(-15:5:15, -14:14);
edin1 = 1:5:length(ugmu);
edin2 = 1:length(ugmu);
sub1 = [1 3 5];
subjects = setdiff(1:length(subjects), sub1);

%% get confirmation bias
gconfb = cell(2,1);
ugmu = union(-15:5:15, -14:14); sigma = 16; % 8 brings around 77%
% average error for each gen mean.
mconfb = cell(2, 1);
nconfb = cell(2, 1);

for c = 1:2 % only here c = 2 is attn, c = 3 is choice   
    for s = 1:length(subjects)
        jobin = {subjects{s},c+1,sigma};
        beh_crr = simul_confbias_estimations_CC(jobin, mdir);
        gmu = beh_crr.gen_means;
        eaccu = beh_crr.estimations;
        %     mind = ~isnan(eaccu);
        %     iind = interm~=99;
        %     eaccu = eaccu(mind&iind);
        zeaccu = zscore(eaccu);
        %     gmu = gmu(mind&iind);
        %     interm = interm(mind&iind);
        interm = beh_crr.interm;
        for g = 1:length(ugmu)
            if sum(gmu==ugmu(g))~=0
                gconfb{c}{g}{1}{s} = zeaccu(gmu==ugmu(g)&interm==-1); % interm: left (-1)
                gconfb{c}{g}{2}{s} = zeaccu(gmu==ugmu(g)&interm==1); % interm: left (-1)
            else
                gconfb{c}{g}{1}{s} = nan;
                gconfb{c}{g}{2}{s} = nan;
            end
        end
    end
    % save bias 
    tmp1 = nan(length(ugmu), length(subjects), 2);
    tmp2 = nan(length(ugmu), length(subjects), 2);
    for g = 1:length(ugmu)
        tmp1(g, :, 1) = cell2mat(cellfun(@(x) nanmean(x), gconfb{c}{g}{1}, 'UniformOutput', false)); % left
        tmp1(g, :, 2) = cell2mat(cellfun(@(x) nanmean(x), gconfb{c}{g}{2}, 'UniformOutput', false)); % right
        tmp2(g, :, 1) = cell2mat(cellfun(@(x) length(~isnan(x)), gconfb{c}{g}{1}, 'UniformOutput', false));
        tmp2(g, :, 2) = cell2mat(cellfun(@(x) length(~isnan(x)), gconfb{c}{g}{2}, 'UniformOutput', false));
    end
    mconfb{c} = tmp1;
    nconfb{c} = tmp2;
    save(fullfile(dat_dir, sprintf('/behav/confbias_simul_Estim_sigma%d_gmu_zscore.mat', sigma)), 'mconfb', 'nconfb')
end


%% plot estimation by genmu bins

subjects = get_sinfo(setdiff(1:34, [1 3 5 12 14 19 21]), 2);
ugmu = union(-15:5:15, -14:14);
gbins = cell(3, 1);
gbins{1} = [12:15 16 17:20]; % near-zero trials
gbins{2} = [7:11 21:25]; % medium trials
gbins{3} = [1:6 26:31]; % easy trials
sigma = 0;
load(fullfile(dat_dir, sprintf('/behav/confbias_simul_Estim_sigma%d_gmu_zscore.mat', sigma)), 'mconfb')

xhic = cell(2,1);
for c = 1:2
    for s = 1:length(subjects)        
        for gb = 1:length(gbins)
            xhic{c}(s, gb) = nanmean(mconfb{c}(gbins{gb}, s, 2)-mconfb{c}(gbins{gb}, s, 1));
        end
    end
end

%% PLOT bias
cols = [0.7.*[1 1 1]; 0 0 0];
vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 5 5];
% annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('N = %d', length(subjects)), 'LineStyle', 'none', 'fontsize', 12)
sh1 = [];
yline(0, ':k'), hold on;
for c = 1:2
    mu =  mean(xhic{c},1); serr = std(xhic{c},[],1)/sqrt(length(subjects));
    plot([1:3; 1:3]+0.05*(c-1), [mu-serr; mu+serr], '-', 'color', cols(c,:), 'LineWidth', 2); hold on;    
    sh1(c) = plot((1:3)+0.05*(c-1), mu, '-', 'color', cols(c,:), 'LineWidth', 2);
    hold on;
end
lh = legend(sh1, 'Cue','Choice', 'FontSize', 7); legend boxoff
% lh = legend(sh1(c), condt{c}, 'FontSize', 7); legend boxoff
lh.Position = lh.Position+[-0.1 0.08 0 0];
xlim([0.5 3.5]); ylim([-0.1 1.1]); 
ylabel(sprintf('\\DeltaSimulated Estimation\n(R-L)'), 'FontName', 'Helvetica', 'FontSize', 9)
set(gca, 'xtick',1:3,'xticklabels', {'bin1','bin2','bin3'}, 'ytick', [0 1],'tickdir', 'out', 'FontName', 'Helvetica')
axis('square')
box off
sname = sprintf('%s/Simul_Estim_sigma%d_BIAS_GENMU_zscore_bins_27subj', fig_dir, sigma);
exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
       'BackgroundColor','none')
% exportgraphics(vis, [sname '.pdf'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, [sname '.png'])

