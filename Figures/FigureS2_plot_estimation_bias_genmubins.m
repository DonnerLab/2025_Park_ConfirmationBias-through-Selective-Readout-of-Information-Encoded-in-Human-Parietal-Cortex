% relevant data: raw behavioral data: Sxx_4decode_rtco.mat

%% Parameters

sess_id = {[1 2 3 4], [1 2], [3 4]};
condt = {'Cue','Choice'};

colslr = {};
colslr{1} = [55,126,210; 45,116,200]./255; % left, light, dark
colslr{2} = [228,26,28; 218,16,18]./255; % right, light, dark

%% get bias
ugmu = -14:14;
gconfb = cell(2,1);

load('subjects.mat')
subj([1 3 5]) = []; % delete subjects with coarse generative means

for c = 2:3
    for s = 1:length(subj)
        fname = sprintf('%s_4decode_rtco.mat', subj{s});
        beh = load(fname); % behavioral data
        % get session indices
        clear tind
        tind = ismember(beh.sess, sess_id{c});
        vind = logical(beh.ref_angles(tind));
        % get interm data
        if c == 2
            interm = beh.interm_cues(tind);
        elseif c == 3
            interm = beh.binary_choices(tind);
        end
        gmu = beh.gen_means(tind);
        eaccu = beh.estimations(tind);
        mind = ~isnan(eaccu);
        iind = interm~=99;
        eaccu = eaccu(mind&iind);
        zeaccu = zscore(eaccu);
        gmu = gmu(mind&iind);
        interm = interm(mind&iind);
        for g = 1:length(ugmu)
            if sum(gmu==ugmu(g))~=0
                gconfb{c-1}{g}{1}{s} = zeaccu(gmu==ugmu(g)&interm==-1); % interm: left (-1)
                gconfb{c-1}{g}{2}{s} = zeaccu(gmu==ugmu(g)&interm==1); % interm: left (-1)
            else
                gconfb{c-1}{g}{1}{s} = nan;
                gconfb{c-1}{g}{2}{s} = nan;
            end
        end
    end
end

%% save bias
% average error for each gen mean.
mconfb = cell(2, 1);
nconfb = cell(2, 1);
for c = 1:2
    tmp1 = nan(length(ugmu), length(subj), 2);
    tmp2 = nan(length(ugmu), length(subj), 2);
    for g = 1:length(ugmu)
        tmp1(g, :, 1) = cell2mat(cellfun(@(x) nanmean(x), gconfb{c}{g}{1}, 'UniformOutput', false));
        tmp1(g, :, 2) = cell2mat(cellfun(@(x) nanmean(x), gconfb{c}{g}{2}, 'UniformOutput', false));
        tmp2(g, :, 1) = cell2mat(cellfun(@(x) length(~isnan(x)), gconfb{c}{g}{1}, 'UniformOutput', false));
        tmp2(g, :, 2) = cell2mat(cellfun(@(x) length(~isnan(x)), gconfb{c}{g}{2}, 'UniformOutput', false));
    end
    mconfb{c} = tmp1;
    nconfb{c} = tmp2;
end

%% Plot estimation and also bins bars

gbins = cell(3, 1);
gbins{1} = [11:19]; % near-zero trials
gbins{2} = [6:10 20:24]; % medium trials
gbins{3} = [1:5 25:29]; % easy trials


%%%%%% GEN MU %%%%%%%%
vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 15 7];
annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('N = %d', length(subj)), 'LineStyle', 'none', 'fontsize', 10)

xlm = [ugmu(1)-0.5 ugmu(end)+0.5];
ylm = [-4 3];

ms = 2;
for c = 1:2  
    ph1 = []; ph2=[];
    subplot(1, 2, c); 
    plot(xlm, [0 0], ':', 'color', 'k', 'LineWidth', 0.1, 'HandleVisibility', 'off'); hold on;
    plot([0 0], ylm, ':', 'color', 'k', 'LineWidth', 0.1, 'HandleVisibility', 'off'); hold on;
    % left
    ph1 = plot(ugmu-0.2, mconfb{c}(:, :, 1), 'o', 'LineStyle', 'none', ... 
        'markersize', ms, 'markerfacecolor', colslr{1}(1,:), 'markeredgecolor', 'none');
    % right
    ph2 = plot(ugmu+0.2, mconfb{c}(:, :, 2), 'o', 'LineStyle', 'none', ... 
        'markersize', ms, 'markerfacecolor', colslr{2}(1, :), 'markeredgecolor', 'none');
    hold all;
    plot(ugmu(gbins{1}), ylm(1)*(ones(1, length(gbins{1}))), 'sk','markerfacecolor','k','markeredgecolor','none', 'HandleVisibility', 'off')
    text([ugmu(gbins{1}(1)) ugmu(gbins{1}(end))]-0.8, (ylm(1)+0.3)*[1 1], [string(ugmu(gbins{1}(1))) string(ugmu(gbins{1}(end)))], 'FontSize',8,'FontName','Helvetica')
    plot(ugmu(gbins{2}), ylm(1)*(ones(1, length(gbins{2}))), 'sk','markerfacecolor',0.5*[1 1 1],'markeredgecolor','none', 'HandleVisibility', 'off')
    text([ugmu(gbins{2}(1)) ugmu(gbins{2}(end))]-0.8, (ylm(1)+0.3)*[1 1], [string(ugmu(gbins{2}(1))) string(ugmu(gbins{2}(end)))], 'FontSize',8,'FontName','Helvetica')
    plot(ugmu(gbins{3}), ylm(1)*(ones(1, length(gbins{3}))), 'sk','markerfacecolor',0.8*[1 1 1],'markeredgecolor','none', 'HandleVisibility', 'off')
    text([ugmu(gbins{3}(1)) ugmu(gbins{3}(end))]-0.8, (ylm(1)+0.3)*[1 1], [string(ugmu(gbins{3}(1))) string(ugmu(gbins{3}(end)))], 'FontSize',8,'FontName','Helvetica')
    hold on;
    
    xlim(xlm);
    ylim(ylm);

    th = title(condt{c});
    if c == 1
        legend([ph1(1);ph2(1)], 'Left Cue', 'Right Cue', 'Location', 'NorthWest')
        legend boxoff
    elseif c == 2
        legend([ph1(1);ph2(1)], 'Left Choice', 'Right Choice', 'Location', 'NorthWest')
        legend boxoff
    end
    axis('square')
    axo = gca;
    axo.Position = [axo.Position(1) axo.Position(2)+0.08 axo.Position(3) axo.Position(4)-0.1];
    [ay, ax] = offsetaxis(axo, 'y', 0.02, 'x', 0.02);
    set(ax, 'tickdir', 'out', 'FontName', 'Helvetica')
    set(ay, 'tickdir', 'out', 'FontName', 'Helvetica')

    ylabel(ay, 'Estimation (zscore)')
    xlabel(ax, 'Generative mean (deg)')
    box off    
end

sname = sprintf('%s/Estim_BIAS_GENMU_zscore_wofit_bins', fig_dir);
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, sname, 'png')



%% plot estimation by genmu bins

ugmu = -14:14;
gbins = cell(3, 1);
gbins{1} = [11:19]; % near-zero trials
gbins{2} = [6:10 20:24]; % medium trials
gbins{3} = [1:5 25:29]; % easy trials

xhic = cell(2,1);
for c = 1:2
    for s = 1:length(subj)        
        for gb = 1:length(gbins)
            xhic{c}(s, gb) = nanmean(mconfb{c}(gbins{gb}, s, 2)-mconfb{c}(gbins{gb}, s, 1));
        end
    end
end

%% plot bias

cols = [0.7.*[1 1 1]; 0 0 0];
vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 5 5];
annotation(vis, 'textbox', [0.02 0.9 0.1 0.1], 'string', sprintf('N = %d', length(subj)), 'LineStyle', 'none', 'fontsize', 9)
sh1 = [];
yline(0, ':k'); hold on;
for c = 1:2
    mu =  nanmean(xhic{c},1); serr = nanstd(xhic{c},[],1)/sqrt(length(subj));
    plot([1:3; 1:3]+0.05*(c-1), [mu-serr; mu+serr], '-', 'color', cols(c,:), 'LineWidth', 2); hold on;    
    sh1(c) = plot((1:3)+0.05*(c-1), mu, '-', 'color', cols(c,:), 'LineWidth', 2);
    hold on;
end
lh = legend(sh1, 'Cue','Choice', 'FontSize', 7); legend boxoff
lh.Position = lh.Position+[-0.1 0.08 0 0];
xlim([0.5 3.5]); ylim([-0.1 1.1]); 
ylabel('\DeltaEstimation (R-L)', 'FontName', 'Helvetica')
set(gca, 'xtick',1:3,'xticklabels', {'bin1','bin2','bin3'}, 'ytick', [0 1],'tickdir', 'out', 'FontName', 'Helvetica')
axis('square')
box off
sname = sprintf('Estim_BIAS_GENMU_zscore_bins_27subj');
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%        'BackgroundColor','none')



%% Do anova and compare between Cue and Choice
y = cat(1, xhic{1}(:), xhic{2}(:));
N = length(subj);
g1 = [ones(N*3, 1); 2*ones(N*3, 1)];
g2 = [ones(N, 1); 2*ones(N, 1); 3*ones(N, 1); ones(N, 1); 2*ones(N, 1); 3*ones(N, 1)];
[~, T] = anovan(y, {g1, g2}, 'model', 'interaction', 'varnames', {'g1','g2'});

cfg = [];
cfg.ttest = 0;
cfg.nperm = 5000;
cfg.perm = 'paired';
cfg.tail = 'both';
p = [];
for gb = 1:3
    p(gb) = perm_test_behav(xhic{1}(:, gb), xhic{2}(:, gb), cfg);
end


%% plot near-zero L/R bar graphs

ugmu = -14:14;
gbins = cell(3, 1);
gbins{1} = [11:19]; % near-zero trials
gbins{2} = [6:10 20:24]; % medium trials
gbins{3} = [1:5 25:29]; % easy trials

xhic = cell(2,1);
for c = 1:2
    for s = 1:length(subj)
        for gb = 1:length(gbins)
            for lr = 1:2
                xhic{c}(lr, s, gb) = nanmean(mconfb{c}(gbins{gb}, s, lr));
            end
        end
    end
end

xpts = [1 2 4 5];
coltmp = [colslr{1}(1,:); colslr{2}(1,:); colslr{1}(1,:); colslr{2}(1,:)];

gb = 1;
N = length(subj);
vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 11 5];
annotation(vis, 'textbox', [0.03 0.92 0.1 0.1], 'string', sprintf('N = %d', length(subj)), 'LineStyle', 'none', 'fontsize', 10)
subplot(1, 2, 1)
yline(0, ':k'), hold on;
data = [xhic{1}(1, :, gb)' xhic{1}(2, :, gb)' xhic{2}(1, :, gb)' xhic{2}(2, :, gb)'];
for xp = 1:length(xpts)
    sh = scatter(ones(N, 1)*xpts(xp), data(:, xp), 16, coltmp(xp,:), 'filled');
    sh.MarkerFaceAlpha = 0.3;
    hold on;
    plot([xpts(xp)-0.2 xpts(xp)+0.2], mean(data(:, xp)).*[1 1], 'LineStyle', '-', 'color', coltmp(xp,:), 'LineWidth', 3);
    hold on;
end
axis('square')
xlim([0 6]); ylim([-1.5 1.5]); ylabel(sprintf('near-zero\nGenerative \\mu'))
set(gca, 'xtick',xpts,'xticklabels', {'L','R','L','R'}, 'ytick', [-1 0 1],'tickdir', 'out')

cfg = [];
cfg.ttest = 0;
cfg.nperm = 5000;
cfg.perm = 'paired';
cfg.tail = 'both';

xpts = [1 2;4 5];
for c = 1:2
    data = [xhic{c}(1, :, gb)' xhic{c}(2, :, gb)'];
    plot([xpts(c,:); xpts(c,:)], max(data(:))+[1.2 1.25; 1.25 1.2]-1, '-k'), hold on;
    plot(xpts(c,:), max(data(:))+[1.25 1.25]-1, '-k')
    
    p = perm_test_behav(data(:, 1), data(:, 2), cfg);
    if p < 0.001
        text(mean(xpts(c,:))-1.2,max(data(:))+0.6,sprintf('p < 10^{-3}'))        
    else
        text(mean(xpts(c,:))-0.8,max(data(:))+0.6,sprintf('p = %1.3f', p))
    end
end

coltmp = [0.6 0.6 0.6; 0 0 0];
xpts = [1 2];
subplot(1, 2, 2)
yline(0, ':k'), hold on;
data = [xhic{1}(2, :, gb)' - xhic{1}(1, :, gb)' xhic{2}(2, :, gb)' - xhic{2}(1, :, gb)'];
for xp = 1:length(xpts)
    sh = scatter(ones(N, 1)*xpts(xp), data(:, xp), 16, coltmp(xp,:), 'filled'); hold on;
    plot([xpts(xp)-0.1 xpts(xp)+0.1], mean(data(:, xp)).*[1 1], 'LineStyle', '-', 'color', coltmp(xp,:), 'LineWidth', 3);
    sh.MarkerFaceAlpha = 0.3;
    hold on;
end
axis('square')
plot([xpts; xpts], max(data(:)).*[1.1 1.15; 1.15 1.1], '-k'), hold on;
plot(xpts, max(data(:)).*[1.15 1.15], '-k')
p = perm_test_behav(data(:, 1), data(:, 2), cfg);
text(mean(xpts)-0.5,max(data(:))+0.7, sprintf('p = %1.3f', p))
xlim([0.5 2.5]); ylim([-0.5 3]); ylabel(sprintf('\\DeltaEstimation\n(Right-Left)'))
set(gca, 'xtick',xpts,'xticklabels', {'Cue','Choice'}, 'ytick', [-1 0 1 2],'tickdir', 'out')
sname = sprintf('Estim_BIAS_GENMU_zscore_singlesub_xnz_27subj')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, [sname '.png'])
