% relevant data: raw behavioral data: Sxx_4decode_rtco.mat

%% Parameters

sess_id = {[1 2 3 4], [1 2], [3 4]};
condt = {'Cue','Choice'};

colslr = {};
colslr{1} = [55,126,210; 45,116,200]./255; % left, light, dark
colslr{2} = [228,26,28; 218,16,18]./255; % right, light, dark

%% get bias
load('subjects.mat')
edges = -27.5:1:27.5;
nsmu = 6; % sample mean based on first 6 samples
gconfb = cell(2,1);
gbinc = cell(2,1);
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
%         smu = mean(beh.samples(tind, :), 2);
        smu = mean(beh.samples(tind, 1:6), 2); % sample mean based on first 6 samples
        eaccu = beh.estimations(tind);
        mind = ~isnan(eaccu);
        iind = interm~=99;
        eaccu = eaccu(mind&iind);
        zeaccu = zscore(eaccu);
        smu = smu(mind&iind);
        interm = interm(mind&iind);
        [N,edges1,ugsmu] = histcounts(smu, edges);
        if find(ugsmu==0)
            fprintf('bin again.')
        end
        for g = 1:length(edges)-1
            if sum(ugsmu==g)~=0
                gconfb{c-1}{g}{1}{s} = zeaccu(ugsmu==g&interm==-1); % interm: left (-1)
                gconfb{c-1}{g}{2}{s} = zeaccu(ugsmu==g&interm==1); 
                gbinc{c-1}{g}{1}{s} = (ugsmu==g)&(interm==-1);                
                gbinc{c-1}{g}{2}{s} = (ugsmu==g)&(interm==1);
                gbinc{c-1}{g}{3}{s} = sum(ugsmu==g);
            else
                fprintf('zero\n')
                gconfb{c-1}{g}{1}{s} = nan;
                gconfb{c-1}{g}{2}{s} = nan;
                gbinc{c-1}{g}{1}{s} = nan;
                gbinc{c-1}{g}{2}{s} = nan;
                gbinc{c-1}{g}{3}{s} = nan;
            end
        end
    end
end

%% plot confimration bias for sample mean
% average error for each sample mean bin.
mconfb = cell(2, 1);
nconfb = cell(2, 1);
for c = 1:2
    tmp1 = nan(length(N), length(subj), 2);
    tmp2 = nan(length(N), length(subj), 2);
    for g = 1:length(N)
        tmp1(g, :, 1) = cell2mat(cellfun(@(x) nanmean(x), gconfb{c}{g}{1}, 'UniformOutput', false));
        tmp1(g, :, 2) = cell2mat(cellfun(@(x) nanmean(x), gconfb{c}{g}{2}, 'UniformOutput', false));
        tmp2(g, :, 1) = cell2mat(cellfun(@(x) length(~isnan(x)), gconfb{c}{g}{1}, 'UniformOutput', false));
        tmp2(g, :, 2) = cell2mat(cellfun(@(x) length(~isnan(x)), gconfb{c}{g}{2}, 'UniformOutput', false));
    end
    mconfb{c} = tmp1;
    nconfb{c} = tmp2;
end

%% Plot estimation and also bins bars

edges = -27.5:1:27.5;
tedges = edges>=-22.5 & edges<=22.5;
tedges(end) = [];

sbins = cell(3, 1);
sbins{1} = [24:27 28 29:32]; % near-zero trials
sbins{2} = [18:24 32:38]; % medium trials
sbins{3} = [6:18 39:51]; % easy trials

%%%%%% SAMPLE MU %%%%%%%%
vis = figure;
vis.Units = 'centimeter';
vis.Position = [2 2 15 7];
annotation(vis, 'textbox', [0.03 0.92 0.1 0.1], 'string', sprintf('N = %d', length(subj)), 'LineStyle', 'none', 'fontsize', 10)

xedges = edges(tedges);
xlm = [xedges(1)-0.5 xedges(end)+0.5];
ylm = [-4 3];

ms = 2;
for c = 1:2  
    ph1 = []; ph2=[];
    subplot(1, 2, c); 
    plot(xlm, [0 0], ':', 'color', 'k', 'LineWidth', 0.1, 'HandleVisibility', 'off'); hold on;
    plot([0 0], ylm, ':', 'color', 'k', 'LineWidth', 0.1, 'HandleVisibility', 'off'); hold on;
    % left
    ph1 = plot(edges(tedges)-0.2, mconfb{c}(tedges, :, 1), 'o', 'LineStyle', 'none', ... 
        'markersize', ms, 'markerfacecolor', colslr{1}(1,:), 'markeredgecolor', 'none');
    hold on;
    % right
    ph2 = plot(edges(tedges)+0.2, mconfb{c}(tedges, :, 2), 'o', 'LineStyle', 'none', ... 
        'markersize', ms, 'markerfacecolor', colslr{2}(1, :), 'markeredgecolor', 'none');
    hold all;
    plot(edges(sbins{1}), ylm(1)*(ones(1, length(sbins{1}))), 'sk','markerfacecolor','k','markeredgecolor','none', 'HandleVisibility', 'off')
    text([edges(sbins{1}(1)) edges(sbins{1}(end))]-0.8, (ylm(1)+0.3)*[1 1], [string(edges(sbins{1}(1))) string(edges(sbins{1}(end)))], 'FontSize',8,'FontName','Helvetica')
    plot(edges(sbins{2}), ylm(1)*(ones(1, length(sbins{2}))), 'sk','markerfacecolor',0.5*[1 1 1],'markeredgecolor','none', 'HandleVisibility', 'off')
    text([edges(sbins{2}(1)) edges(sbins{2}(end))]-0.8, (ylm(1)+0.3)*[1 1], [string(edges(sbins{2}(1))) string(edges(sbins{2}(end)))], 'FontSize',8,'FontName','Helvetica')
    plot(edges(sbins{3}), ylm(1)*(ones(1, length(sbins{3}))), 'sk','markerfacecolor',0.8*[1 1 1],'markeredgecolor','none', 'HandleVisibility', 'off')
    text([edges(sbins{3}(1)) edges(sbins{3}(end))]-0.8, (ylm(1)+0.3)*[1 1], [string(edges(sbins{3}(1))) string(edges(sbins{3}(end)))], 'FontSize',8,'FontName','Helvetica')
    hold on;

    th = title(condt{c});
%     if c == 1
%         legend([ph1(1);ph2(1)], 'Left Cue', 'Right Cue', 'Location', 'NorthWest')
%         legend boxoff
%     elseif c == 2
%         legend([ph1(1);ph2(1)], 'Left Choice', 'Right Choice', 'Location', 'NorthWest')
%         legend boxoff
%     end
    axis('square')
    axo = gca;
    axo.Position = [axo.Position(1) axo.Position(2)+0.08 axo.Position(3) axo.Position(4)-0.1];
    [ay, ax] = offsetaxis(axo, 'y', 0.02, 'x', 0.02);
    set(ax, 'tickdir', 'out', 'FontName', 'Helvetica')
    set(ay, 'tickdir', 'out', 'FontName', 'Helvetica')

    ylabel(ay, 'Estimation (zscore)')
    xlabel(ax, 'Sample mean (deg)')
    box off    
end

sname = sprintf('Estim_BIAS_SAMPMU%d_zscore_wofit_bins', nsmu);
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, sname, 'png')



%% plot estimation by sampmu bins
xhic = cell(2,1);
for c = 1:2
    for s = 1:length(subj)
        for gb = 1:length(sbins)
            xhic{c}(s, gb) = nanmean(mconfb{c}(sbins{gb}, s, 2)-mconfb{c}(sbins{gb}, s, 1));
        end
    end
end

%% Plot bias
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
sname = sprintf('Estim_BIAS_SAMPMU%d_zscore_bins', nsmu);
exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
       'BackgroundColor','none')
% exportgraphics(vis, [sname '.pdf'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, [sname '.png'])

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

xhic = cell(2,1);
for c = 1:2
    for s = 1:length(subj)
        for gb = 1:length(sbins)
            for lr = 1:2
                xhic{c}(lr, s, gb) = nanmean(mconfb{c}(sbins{gb}, s, lr));
            end
        end
    end
end

N = length(subj);

xpts = [1 2 4 5];
coltmp = [colslr{1}(1,:); colslr{2}(1,:); colslr{1}(1,:); colslr{2}(1,:)];

gb = 1;

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
xlim([0 6]); ylim([-1.5 1.5]); ylabel(sprintf('near-zero\nSample mean'))
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
        text(mean(xpts(c,:))-0.2,max(data(:))+0.6,sprintf('p = %1.3f', p))
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
xlim([0.5 2.5]); ylim([-1 3]); ylabel(sprintf('\\DeltaEstimation\n(Right-Left)'))
set(gca, 'xtick',xpts,'xticklabels', {'Cue','Choice'}, 'ytick', [-1 0 1 2],'tickdir', 'out')
sname = sprintf('%s/Estim_BIAS_SAMPMU%d_zscore_singlesub_xnz', fig_dir, nsmu)
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, [sname '.png'])