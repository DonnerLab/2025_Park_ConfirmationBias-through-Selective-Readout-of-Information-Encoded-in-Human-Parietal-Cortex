% relevant data: raw behavioral data: Sxx_4decode_rtco.mat

%% Parameters
conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};
cic = {'con','incon', 'all'};
sess_id = {[1 2 3 4], [1 2], [3 4]};

sf = 200;
etimep = -0.1:1/sf:1;
indx = [1:26];
rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

load('subjects.mat')

cols1 = {};
cols1{1}{1} = cat(1, repmat([228,26,28], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{1}{2} = repmat([255 184 10], 6, 1)./255;
cols1{2}{1} = cat(1, repmat([55,126,210], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{2}{2} = repmat([255 184 10], 6, 1)./255;
cols1{3}{1} = [0 0 0]./255;



%% plot psychometric curve (choice & cue)
colslr = {};
colslr{1} = [178,171,210; 94,60,153]./255; % left: light, dark, blueish
colslr{2} = [253,184,99; 230,97,1]./255; % right: light, dark, redish


% get psignifit toolbox: https://github.com/wichmann-lab/psignifit


%% get behavioral data
% get behavioral data for gen mean
ugmu = union(-15:5:15, -14:14);
gbehav = nan(2, length(ugmu), length(subj));
gbinary = cell(2, length(ugmu), length(subj));
for c = 2:3
    for s = 1:length(subj)
        fname = sprintf('%s_4decode_rtco.mat', subj{s});
        beh = load(fname); % behavioral data
        clear tind
        tind = ismember(beh.sess, sess_id{c});
        % get estimation data by gen mean
        estim = beh.estimations(tind);
        gmu = beh.gen_means(tind);
        egmu = gmu(~isnan(estim));
        estim = estim(~isnan(estim));
        zestim = zscore(estim);
        for g = 1:length(ugmu)
            if sum(gmu==ugmu(g))~=0
                gbehav(c-1, g, s) = mean(zestim(egmu==ugmu(g)));
                if c == 3
                    % get binary choice
                    ichoice = beh.binary_choices(tind); % -1, 1, 0 (for cue condition, correct button press), 99 (miss)
                    tmp = ichoice(gmu==ugmu(g));
                    tmp(tmp==99) = nan;
                    gbinary{c-1, g, s} = tmp;
                elseif c == 2
                    % get interm cue
                    icue = beh.interm_cues(tind); % -1, 1, NaN
                    tmp = icue(gmu==ugmu(g));
                    gbinary{c-1, g, s} = tmp;
                end
            else
                continue;
            end
        end
    end
end

nsub = nansum(~isnan(permute(gbehav(1, :, :), [2 3 1])),2);

%% get estimation error
ugmu = union(-15:5:15, -14:14);
geerr = nan(2, length(ugmu), length(subj)); 

for c = 2:3
    for s = 1:length(subj)
        fname = sprintf('%s_4decode_rtco.mat', subj{s});
        beh = load(fname); % behavioral data
        clear tind
        tind = ismember(beh.sess, sess_id{c});
        % get estimation data by gen mean
        eerr = beh.estim_resp_accu(tind);
        gmu = beh.gen_means(tind);
        egmu = gmu(~isnan(eerr));
        eerr = eerr(~isnan(eerr));
        for g = 1:length(ugmu)
            if sum(gmu==ugmu(g))~=0
                geerr(c-1, g, s) = nanmean(abs(eerr(egmu==ugmu(g))));
            else
                continue;
            end
        end
    end
end

% standard error of mean for estimation
squared_error = geerr.^2;
nnsub = repmat(nsub, 1, 2); nnsub = nnsub';
mean_squared_error = nansum(squared_error, 3)./nnsub;
fprintf('Cue: %1.2f deg\n', mean(sqrt(mean_squared_error(1, :))))
fprintf('Choice: %1.2f deg\n', mean(sqrt(mean_squared_error(2, :))))


% mean accuracy of binary choice
gbina = nan(length(ugmu), length(subj)); 
c = 3;
for s = 1:length(subj)
    fname = sprintf('%s_4decode_rtco.mat', subj{s});
    beh = load(fname); % behavioral data
    clear tind
    tind = ismember(beh.sess, sess_id{c});
    % get binary choice data by gen mean
    bina = beh.interm_resp_accu(tind);
    gmu = beh.gen_means(tind);
    egmu = gmu(~isnan(bina));
    bina = bina(~isnan(bina));
    for g = 1:length(ugmu)
        if sum(gmu==ugmu(g))~=0
            gbina(g, s) = nanmean(bina(egmu==ugmu(g)));
        else
            continue;
        end
    end
end

fprintf('Choice: %1.2f %%\n', 100*nanmean(nanmean(gbina, 2)))

%% get behavioral data for sample mean
edges = -27.5:1:27.5;
gsbehav = nan(2, length(edges), length(subj)); 
gsbinary = cell(2, length(edges), length(subj));
tmp = [];
for c = 2:3
    for s = 1:length(subj)
        fname = sprintf('%s_4decode_rtco.mat', subj{s});
        beh = load(fname); % behavioral data
        clear tind
        tind = ismember(beh.sess, sess_id{c});
        % get estimation data by sample mean
        estim = beh.estimations(tind);
        smu = beh.sample_means(tind);
        smub = mean(beh.samples(tind, 1:6), 2);
        esmu = smu(~isnan(estim));
        estim = estim(~isnan(estim));
        zestim = zscore(estim);
        [N,~,ugsmu] = histcounts(esmu, edges);

        if find(ugsmu==0)
            fprintf('bin again.')
        end
        for g = 1:length(edges)
            if sum(ugsmu==g)~=0
                gsbehav(c-1, g, s) = mean(zestim(ugsmu==g));
                [~,~,ugismu] = histcounts(smub, edges);
                if c == 3
                    % get binary choice
                    ichoice = beh.binary_choices(tind);
                    tmp = ichoice(ugismu==g);
                    tmp(tmp==99) = nan;
                    gsbinary{c-1, g, s} = tmp;
                elseif c == 2
                    % get interm cue
                    icue = beh.interm_cues(tind);
                    tmp = icue(ugismu==g);
                    gsbinary{c-1, g, s} = tmp;
                end
            else
                gsbinary{c-1, g, s} = [];
                continue;
            end
        end
    end
end

%% 2AFC choice & interm cue
% average error for each gen mean.
gpsy = cell(2, 1); %nan(2, length(ugmu), length(subj), 2, 2);
for c = 1:2
    for g = 1:length(ugmu)
        for s = 1:length(subj)
            if ~isempty(gbinary{c, g,s})
                gpsy{c}(g, s, 1, 1) = nansum(gbinary{c, g, s}==-1); % left
                gpsy{c}(g, s, 1, 2) = nansum(~isnan(gbinary{c, g, s})); % total
                gpsy{c}(g, s, 2, 1) = nansum(gbinary{c, g, s}==1); % right
                gpsy{c}(g, s, 2, 2) = nansum(~isnan(gbinary{c, g, s})); % total
            else
                gpsy{c}(g, s, :, :) = nan;
            end
        end
    end
end

%% plot psy curve Generative mean (both cue and choice)
plotOptions = [];
plotOptions.extrapolLength = 0.2;

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 3.5 3];

result = cell(2, 1);
xlm = [ugmu(1)-0.5 ugmu(end)+0.5];
ylm = [0 1];
for c = 1:2
    for lr = 2%1:2
        plot([0 0], ylm, ':', 'color', 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off')
        hold on;
        plot(ugmu, nanmean(gpsy{c}(:, :, lr, 1)./gpsy{c}(:, :, lr, 2), 2), 'ok', 'color', cols1{c}{1}(1, :),'markerfacecolor', cols1{c}{1}(1, :), 'markersize', 2, 'LineStyle', 'none', 'markeredgecolor', 'none')
        data = [ugmu' nanmean(gpsy{c}(:, :, lr, 1), 2) nanmean(gpsy{c}(:, :, lr, 2), 2)];
        if c == 2            
            options = [];
            options.sigmoidName = 'norm';
            options.expType = 'YesNo';
        elseif c == 1
            options = [];
            options.sigmoidName = 'norm';
            options.expType = 'nAFC';
            options.expN = 4;
            options.threshPC = 0.75;
        end
        result{c} = psignifit(data, options);
        
        xlength   = max(result{c}.data(:,1))-min(result{c}.data(:,1));
        x         = linspace(min(result{c}.data(:,1)),max(result{c}.data(:,1)),1000);
        xLow      = linspace(min(result{c}.data(:,1))-plotOptions.extrapolLength*xlength,min(result{c}.data(:,1)),100);
        xHigh     = linspace(max(result{c}.data(:,1)),max(result{c}.data(:,1))+plotOptions.extrapolLength*xlength,100);
        
        fitValuesLow    = (1-result{c}.Fit(3)-result{c}.Fit(4))*arrayfun(@(x) result{c}.options.sigmoidHandle(x,result{c}.Fit(1),result{c}.Fit(2)),xLow)+result{c}.Fit(4);
        fitValuesHigh   = (1-result{c}.Fit(3)-result{c}.Fit(4))*arrayfun(@(x) result{c}.options.sigmoidHandle(x,result{c}.Fit(1),result{c}.Fit(2)),xHigh)+result{c}.Fit(4);
        
        fitValues = (1-result{c}.Fit(3)-result{c}.Fit(4))*arrayfun(@(x) result{c}.options.sigmoidHandle(x,result{c}.Fit(1),result{c}.Fit(2)),x)+result{c}.Fit(4);
        hold on;
        hline = plot(x, fitValues, 'Color', cols1{c}{1}(1, :),'LineWidth',1);
        hold on;
        plot(xLow,  fitValuesLow,'--',  'Color', cols1{c}{1}(1, :),'LineWidth',1)
        hold on;
        plot(xHigh, fitValuesHigh,'--', 'Color', cols1{c}{1}(1, :),'LineWidth',1)
        hold on;
        if c == 1 % plot asymptotes
            plot(repmat([ugmu(1); ugmu(end)], 1, 2), [0.25 0.75; 0.25 0.75], ':m', 'LineWidth', 1)
        end
        xlim(xlm);
        ylim(ylm);
        hold on;        
    end
    hold on;
end
% sedin1 = 3:5:24;
ax = gca;
set(ax, 'FontName', 'Helvetica')
ax.Position = [ax.Position(1)+0.1 ax.Position(2)+0.1 ax.Position(3)-0.1 ax.Position(4)-0.1];
% [ay, ax] = offsetaxis(axo, 'y', 0.04, 'x', 0.04);
set(ax, 'xtick', ugmu(1:5:end),'xticklabels', string(ugmu(1:5:end)), 'tickdir', 'out', 'FontSize', 8, 'FontName', 'Helvetica')
set(ax, 'ytick', [0 1], 'tickdir', 'out', 'FontSize', 8, 'FontName', 'Helvetica')
ax.XAxis.TickLength = [0.03 0];
ax.YAxis.TickLength = [0.03 0];
% ax.XAxis.LineWidth = 1.5;
% ay.YAxis.LineWidth = 1.5;

ylabel(ax, sprintf('Portion\nright choice'), 'FontSize', 8)
xh = xlabel(ax, 'gen.\mu ({\circ})', 'FontSize', 8);
xh.Position = [-7 xh.Position(2) xh.Position(3)];
%     legend(axo, 'Left','Right', 'FontSize', 14)
%     legend(axo, 'boxoff')
box off
%     title(conditions{c})

fig_file = sprintf('confbias_gmu_cc_psycurve_right.eps');
hgexport(vis, fig_file)