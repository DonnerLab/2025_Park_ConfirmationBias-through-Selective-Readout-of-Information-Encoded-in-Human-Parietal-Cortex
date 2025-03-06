% relevant data: raw behavioral data: Sxx_4decode_rtco.mat

load('subjects.mat')

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};
cic = {'con','incon', 'all'};
sess_id = {[1 2 3 4], [1 2], [3 4]};


cols1 = {};
cols1{1}{1} = cat(1, repmat([228,26,28], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{1}{2} = repmat([255 184 10], 6, 1)./255;
cols1{2}{1} = cat(1, repmat([55,126,210], 7, 1)./255, repmat([10 184 0], 6, 1)./255);
cols1{2}{2} = repmat([255 184 10], 6, 1)./255;
cols1{3}{1} = [0 0 0]./255;

%% find gen means, normalize data before combining across subjects (zscore)
% get behavioral data for gen mean
ugmu = union(-15:5:15, -14:14);
gbehav = cell(2,1);
edges = -3:0.1:3;

for c = 1:2
    for s = 1:length(subj)
        fname = sprintf('%s/P03/behav/%s_4decode_rtco.mat', mdir, subj{s});
        beh = load(fname); % behavioral data
        clear tind
        tind = ismember(beh.sess, sess_id{c+1});
        % get estimation data by gen mean
        estim = beh.estimations(tind);
        gmu = beh.gen_means(tind);
        egmu = gmu(~isnan(estim));
        estim = estim(~isnan(estim));
        zestim = zscore(estim);
        for g = 1:length(ugmu)
%             if sum(gmu==ugmu(g))~=0
                [gbehav{c}(:, g, s), edges1] = histcounts(zestim(egmu==ugmu(g)), edges);
%             else
%                 continue;
%             end
        end
    end
end

%% PLOT Heatmap for generative mean: Normalized
gry = colormap('gray');
gry = flipud(gry);

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 18 9];
for c = 1:2
    subplot(1,2,c);
    data = gbehav{c};
    data = sum(data, 3);
    imagesc(fliplr(ugmu), edges, data)
    set(gca, 'ytick', [-2 0 2], 'yticklabel', string([2 0 -2]), 'tickdir', 'out')
    colormap(gry)
    xlabel('Generative mean (\circ)')
    title(condt{c})
    if c == 1
        ylabel('Estimation (zscore)')
    end
    if c == 2
        cb = colorbar;
        cb.Position = [cb.Position(1)+0.05 cb.Position(2)+0.12 cb.Position(3) 0.8*cb.Position(4)];
        cb.Label.String = 'Counts';
    end
    box off
    ax = gca;
    ax.Position = [ax.Position(1)-0.02*c ax.Position(2)+0.12 ax.Position(3) 0.8*ax.Position(4)];
end

sname = sprintf('%s/Estimation_heatmap_gmu_zscore', fig_dir);
hgexport(vis, sname)
saveas(vis, sname, 'png')


%% get behavioral data for sample mean for heatmap: Normalized
sbins = -25:2:25;
edges = -3:0.1:3;
gsbehav = cell(2, 1); 

for c = 2:3
    for s = 1:length(subj)
        fname = sprintf('%s/P03/behav/%s_4decode_rtco.mat', mdir, subj{s});
        beh = load(fname); % behavioral data
        clear tind
        tind = ismember(beh.sess, sess_id{c});
        % get estimation data by sample mean
        estim = beh.estimations(tind);
        smu = beh.sample_means(tind);
        esmu = smu(~isnan(estim));
        estim = estim(~isnan(estim));
        zestim = zscore(estim);
        % get where each trial lies withing the smu bin
        [N,~,ugsmu] = histcounts(esmu, sbins);        
        ugsmu = ugsmu + 1;
        for g = 1:length(sbins)
            if sum(ugsmu==g)~=0
                gsbehav{c-1}(:, g, s) = histcounts(zestim(ugsmu==g), edges);
            else
                continue;
            end
        end
    end
end


%% PLOT Heatmap with sample mean: zscore
% gry = colormap('gray');
% gry = flipud(gry);

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 18 9];

for c = 1:2
    subplot(1,2,c);
    data = gsbehav{c};
    data = sum(data, 3);
    imagesc(fliplr(sbins), edges, data)
    set(gca, 'ytick', [-2 0 2], 'yticklabel', string([2 0 -2]), 'tickdir', 'out')
    colormap(gry)
    xlabel('Sample mean (\circ)')
    title(condt{c})
    if c == 1
        ylabel('Estimation (zscore)')
    end
    if c == 2
        cb = colorbar;
        cb.Position = [cb.Position(1)+0.08 cb.Position(2)+0.12 cb.Position(3) 0.8*cb.Position(4)];
        cb.Ticks = [0 50 100];
        cb.Label.String = 'Counts';        
    end
    box off
    ax = gca;
    ax.Position = [ax.Position(1)-0.02*c ax.Position(2)+0.12 ax.Position(3) 0.8*ax.Position(4)];
end

sname = sprintf('%s/Estimation_heatmap_smu_zscore', fig_dir);
hgexport(vis, sname)
saveas(vis, sname, 'png')

