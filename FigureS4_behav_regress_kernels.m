% relevant data: raw behavioral data: Sxx_4decode_rtco.mat


dots = {'s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','s_{11}','s_{12}'};
pdots = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

load('subjects.mat')

sess_id = {[1 2 3 4], [1 2], [3 4]};
condt = {'Cue','Choice'};

% colormap gray;
qtl_vals = [0.2, 0.6, 1.0];

udpt = '_udpool';

%% kernels with only one sample at a time
% original experimental data
estimbias = false;
if estimbias
    nr = 26;
    spr = 3;
else
    nr = 25;
    spr = 2;
end
nz = true;
if nz
    fgmu = 4;
    zgmu = -fgmu:fgmu;
else
    fgmu = 14;
end
ng = true; % if true, no gen mean regressor. to do this estimbias also has to be false
fpkse = cell(2, 1);
for c = 1:2    
    estim_kernels = nan(24, length(subjects)); estim_kernels_pval = nan(length(subjects));genmub = [];estimb = [];
    for sj = 1:length(subjects)
        
        data = load(sprintf('%s/behav/%s_4decode_rtco.mat', dat_dir, subjects{sj}));
        tind = ismember(data.sess, sess_id{c+1});
        if nz
            zgind = ismember(data.gen_means(tind), zgmu);
        end
        samples = [data.samples(tind, 1:6) data.samples(tind, 1:6) data.samples(tind, 7:12) data.samples(tind, 7:12)];
        genmu = data.gen_means(tind);
        if ~ng
            samples = samples-genmu;
        end
        samp_indx = [];
        for int = 1:2
            % Choice to sample consistency            
            if c == 2
                interm_sign = sign(data.binary_choices(tind));
                correct_trls = data.binary_choices(tind).*genmu>0;
            elseif c == 1
                interm_sign = sign(data.interm_cues(tind));
                correct_trls = data.interm_cues(tind).*genmu>0;
            end
                       
            consistent_samples = data.samples(tind, 1+6*(int-1):6*int).*repmat(interm_sign, 1, 6)>0;
            inconsistent_samples = data.samples(tind, 1+6*(int-1):6*int).*repmat(interm_sign, 1, 6)<0;
            
            samp_indx = [samp_indx consistent_samples inconsistent_samples];
        end
        b = []; 
        for d = 1:24
            % all trials
            if nz
                trls2use = ~isnan(data.estim_rt(tind)) & ~isnan(data.interm_rt(tind))&zgind&samp_indx(:, d);
            else
                trls2use = ~isnan(data.estim_rt(tind)) & ~isnan(data.interm_rt(tind))&samp_indx(:, d);
            end
            estim_responses = data.estimations(tind);
            estim_responses = estim_responses(trls2use);
            smpls = samples(trls2use, d);
            interm = interm_sign(trls2use);
            gmu = genmu(trls2use);

            if estimbias && ~ng
                [b(d,:),~,~,~,stats] = regress(estim_responses, [ones(length(estim_responses),1) smpls gmu interm]);
            elseif ~estimbias && ng
                [b(d,:),~,~,~,stats] = regress(estim_responses, [ones(length(estim_responses),1) smpls]);
            elseif ~estimbias && ~ng
                [b(d,:),~,~,~,stats] = regress(estim_responses, [ones(length(estim_responses),1) smpls gmu]);            
            end
        end
        estim_kernels_pval(sj, 1) = stats(3);
        estim_kernels(:, sj) = b(:,2);
        if ~ng
            genmub(:, sj) = b(:, 3); 
        end
        if estimbias
            estimb(:, sj) = b(:, 4); 
        end
    end
    for int = 1:2
        for cs = 1:2
            fpkse{c}(1, cs, 1+6*(int-1):6*int, :) = estim_kernels(12*(int-1)+(1+6*(cs-1):6*cs), :);
        end
    end
end
save(sprintf('%s/gFTPK_SE_int12_cic_pm4_np%d%s_CC.mat', dat_dir, npl, udpt), 'fpkse')
%% plot
fgmu = 4;

alp = 0.5;
mrs = 2;
lnw1 = 1;
lnw2 = 0.25;
bv1 = -0.08;
bv2 = -0.12;
bt = 0.04;

cols = [10 184 0; 255 184 10]./255; % con, incon
fdrp = 0.05;
npl = 1; f = 1;
tnpl = 12/npl;
tnp = tnpl/2; % divide by number of intervals(=2)
load(sprintf('%s/P03/behav/gFTPK_SE_int12_cic_pm4_np%d%s_CC.mat', mdir, npl, udpt), 'fpkse')
fw = 16;
fh = 5;
vis = figure;
vis.Units = 'centimeters';
vis.Position = [1 1 fw fh];
[x_st, y_st, x_w, y_h] = multi_axes(fw, fh, 2, 1, [1 2.3], [0.1 0.02 0.005 0.22]);
% annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('Gen\\mu: \\pm%d', fgmu), 'LineStyle', 'none', 'fontsize', 12)

cfg = []; 
cfg.nperm = 5000; 
cfg.ttest = false; 
cfg.perm = 'unpaired';
cfg.tail = 'both';
for c = 1:2
    % plot
    ax = axes('Units', 'centimeters', 'Position', [x_st(c), y_st, x_w, y_h]);
    for int = 1:2
        for cs = 1:2
            % absolute values
            plot([1+tnp*(int-1):tnp*int; 1+tnp*(int-1):tnp*int], [squeeze(mean(abs(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :)), 4))'-squeeze(std(abs(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :)), [], 4))'./sqrt(length(subjects)-1); ...
                squeeze(mean(abs(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :)), 4))'+squeeze(std(abs(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :)), [], 4))'./sqrt(length(subjects)-1)], 'LineStyle', '-','color', cols(cs, :));
            hold on;
            
            plot(1+tnp*(int-1):tnp*int, squeeze(mean(abs(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :)), 4)), 'marker', 'o', 'markersize', mrs, 'color', cols(cs, :), ...
                'markerfacecolor', cols(cs, :), 'markeredgecolor', 'none', 'LineWidth', lnw1)
            hold on;
        end
    end
    % do stats
    pval1 = [];
    for cs = 1:2           
        % from baseline
        for pd = 1:tnpl
            x = permute(abs(fpkse{c}(f, cs, pd, :)), [4 1 2 3]); y = zeros(length(subjects), 1);
            pval1(pd, cs) = perm_test_behav(x,y,cfg);
        end
        pval1(pval1(:, cs)==0, cs) = 10^-30;
    end

    for cs = 1:2
        h = fdr(pval1(:,cs), fdrp);

        bh = bar(find(h)-0.2*(cs-1), bv1+bt.*sign(pval1(h,cs)));
        bh.FaceColor = cols(cs, :);
        bh.EdgeColor = 'none';
        bh.BarWidth = 0.2;
        bh.BaseValue = bv1;
        hold on;
    end
    hold on;
    % con-incon difference
    pval = [];
    for pd = 1:tnpl
        x = permute(abs(fpkse{c}(f, 1, pd, :)), [4 1 2 3]); y = permute(abs(fpkse{c}(f, 2, pd, :)), [4 1 2 3]);       
        pval(pd,1) = perm_test_behav(x,y,cfg);
    end
    pval(pval==0) = 10^-15;
    h1 = fdr(pval, fdrp);
    b2 = bar(find(h1(1:12))+0.2, bv1+bt.*sign(pval(h1(1:12))));
    b2.FaceColor = [0 0 0];
    b2.EdgeColor = 'none';
    b2.BarWidth = 0.2;
    b2.BaseValue = bv1;
    
    hold on;
    plot([0.7 tnpl+0.3], [0 0], ':k', 'LineWidth', 0.5)
    hold on;
    axo = gca;
    if npl == 1
        set(axo, 'xtick', [1 6 7 12], 'xticklabels', dots([1, 6, 7, 12]), 'ytick', 0:0.2:0.4, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    else
        set(axo, 'xtick', 1:tnpl, 'xticklabels', pdots(1:tnpl), 'ytick', 0:0.2:0.4, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    end
    set(axo, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    axo.XAxis.TickLength = [0.02 0];
    axo.YAxis.TickLength = [0.02 0];
    xlim(axo, [0.5 tnpl+0.5]); ylim(axo, [bv2 0.51])
    if c == 1
        ylabel(axo, sprintf('Regression\ncoefficients'))
    end
    title(condt{c})
    box off
end
% 
% sname = sprintf('%s/Figure1_PKSE_int12_np%d%s_CC_pm%d_unpairedperm', fig_dir, npl, udpt, fgmu);
% hgexport(vis, sname)
% 
% saveas(vis, sname, 'png');

%% plot bar graph only for samples 7 and 8

diff_bias = cell(2, 1);
goodsubs = 1:30;
for c = 1:2
    diff_bias{c} = squeeze(abs(fpkse{c}(1, 1, :, goodsubs))-abs(fpkse{c}(1, 2, :, goodsubs)));
end

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 5 4];
data = [[diff_bias{1}(7, :) diff_bias{1}(8,:)]' [diff_bias{2}(7, :) diff_bias{2}(8,:)]'];
hold all;
bh = bar([1 2], mean(data, 1));
bh.FaceColor = 'flat';
bh.CData = [0.6 0.6 0.6; 0 0 0];
bh.EdgeColor = 'none';
plot([1 1], [mean(data(:, 1))-std(data(:, 1))/sqrt(60) mean(data(:, 1))+std(data(:, 1))/sqrt(60)], '-k')
plot([2 2], [mean(data(:, 2))-std(data(:, 2))/sqrt(60) mean(data(:, 2))+std(data(:, 2))/sqrt(60)], '-k')
box off; [h,p] = ttest(data(:, 1), data(:, 2));
upl = 1.2*mean(data(:, 2))+std(data(:, 2))/sqrt(60);
plot([1 2], upl*[1 1], '-k', 'LineWidth', 1); plot([1 1], [upl-0.02 upl], '-k'); plot([2 2], [upl-0.02 upl], '-k')
text(0.93, 1.3*mean(data(:, 2))+std(data(:, 2))/sqrt(60), sprintf('p=%1.4f',p), 'FontSize', 8, 'FontName', 'Helvetica')
set(gca, 'xtick', [1 2], 'xticklabels', {'Cue','Choice'}, 'ytick', 0:0.1:0.3, 'tickdir', 'out', 'FontName', 'Helvetica')
xlim([0 3]); ylim([0 0.3])
ylabel(sprintf('\\Delta|Regression \\beta|\nConsistent-Inconsistent'), 'FontSize', 8, 'FontName', 'Helvetica')

sname = sprintf('%s/diff_CC_PK_cic_s7s8.eps', fig_dir);
hgexport(vis, sname)

%% plot signed plot
fgmu = 4;

alp = 0.5;
mrs = 2;
lnw1 = 1;
lnw2 = 0.25;
bv1 = -0.37;
bv2 = -0.41;
bt = 0.04;

cols = [10 184 0; 255 184 10]./255; % con, incon
fdrp = 0.05;
npl = 1; f = 1;
tnpl = 12/npl;
tnp = tnpl/2; % divide by number of intervals(=2)
load(sprintf('%s/P03/behav/gFTPK_SE_int12_cic_pm4_np%d%s_CC.mat', mdir, npl, udpt), 'fpkse')

fw = 16;
fh = 5;
vis = figure;
vis.Units = 'centimeters';
vis.Position = [1 1 fw fh];
[x_st, y_st, x_w, y_h] = multi_axes(fw, fh, 2, 1, [1 2.3], [0.1 0.02 0.005 0.22]);
% annotation(vis, 'textbox', [0.05 0.9 0.1 0.1], 'string', sprintf('Gen\\mu: \\pm%d', fgmu), 'LineStyle', 'none', 'fontsize', 12)

cfg = []; 
cfg.nperm = 5000; 
cfg.ttest = false; 
cfg.perm = 'unpaired';
cfg.tail = 'both';
for c = 1:2
    % plot
    ax = axes('Units', 'centimeters', 'Position', [x_st(c), y_st, x_w, y_h]);
    for int = 1:2
        for cs = 1:2           
            % signed values
            plot([1+tnp*(int-1):tnp*int; 1+tnp*(int-1):tnp*int], [squeeze(mean(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :), 4))'-squeeze(std(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :), [], 4))'./sqrt(length(subjects)-1); ...
                squeeze(mean(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :), 4))'+squeeze(std(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :), [], 4))'./sqrt(length(subjects)-1)], 'LineStyle', '-','color', cols(cs, :));
            hold on;
            
            plot(1+tnp*(int-1):tnp*int, squeeze(mean(fpkse{c}(f, cs, 1+tnp*(int-1):tnp*int, :), 4)), 'marker', 'o', 'markersize', mrs, 'color', cols(cs, :), ...
                'markerfacecolor', cols(cs, :), 'markeredgecolor', 'none', 'LineWidth', lnw1)            
            hold on;
        end
    end
    % do stats
    pval1 = [];
    for cs = 1:2           
        % from baseline
        for pd = 1:tnpl
            x = permute(fpkse{c}(f, cs, pd, :), [4 1 2 3]); y = zeros(length(subjects), 1);
            pval1(pd, cs) = perm_test_behav(x,y,cfg);
        end
        pval1(pval1(:, cs)==0, cs) = 10^-30;
    end

    for cs = 1:2
        h = fdr(pval1(:,cs), fdrp);

        bh = bar(find(h)-0.2*(cs-1), bv1+bt.*sign(pval1(h,cs)));
        bh.FaceColor = cols(cs, :);
        bh.EdgeColor = 'none';
        bh.BarWidth = 0.2;
        bh.BaseValue = bv1;
        hold on;
    end
    hold on;
    % con-incon difference
    pval = [];
    for pd = 1:tnpl
        x = permute(fpkse{c}(f, 1, pd, :), [4 1 2 3]); y = permute(fpkse{c}(f, 2, pd, :), [4 1 2 3]);
        pval(pd,1) = perm_test_behav(x,y,cfg);
    end
    pval(pval==0) = 10^-30;
    h1 = fdr(pval, fdrp);
    b2 = bar(find(h1(1:12))+0.2, bv1+bt.*sign(pval(h1(1:12))));
    b2.FaceColor = [0 0 0];
    b2.EdgeColor = 'none';
    b2.BarWidth = 0.2;
    b2.BaseValue = bv1;
    
    hold on;
    plot([0.7 tnpl+0.3], [0 0], ':k', 'LineWidth', 0.5)
    hold on;
    axo = gca;
    if npl == 1
        set(axo, 'xtick', [1 6 7 12], 'xticklabels', dots([1, 6, 7, 12]), 'ytick', -0.4:0.2:0.4, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    else
        set(axo, 'xtick', 1:tnpl, 'xticklabels', pdots(1:tnpl), 'ytick', -0.4:0.2:0.4, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    end
    set(axo, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 7)
    axo.XAxis.TickLength = [0.02 0];
    axo.YAxis.TickLength = [0.02 0];
    xlim(axo, [0.5 tnpl+0.5]); ylim(axo, [bv2 0.41]); 
    if c == 1
        ylabel(axo, sprintf('Regression\ncoefficients'))
    end
    title(condt{c})
    box off
end

sname = sprintf('%s/Figure1_PKSE_int12_np%d%s_CC_pm%d_unpairedperm_signed', fig_dir, npl, udpt, fgmu);
% hgexport(vis, sname)
% 
% saveas(vis, sname, 'png');


