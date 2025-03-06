%% Parameters

load('subjects.mat')

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};

dots = {'s_1','s_2','s_3','s_4','s_5','s_6','s_7','s_8','s_9','s_{10}','s_{11}','s_{12}'};
pdots = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOT FULL TRIAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alp = 0.5;
mrs = 2;
lnw1 = 1;
lnw2 = 0.25;
bv1 = -0.08;
bv2 = -0.12;
bt = 0.04;

cols = [10 184 0; 255 184 10]./255; % con, incon


fw = 30;
fh = 7.5;
vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 fw fh];

fgmu = 4;
npl = 1;
tnpl = 12/npl;
tnp = tnpl/2; % divide by number of intervals(=2)
pl = true;
if pl
    udpt = '_udpool';
elseif ~pl
    udpt = []; 
end

nperm = 5000;

load(sprintf('gFTMI_SE_int12_cic_pm4_np%d%s_CC.mat', npl, udpt), 'fmise')

for c = 1:2
    for f = 1:length(fgmu)
        [x_st, y_st, x_w, y_h] = multi_axes(fw, fh, 2, 1, [1 2.4], [0.1 0.02 0.002 0.22]);
        ax = axes('Units', 'centimeters', 'Position', [x_st(c), y_st(f), x_w, y_h]);
        for int = 1:2
            for cs = 1:2
                plot([1+tnp*(int-1):tnp*int; 1+tnp*(int-1):tnp*int], [squeeze(nanmean(fmise{c}(f, cs, 1+tnp*(int-1):tnp*int, :), 4))'-squeeze(nanstd(fmise{c}(f, cs, 1+tnp*(int-1):tnp*int, :), [], 4))'./sqrt(length(subj)-1); ... 
                    squeeze(nanmean(fmise{c}(f, cs, 1+tnp*(int-1):tnp*int, :), 4))'+squeeze(nanstd(fmise{c}(f, cs, 1+tnp*(int-1):tnp*int, :), [], 4))'./sqrt(length(subj)-1)], 'LineStyle', '-','color', cols(cs, :));
                hold on;

                plot(1+tnp*(int-1):tnp*int, squeeze(nanmean(fmise{c}(f, cs, 1+tnp*(int-1):tnp*int, :), 4)), 'marker', 'o', 'markersize', mrs, 'color', cols(cs, :), ...
                    'markerfacecolor', cols(cs, :), 'markeredgecolor', 'none', 'LineWidth', lnw1)
                hold on;
            end
            hold on;
        end
        % do stats        
        for cs = 1:2
            pval = []; 
            % from baseline                
            for pd = 1:tnpl
                tmp1e = squeeze(fmise{c}(f, cs, pd, :));
                mu_orig = mean(tmp1e);
                mu_boot = zeros(size(tmp1e, 2), nperm);
                for b=1:nperm
                    % randomly assign signs
                    rs = sign(rand(size(tmp1e, 1), 1)-0.5);
                    mu_boot(:,b) = mean(rs.*tmp1e);
                end
                pval(pd) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
            end
            pval(pval==0) = 10^-30;
            h = fdr(pval, 0.01);
            bh = bar(find(h)-0.2*(cs-1), bv1+bt.*sign(pval(h)));
            bh.FaceColor = cols(cs, :);
            bh.EdgeColor = 'none';
            bh.BarWidth = 0.2;
            bh.BaseValue = bv1;
            hold on;
        end
        hold on;
        % con-incon difference
        pvalcic = [];
        for pd = 1:tnpl
            tmp1e = squeeze(fmise{c}(f, 1, pd, :));
            tmp2e = squeeze(fmise{c}(f, 2, pd, :));
            tmpe = tmp1e-tmp2e;
            mu_orig = mean(tmpe);
            tmpce = cat(1, tmp1e, tmp2e);
            mu_boot = zeros(size(tmpe, 2), nperm);
            for b=1:nperm
                % randomly switch labels between con and incon
                rs = cat(1, ones(size(tmpe,1), 1), 2*ones(size(tmpe,1), 1));
                rs = logical(rs(randperm(length(rs)))-1);
                tmp1v = tmpce(rs, :);
                tmp2v = tmpce(~rs, :);
                tmpv = tmp1v-tmp2v;
                mu_boot(:,b) = mean(tmpv);
            end
            pvalcic(pd) = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
        end    
        pvalcic(pvalcic==0) = 10^-12;
        h1 = fdr(pvalcic, 0.01);
        b2 = bar(find(h1)+0.2, bv1+bt.*sign(pvalcic(h1)));
        b2.FaceColor = [0 0 0];
        b2.EdgeColor = 'none';
        b2.BarWidth = 0.2;
        b2.BaseValue = bv1;       

        hold on;
        plot([0.7 tnpl+0.3], [0 0], ':k', 'LineWidth', 0.5)
        hold on;
        axo = gca;
        if npl == 1
            set(axo, 'xtick', [1 6 7 12], 'xticklabels', dots([1, 6, 7, 12]), 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 21)
        else
            set(axo, 'xtick', 1:tnpl, 'xticklabels', pdots(1:tnpl), 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 21)
        end
        set(axo, 'tickdir', 'out', 'FontName', 'Helvetica', 'FontSize', 21)
        axo.XAxis.TickLength = [0.02 0];
        axo.YAxis.TickLength = [0.02 0];
        xlim(axo, [0.5 tnpl+0.5]); ylim(axo, [bv2 0.51])
        ylabel(axo, 'I(S;E) (bit)')
        title(condt{c})
        box off
    end
end
sname = sprintf('ISE_int12_pm%d_np%d%s_CC', fgmu(f), npl, udpt);
% hgexport(vis, sname)
% saveas(vis, sname, 'png')
