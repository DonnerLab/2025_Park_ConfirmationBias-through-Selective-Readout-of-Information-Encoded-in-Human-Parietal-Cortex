% relevant data: gFTII_%s_parcels_pm%d_%s_np%d_msh.mat, gFTII_%s_parcels_pm%d_%s_np%d_msh_PC1.mat

conditions = {'all', 'attn', 'choice'};
cic = {'con','incon'};
cict = {'Consistent','Inconsistent'};
condt = {'Cue','Choice'};

load('subjects.mat')

cols = {};
cols{1} = [10 184 0]./255;
cols{2} = [255 184 10]./255;

psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};

%% t-value
BM_params
alp = [0.01 0.05];
fgmu = 4; udpt = '_udpool'; npl = 2;
axw = 4;
axh = 1.8*axw;
load(sprintf('%s/P03/behav/gminsp_un_pooled.mat', mdir), 'gminsp_udpooled')
test_type = 'P'; 

clim = [-6.5 6.5];
ctick = [-6 0 6];
pvalt = cell(2, 1);
tvalt = cell(2, 1);
hh = cell(2, 1);

lgg = logical([1 0]);
msk = false; mask = {};
hemv = {'lateral','medial'};
tails = {'both','right'};
for al = 2%1:length(alp)
    for tl = 2%:length(tails)
        for ll = 1%:2
            thflag = lgg(ll); % threhold
            if thflag
                sup = 'th';
                cbt = 'II (bit,fdr)';
                lp = [4 4 0];
            else
                sup = [];
                cbt = 'II';
                lp = [4 6 0];
            end
            for c = 2:3
                vis = figure;
                vis.Units = 'centimeters';
                vis.Position = [2 2 axw axh]; % - - width height
                % axes matrix
                [x_st, y_st, x_w, y_h] = multi_axes(axw, axh, 1, 2, [3 4], [-0.1, 0.04, 0.18, 0], [0.001 0.001]);
                
                load(sprintf('gFTII_%s_parcels_pm%d_%s_np%d_msh_PC1.mat', conditions{c}, fgmu, udpt, npl), 'gFTIIpc1')
                % if PC1
                gFTII = gFTIIpc1;
                npd = size(gFTII{1,1}, 3);
                for pd = 4                    
                    if pd == 4
                        dtxt = 'onlys7s8';
                    elseif pd == 1
                        dtxt = 'onlys1s2';
                    end
                    minsp = min(gminsp_udpooled{fgmu+2, c-1}(:, :, pd), [], 2);
                    goodsub = minsp >= 81;
                    for a = 1:22
                        data22 = nan(sum(goodsub), 2);
                        areas = glasser_group(a);
                        aind = ismember(gFTII(:, 2), areas);
                        for cs = 1:2
                            tmp = cellfun(@(x) permute((x(cs, :, pd, goodsub)), [1 4 2 3]), gFTII(aind, 1), 'UniformOutput', false);
                            tmp = cell2mat(tmp);
                            tmp = nanmean(nanmean(tmp(:, :, d1:d2), 3), 1);
                            data22(:, cs) = tmp;
                        end

                        if strcmp(test_type, 'P')
                            ttt = 'Permutation test of mean difference';
                            cfg = []; cfg.nperm = 5000; cfg.perm = 'unpaired'; cfg.ttest = false; cfg.tail = 'right';                            
                            % permutation test                            
                            [pvalt{c-1}(a, 1), tvalt{c-1}(a, 1)] = perm_test_behav(data22(:, 1), data22(:, 2), cfg); 
                            % if PC1
                            clim = 10^-3.*[-2 2];
                            ctick = 10^-3.*[-2 2];                            
%                             clim = 10^-3.*[-1 1];
%                             ctick = 10^-3.*[-1 1];
                            cbt = '\DeltaII (bit)';
                        end
                    end
                    % fdr threshold
                    th_tval = tvalt{c-1};
                    if thflag && ~msk
                        hh{c-1} = fdr(pvalt{c-1}, alp(al));
                        th_tval(~hh{c-1}) = nan;
                    elseif thflag && msk
                        tmp = ones(22, 1);
                        hh{c-1} = fdr(pvalt{c-1}, alp(al));
                        tmp(hh{c-1}) = nan;
                    end
                    % convert back to 180 parcels
                    th_tval180 = nan(180, 1); mask180 = nan(180, 1);
                    for a = 1:22
                        areas = glasser_group(a);
                        aind = ismember(gFTII(:, 2), areas);
                        th_tval180(aind, 1) = th_tval(a);
                        mask180(aind, 1) = tmp(a);
                    end
                    for hi = 1:2
                        ax = axes(vis, 'Units', 'centimeters', 'Position', [x_st y_st(hi) x_w, y_h]); hold all;
                        cfg = BM_cfg(ax, 'L', clim, vws(hi, :), tmap22, cbt, 'on', mdir);
                        cfg.mask_alpha = 0.5;
                        vwt = hemv{hi};
%                         hm = plot_bm_opacity(cfg, th_tval180, mask180);
                        hm = plot_bm(cfg, th_tval180);
                        camlight(10, 10)
                        hm.Label.Position = [2 0 0.5];
                        hm.FontSize = 8;
                        hm.Position = [hm.Position(1)+0.13 0.08+hm.Position(2) 0.9*hm.Position(3) 0.5*hm.Position(4)];
                        hm.Ticks = ctick;
                        if hi == 1
                            title(sprintf('%s\nN_{subj} = %d',psmpt{pd},sum(goodsub)))
                        end
                    end
                end
                annotation(vis, 'textbox',[0.05 0.9 0.1 0.1], 'String', ...
                    sprintf('%s, %s (fdr %1.2f), %s tail, Consistent-Inconsistent. gen means: \\pm(%d:0), %s', condt{c-1}, ttt, alp(al), tails{tl}, fgmu, twt), 'LineStyle', 'none')
            
                set(findall(vis, '-property','FontName'),'FontName','Helvetica')
                set(findall(vis, '-property','FontSize'),'FontSize',8)
                
                sname = sprintf('GRP_%s_II_BM_pm%d%s_np%d%s_%s_a%d_%s_%s_%s_right_msh_PC1', conditions{c}, fgmu, udpt, npl, sup, tw, 100*(1-alp(al)), tails{tl}, test_type, dtxt)
%                 print(vis, sname, '-dpng', '-r300')
            end
        end
    end
end

