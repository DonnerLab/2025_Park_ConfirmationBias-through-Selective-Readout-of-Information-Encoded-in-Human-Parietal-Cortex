% relevant data: %s_pm%d_btS_btE_btSE_abs.mat, gFTII_%s_parcels_pm%d_%s_np%d_msh.mat
addpath('plot_maps')
addpath('BrewerMap')

conditions = {'attn', 'choice'};
condt = {'Cue','Choice'};
cic = {'Consistent','Inconsistent'};

load('subjects.mat')

%% parameters
sf = 160; 

xtime = -0.1:1/sf:0.7;
bs0 = find(xtime==0);
d1 = find(xtime >= 0.1-0.1/sf & xtime <= 0.1+0.1/sf);
d2 = find(xtime >= 0.5-0.1/sf & xtime <= 0.5+0.1/sf);

ARG.pcrit = 0.05;  
ARG.Nperm = 2000;
ARG.T_Range =[xtime(1) xtime(end)]; 
ARG.tw = [d1 d2];

vws = [-50 7; 40 2];

gmapp = brewermap(90, '*PuRd'); gmapp = flip(gmapp); 
gmapn = brewermap(90, '*PuBuGn'); gmap = cat(1, gmapn, gmapp);
gmap=[0.5.*[1 1 1]; gmap; 0.4.*[1 1 1]]; 
pmap = brewermap(180, '*PuBu'); pmap = flip(pmap);
pmap=[0.5.*[1 1 1]; pmap; 0.4.*[1 1 1]]; 
tmap = brewermap(180, '*RdYlBu'); tmap=[0.5.*[1 1 1]; tmap; 0.4.*[1 1 1]];
ptmap = brewermap(180, '*YlGnBu'); ptmap = flip(ptmap); ptmap=[0.5.*[1 1 1]; ptmap; 0.4.*[1 1 1]];



psmpt = {'s_1:s_2','s_3:s_4','s_5:s_6','s_7:s_8','s_9:s_{10}','s_{11}:s_{12}'};


%% PLOT raw
pvalt = cell(2, 2);
tvalt = cell(2, 2);
hh = cell(2, 2);
climt = 10^-2.*[-1.5 1.5];
ctick = 10^-2.*[-1.5 0 1.5];

BM_params
hemv = {'lateral','medial'};

clim = [-0.03 0.03];
pds = 4;
fgmu = 4;
msk = true;

hw = 4;
vw = hw*3/4;
[x_st, y_st, x_w, y_h] = multi_axes(hw, vw, 2, 2, [3 4], [-0.04, 0.04, 0.18 0], [0.001 0.001]);

hw1 = 4;
vw1 = 1.8*vw;
[x_st1, y_st1, x_w1, y_h1] = multi_axes(hw1, vw1, 1, 2, [3 4], [0.02, 0.04, 0.18 0], [0.001 0.001]);

for c = 1:2
    load(sprintf('%s_pm%d_btS_btE_btSE_abs.mat', condt{c}, fgmu), 'gCORR','ROI')
    datao = nan(180, 30, length(pds));
    for ii=1:size(gCORR,1)
        for cs = 1:2
            for pd = 1:length(pds)
                datao(ii, :, cs, pd) = permute(nanmean(gCORR{ii,1}{cs, pds(pd)}(3, d1:d2, :), 2), [3 1 2]);
            end
        end
    end
    for pd = 1:length(pds)
        data22 = nan(22, 30, 2); % ROI x subj x cic
        %%% get ROIs %%%
        for a = 1:22
            areas = glasser_group(a);
            aind = ismember(ROI, areas);
            data22(a, :, :) = permute(nanmean(datao(aind, :, :, pd), 1), [2 3 1]);
        end
        
        %%%%%% do con-incon %%%%%%
        vis1 = figure;
        vis1.Units = 'centimeters';
        vis1.Position = [2 2 hw1 vw1]; % - - width height
        annotation(vis1, 'textbox', 'Position',[0.05 0.9 0.1 0.1],'String', sprintf('%s, %s', condt{c}, psmpt{pds(pd)}), 'LineStyle', 'none', 'FontSize', 7)
        
        
        cfg = []; cfg.nperm = 5000; cfg.ttest = false; cfg.perm = 'unpaired'; cfg.tail = 'right';         
        % permutation test
        for a = 1:22
            con = permute(data22(a, :, 1), [2 3 1]); % this makes it subj x 1
            incon = permute(data22(a, :, 2), [2 3 1]); % this makes it subj x 1
            [pvalt{c,pd}(a,1), tvalt{c,pd}(a,1)] = perm_test_behav(abs(con), abs(incon), cfg);
        end
        % fdr threshold
        th_tval = tvalt{c,pd};
        hh{c,pd} = fdr(pvalt{c,pd}, 0.05);
%         th_tval(~hh{c,pd}) = nan;
        tmp = ones(22, 1);
        tmp(hh{c,pd},1) = nan;

        % convert back to 180 parcels
        th_tval180 = nan(180, 1); mask180 = nan(180, 1);
        for a = 1:22
            areas = glasser_group(a);
            aind = ismember(ROI, areas);
            th_tval180(aind, 1) = th_tval(a);
            mask180(aind, 1) = tmp(a);
        end
        cbt = '\DeltaRho';
        for hi = 1:2
            ax = axes(vis1, 'Units', 'centimeters', 'Position', [x_st1 y_st1(hi) x_w1, y_h1]); hold all;
            cfg = BM_cfg(ax, 'L', climt, vws(hi, :), tmap22, cbt, 'on', mdir);
            cfg.mask_alpha = 0.5;
            vwt = hemv{hi};
%             hm = plot_bm_opacity(cfg, th_tval180, mask180);
            hm = plot_bm(cfg, th_tval180);
            camlight(10, 10)
            hm.Label.Position = [4 0 0.5];
            hm.FontSize = 8;
            hm.FontName = 'Helvetica';
            hm.Position = [hm.Position(1)+0.13 0.08+hm.Position(2) hm.Position(3) 0.5*hm.Position(4)];
            hm.Ticks = ctick;
            hm.Ruler.Exponent = -2;
            if hi == 1
                th2 = title('Con-Incon');
                th2.Position = [th2.Position(1)+50 th2.Position(2:3)];
            end
        end
        sname = sprintf('BM_Decoding_%s_pd%d_cic_muabs_fdr95_right', conditions{c}, pds(pd));
        print(vis1, sname, '-dpng', '-r300')
%         exportgraphics(vis1, [sname '.eps'],'ContentType','vector',...
%         'BackgroundColor','none')
    end
end

