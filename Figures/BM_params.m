
addpath(sprintf('%s/P03/P03-code/plot_maps', mdir))
addpath(sprintf('%s/BrewerMap', mdir))


sf = 160;
etimep = -0.5:1/sf:4.5; % full length of ERF
xtime = -0.1:1/sf:0.7;
d0 = find(xtime == 0);

d1s = 0.1;
d2s = 0.5;
twt = sprintf('%1.1f s - %1.1f s', d1s, d2s);
tw = sprintf('%dmsto%dms',1000*d1s,1000*d2s);
d1 = find(xtime <= d1s + 0.5/sf & xtime >= d1s - 0.5/sf);
d2 = find(xtime <= d2s + 0.5/sf & xtime >= d2s - 0.5/sf);

vws = [-42 7; 40 2];
% vws = [-140 10; 40 2]; % frontal view


rmapp = brewermap(90, '*PuRd'); rmapp = flip(rmapp); 
rmapn = brewermap(90, '*PuBuGn'); rmap = cat(1, rmapn, rmapp);
rmap=[0.5.*[1 1 1]; rmap; 0.4.*[1 1 1]]; 

tmap = brewermap(180, '*RdYlBu'); tmap=[0.5.*[1 1 1]; tmap; 0.4.*[1 1 1]];
ptmap = brewermap(180, '*Purples'); ptmap = flip(ptmap); ptmap=[0.5.*[1 1 1]; ptmap; 0.4.*[1 1 1]];


rmapp11 = brewermap(13, '*PuRd'); rmapp11 = flip(rmapp11); 
rmapn11 = brewermap(13, '*PuBuGn'); rmap22 = cat(1, rmapn11, rmapp11);
rmap22=[0.5.*[1 1 1]; rmap22; 0.4.*[1 1 1]]; 

tmap22 = brewermap(26, '*RdYlBu'); tmap22=[0.5.*[1 1 1]; tmap22; 0.4.*[1 1 1]];
ptmap22 = brewermap(26, '*Reds'); ptmap22 = flip(ptmap22); ptmap22=[0.5.*[1 1 1]; ptmap22; 0.4.*[1 1 1]];
