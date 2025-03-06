function cfg = BM_cfg(ax, hemi, clim, vws, cmap, cbt, cbar, mdir)

cfg = [];
cfg.mom_dir = mdir;
cfg.ah = ax;
cfg.hemi = hemi;
cfg.clim = clim;
cfg.view = vws;
cfg.cmap = cmap;
cfg.colorbar = cbar;
cfg.cbtext = cbt;
cfg.amb = 0.7; % strength of ambient light D: 0.3, [0 1]
cfg.diff = 0.4; % strength of diffuse light D: 0.6, [0 1]
cfg.spec = 0; % strength or mirror like reflection D: 0.9, [0 1]