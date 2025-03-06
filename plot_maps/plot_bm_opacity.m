function hm = plot_bm_opacity(cfg, data, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot brain on top of anatomical template brain.
% cfg.ah = axes handle
% cfg.hemi = 'L' or 'R'; hemisphere
% cfg.clim = two element row vector of the color limits
% cfg.view = angles of view (ex. [270 0]; lateral, [-5 0]; caudal, [90
% 0]; medial
% cfg.cmap = colormap. put grey values on both extremes of map to color the
%               template brain. 
% cfg.colorbar = 'on' or 'off'
% cfg.cbtext = a string for the colorbar title
% cfg.amb = strength of ambient light D: 0.3, [0 1]
% cfg.diff = strength of diffuse light D: 0.6, [0 1]
% cfg.spec = strength or mirror like reflection D: 0.9, [0 1]
%
% data: a 180 element vector with the data to be plotted. NaNs are not
% plotted.
%
% returns handle of colorbar is colorbar is set to 'on'
% hm.Label.Position = [5 0 0.5];
% hm.Ruler.Exponent = -3;
% hm.FontSize = 9;
% hm.Position = [hm.Position(1)+0.05 hm.Position(2)+0.01 0.8*hm.Position(3) 0.5*hm.Position(4)];
% hm.Ticks = 10^-3.*[0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% addpath('/home/hamepark/P03/P03-code')
% addpath('/home/hamepark/P03/P03-code/plot_maps')
% addpath('/home/hamepark/P03/BrewerMap')

% load this instead of the Glasser label thing
% this at.cot.struct_names' has the label names with the hyphen instead of the underscore so no need to correct the ROI names. equivalent of Glasser_labels.csv.
load(sprintf('%s/P03/P03-code/plot_maps/atlas_%s.mat', cfg.mom_dir, cfg.hemi), 'at')
% this is for the face and vertices
load(sprintf('%s/P03/P03-code/plot_maps/template_inflated_%s.mat', cfg.mom_dir, cfg.hemi), 'tinf')    
% this is the grey anatomical brain
load(sprintf('%s/P03/P03-code/plot_maps/greybrain_%s.mat', cfg.mom_dir, cfg.hemi), 'gb')

n_vertices = size(tinf.vertices, 1);
areas = at.cot.struct_names(2:end);
strid = at.cot.table(2:end, 5); % size of # ROI, this is the equivalent of the atlas.indexmax
labels = at.label;
% assign values at each far end, so that we can use one colormap
gb = sign(gb);
gb(gb < 0) = -150;
gb(gb > 0) = 150;

% ROIs=Glasser_roi; % list with 180 names of Glasser ROIs
ROIs =  {};
for a = 1:22
    ROIs = cat(2, ROIs, glasser_group(a));
end
% create a single list with 360 ROIs (180 left, 180 right hemisphere)
Areas_labels=cell(length(ROIs), 1);
for r=1:length(ROIs)
    Areas_labels{r}=[cfg.hemi '_' ROIs{r} '_ROI'];
end

% xlm=[];ylm=[];zlm=[];
% [xlm(1) xlm(2)] = bounds(tinf.vertices(:, 1));
% [ylm(1) ylm(2)] = bounds(tinf.vertices(:, 2));
% [zlm(1) zlm(2)] = bounds(tinf.vertices(:, 3));
% xlm = ceil(1.1*xlm); ylm = ceil(1.1*ylm); zlm = ceil(1.1*zlm);

%% plot t-values -> on same axes -> this works
% put data into vertices
Vertices_values = NaN(n_vertices, 1);
nidAr = cell(180, 1);
for ii=1:size(Areas_labels,1)
    ROI = Areas_labels{ii};
    id_Area = strid(strcmp(areas, ROI));
    if ~isempty(id_Area)
        idAr = find(labels==id_Area);
        nidAr{ii} = idAr;
        Vertices_values(idAr) = data(ii);
    elseif isempty(id_Area)
        fprintf('empty:%d, roi: %s\n', ii, ROI)
    end
end

% underlay the grey template brain
h1 = trisurf(tinf.faces, tinf.vertices(:,1), tinf.vertices(:,2), tinf.vertices(:,3), gb, 'edgealpha', 0);
if isfield(cfg, 'amb')
    h1.AmbientStrength = cfg.amb;
else
    h1.AmbientStrength = 0.3;
end
if isfield(cfg, 'diff')
    h1.DiffuseStrength = cfg.diff;
else
    h1.DiffuseStrength = 0.6;
end
if isfield(cfg, 'spec')
    h1.SpecularStrength = cfg.spec;
else
    h1.SpecularStrength = 0.9;
end
view(cfg.ah, cfg.view); 
xlim(cfg.ah, [-100 100]); zlim(cfg.ah, [-100 100]); ylim(cfg.ah, [-120 120]) % this is quite specific to the template we are using
% xlim(cfg.ah, xlm); zlim(cfg.ah, zlm); ylim(cfg.ah, ylm) % this is quite specific to the template we are using
axis off

hold all 

% plot the data map 
h2 = trisurf(tinf.faces, tinf.vertices(:,1), tinf.vertices(:,2), tinf.vertices(:,3), Vertices_values, 'edgealpha', 0);
if isfield(cfg, 'amb')
    h2.AmbientStrength = cfg.amb;
else
    h2.AmbientStrength = 0.3;
end
if isfield(cfg, 'diff')
    h2.DiffuseStrength = cfg.diff;
else
    h2.DiffuseStrength = 0.6;
end
if isfield(cfg, 'spec')
    h2.SpecularStrength = cfg.spec;
else
    h2.SpecularStrength = 0.9;
end

colormap(cfg.ah, cfg.cmap)
caxis(cfg.ah, cfg.clim)
axis off

view(cfg.ah, cfg.view);
xlim(cfg.ah, [-100 100]); zlim(cfg.ah, [-100 100]); ylim(cfg.ah, [-120 120]) % this is quite specific to the template we are using
% xlim(cfg.ah, xlm); zlim(cfg.ah, zlm); ylim(cfg.ah, ylm) % this is quite specific to the template we are using

% plot the mask 
% put mask into vertices
mask_values = NaN(n_vertices, 1);
nidAr = cell(180, 1);
for ii=1:size(Areas_labels,1)
    ROI = Areas_labels{ii};
    id_Area = strid(strcmp(areas, ROI));
    if ~isempty(id_Area)
        idAr = find(labels==id_Area);
        nidAr{ii} = idAr;
        mask_values(idAr) = mask(ii);
    elseif isempty(id_Area)
        fprintf('empty:%d, roi: %s\n', ii, ROI)
    end
end
hold all;
h3 = trisurf(tinf.faces, tinf.vertices(:,1), tinf.vertices(:,2), tinf.vertices(:,3), mask_values, 'edgealpha', 0);
h3.FaceAlpha = cfg.mask_alpha;
if isfield(cfg, 'amb')
    h3.AmbientStrength = cfg.amb;
else
    h3.AmbientStrength = 0.3;
end
if isfield(cfg, 'diff')
    h3.DiffuseStrength = cfg.diff;
else
    h3.DiffuseStrength = 0.6;
end
if isfield(cfg, 'spec')
    h3.SpecularStrength = cfg.spec;
else
    h3.SpecularStrength = 0.9;
end

% colormap(h3, [1 1 1])
% caxis(h3, [0 1])
% axis off
% 
% view(cfg.ah, cfg.view);
% xlim(cfg.ah, [-100 100]); zlim(cfg.ah, [-100 100]); ylim(cfg.ah, [-120 120]) % this is quite specific to the template we are using

axis off
hm = [];
if strcmp(cfg.colorbar, 'on')
    hm = colorbar;
    hm.Box = 'off';
    hm.LineWidth = 0.1;
    hm.TickDirection = 'out';
    hm.TickLength = 0;
    hm.Label.String = cfg.cbtext;
    hm.Limits = cfg.clim;
end


end




