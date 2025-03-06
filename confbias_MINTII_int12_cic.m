function confbias_MINTII_int12_cic(jobin)


% addpath(genpath('MINT'))

%% Parameters
sn = jobin{1};
% convert subject index to string
s = sprintf('S%s', sn);

if strcmp(jobin{2}, 'attn')
    c = 2;
elseif strcmp(jobin{2}, 'choice')
    c = 3;
end
a = jobin{3};
areas = glasser_group(a);

fgmu = jobin{4}; % filter generative means 

conditions = {'all', 'attn', 'choice'};

sess_id = {[1 2 3 4], [1 2], [3 4]}; % 1,2 is Cue(attn), 3,4 is Choice sessions


sf = 160;
etimep = -0.5:1/sf:4.5; % full length of ERF

%%%% define epochs
d0s = {0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 2.70, 2.85, 3.00, 3.15, 3.30, 3.45};
% get samples of onset times
d0 = cellfun(@(x) find(etimep>x-0.5/sf & etimep<x+0.5/sf), d0s, 'UniformOutput', false); 
d0 = cell2mat(d0);

epoch = [];
% for all dots, -0.1s ~ 0.7 s
for d = 1:12
    epoch(d, :) = [d0(d)-round(sf*0.1) d0(d) + round(sf*0.7)];    
end

% number of samples to pool
npl = 2;
tnpl = round(12/npl);

%%%% get behavioral data
fname = [s '_4decode_rtco.mat'];
beh = load(fname); % behavioral data
clear tind vind
tind = ismember(beh.sess, sess_id{c});

% for S
cind = []; tmp = [];
samples = beh.samples(tind, :);
% pool samples by npl
for pd = 1:tnpl
    tmp = samples(:, 1+npl*(pd-1):npl*pd);
    cind(:, pd) = tmp(:);
end

if c == 2
    interm = beh.interm_cues(tind);
elseif c == 3
    interm = beh.binary_choices(tind);
end

tiind = {}; iind = {}; tmp = [];
tiind{1} = [beh.samples(tind, 1:6).*repmat(interm, 1, 6) > 0 ~isnan(beh.consistent_samples(tind, 7:12))];
tiind{2} = [beh.samples(tind, 1:6).*repmat(interm, 1, 6) < 0 ~isnan(beh.inconsistent_samples(tind, 7:12))]; % inconsistent
for cs = 1:2
    for pd = 1:tnpl
        tmp = tiind{cs}(:, 1+npl*(pd-1):npl*pd);
        iind{cs}(:, pd) = tmp(:);
    end
end

% valid trials based on interm button press
ind1 = beh.binary_choices(tind)~=99; % missed interm button press (includes both intervals)
ind1 = repmat(ind1, 1, npl);
ind1 = ind1(:);
ind2 = ~isnan(beh.estimations(tind)); % missed estimation
ind2 = repmat(ind2, 1, npl);
ind2 = ind2(:);

% get estimation data
tmp = beh.estimations(tind); 
tmp = repmat(tmp, 1, npl);
estim = tmp(:);


if fgmu > 0
    % valid trials2: based on generative mean:
    zgmu = -fgmu:fgmu;
    tmp = ismember(beh.gen_means(tind), zgmu);
    tmp = repmat(tmp, 1, npl);
    zgind = tmp(:);
elseif fgmu == 0 % use all trials
    tmp = ~isnan(beh.gen_means(tind));
    tmp = repmat(tmp, 1, npl);
    zgind = tmp(:);
end

% finding the smallest number of dataset to compute
minsp = [];
for cs = 1:2
    minsp(cs, :) = sum(iind{cs}&ind1&ind2&zgind, 1);
end
minsp = min(minsp, [], 1);

%% compute II
nshuff = 100; % number of permutations
subsmp = 5; % number of sub-sampling
bi = 3;

opts = [];
opts.bin_method= {'eqpop'};
opts.n_bins = {bi,bi,bi};
opts.bias= 'naive';
opts.computeNulldist= true;
opts.redundancy_measure = 'I_min';
opts.shuffling ={'B'};
opts.n_samples = nshuff;
opts.dim_shuffle = {'Trials'};
opts.parallel_sampling = false;
opts.supressWarnings = true;
% opts.pid_constrained = true; 
opts.minsp = minsp;
p = 1;
t1 = tic;
for p = 1:length(areas)
    clear area
    if a ~= 26
        area = sprintf('HCPMMP1_%s',areas{p});
    else
        area = areas{p};
    end
    area = sprintf('HCPMMP1_%s',areas);
    % load PCA compon
    filename = sprintf('PCA_%s_%s_%s.mat', s, conditions{c}, area)
    load(fullfile(sprintf('/home/hamepark/P03/meg/PCA/%s/',s), filename), 'score')
    
    data1 = cellfun(@(x) permute(x, [3 1 2]), score, 'UniformOutput', false);
    data1 = cell2mat(data1);
    data1 = permute(data1, [2 3 1]); % ntrials x ncomp x ntimes
    data = nan(size(data1, 1), size(data1, 2), length(etimep));
    for tr = 1:size(data1, 1)
        [Y, Yx] = resample(permute(data1(tr, :, :), [3 2 1]), -0.5:1/200:4.5, 160);
        data(tr, :, :) = permute(Y, [2 1]);
    end
    
    % PCA comps are already downsampled to 200 Hz
    % PCA comps are downsampled to 160 Hz (17.4.2024)
    assert(size(data, 3)==length(etimep))
    % baseline correction (-500 ms to 0 ms)
    baseline = mean(data(:, :, 1:d0(1)), 3);
    baseline = repmat(baseline, [1, 1, size(data, 3)]);
    data = data-baseline;
    clear baseline
    tdata = [];
    for d = 1:12
        tdata(d, :, :, :) = data(:, :, epoch(d, 1):epoch(d, 2));
    end
    tdata = permute(tdata, [2 1 3 4]); % trials x samples x PC x time
    clear data
    IIs = cell(2, tnpl);
    for pd = 1:tnpl
        tmp = tdata(:, 1+npl*(pd-1):npl*pd, :, :);
        data = reshape(tmp, size(tmp, 1)*size(tmp, 2), size(tmp, 3), size(tmp, 4));
        for cs = 1:2
            vid = iind{cs}(:, pd)&ind1&ind2&zgind; % up/down pooled
            R = data(vid, :, :);
            % get stimulus
            S = cind(vid, pd);
            % get estimation
            E = estim(vid);
            [~, Npc, Nt] = size(R);
            IItmp = nan(Npc, Nt, nshuff+1, subsmp);
            t2 = tic;
            for sb = 1:subsmp
                % match number of samples for con/incon
                % get number of excess consistent samples
                if length(S) >= minsp(pd)
                    did = randperm(length(S), minsp(pd));
                else
                    did = randperm(length(S));
                end
                did = sort(did);
                did = did';
                Rt = R(did,:,:);
                St = S(did)';
                Et = E(did)';
                
                [Ntr, Npc, Nt] = size(Rt);
                %% Compute II (binned method)
                tic
                for i = 1 : Npc
                    for t = 1 : Nt
                        % S, R, E: nDimsX x nTrials
                        [II_biased, ~,IInull] = II({St; Rt(:, i, t)'; Et}, {'II(A,B,C)'},opts);
                        IItmp(i, t, 1, sb) = cell2mat(II_biased);
                        IItmp(i, t, 2:end, sb) = IInull;
                    end
                end
                toc
            end
            fprintf('\n')
            tm = toc(t2);
            IIs{cs, pd} = nanmean(IItmp, 4);
            fprintf('one condition: %1.2f', tm)
            fprintf('\n')
        end
    end
    sname = sprintf('II_%s_%s_%s_int12_pm%d_rtco_udpool_np%d_msh.mat', s, s, conditions{c}, areas{p}, fgmu, npl)
    save(sname, 'IIs', 'opts')
    fprintf('\n')
end
ttt = toc(t1);
fprintf('total computation time: %1.2f hours\n', ttt/3600)
end





