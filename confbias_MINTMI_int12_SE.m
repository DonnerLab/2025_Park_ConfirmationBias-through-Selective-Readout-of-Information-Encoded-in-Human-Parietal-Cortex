function confbias_MINTMI_int12_SE(jobin)

%% Parameters
sn = jobin{1};
% convert subject index to string
s = sprintf('S%s', sn);

if strcmp(jobin{2}, 'attn')
    c = 2;
elseif strcmp(jobin{2}, 'choice')
    c = 3;
end
fgmu = jobin{3}; % filter sample means

conditions = {'all', 'attn', 'choice'};

sess_id = {[1 2 3 4], [1 2], [3 4]}; % 1,2 is Cue(attn), 3,4 is Choice sessions

% number of samples to pool
npl = 1;
tnpl = round(12/npl);

% updown is pooled. 

%% compute MI between dots and estimation
fname = [s '_4decode_rtco.mat'];
beh = load(fname); % behavioral data
clear tind vind
tind = ismember(beh.sess, sess_id{c});

% for S
samples = beh.samples(tind, :); cind = nan(npl*size(samples, 1), tnpl);
% pool samples by npl
for pd = 1:tnpl
    tmp = samples(:, 1+npl*(pd-1):npl*pd);
    cind(:, pd) = tmp(:);
end
clear tmp

if c == 2
    interm = beh.interm_cues(tind);
elseif c == 3
    interm = beh.binary_choices(tind);
end

tiind = cell(3, 1); iind = cell(3, 1);
tiind{1} = [beh.samples(tind, 1:6).*repmat(interm, 1, 6) > 0 ~isnan(beh.consistent_samples(tind, 7:12))];
tiind{2} = [beh.samples(tind, 1:6).*repmat(interm, 1, 6) < 0 ~isnan(beh.inconsistent_samples(tind, 7:12))]; % inconsistent
tiind{3} = ~isnan(beh.samples(tind, :)); % all
for cs = 1:3
    for pd = 1:tnpl
        tmp = tiind{cs}(:, 1+npl*(pd-1):npl*pd);
        iind{cs}(:, pd) = tmp(:);
    end
end    
clear tiind

% valid trials based on interm button press
ind1 = beh.binary_choices(tind)~=99; % missed interm button press (includes both intervals)
ind1 = repmat(ind1, 1, npl);
ind1 = ind1(:);
% valid trials based on estimation
ind2 = ~isnan(beh.estimations(tind)); % missed estimation
% % here add one for the correct trials
% tpind = interm.*beh.gen_means(tind)>0;
% ind2 = ind2&tpind;
ind2 = repmat(ind2, 1, npl);
ind2 = ind2(:);

% get estimation data
tmp = beh.estimations(tind);
tmp = repmat(tmp, 1, npl);
estim = tmp(:);
clear tmp

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

% match number of trials for both Cue and Choice conditions
if npl == 1
    load('gminsp_un_pooled_np1.mat', 'gminsp_udpooled')
elseif npl == 2
    load('gminsp_un_pooled.mat', 'gminsp_udpooled')
end

load('subjects.mat')
sid = find(ismember(subj, s));
clear tmp
if fgmu > 0
    tmp(:, :, 1) = squeeze(gminsp_udpooled{fgmu+2, 1}(sid, :, :)); 
    tmp(:, :, 2) = squeeze(gminsp_udpooled{fgmu+2, 2}(sid, :, :)); 
elseif fgmu == 0
    tmp(:, :, 1) = squeeze(gminsp_udpooled{1, 1}(sid, :, :)); 
    tmp(:, :, 2) = squeeze(gminsp_udpooled{1, 2}(sid, :, :)); 
end
minsp = min(tmp, [], 3);
minsp = min(minsp, [], 1); minsp = min(minsp);


%% compute MI
nshuff = 100;
bi = 3;
opts = [];
opts.bias = 'naive';
opts.n_bins = {bi,bi};
opts.bin_method = {'eqpop'};
opts.computeNulldist= true;
opts.n_samples = nshuff;
opts.shuffling = {'A'};
opts.dim_shuffle = {'Trials'};
opts.parallel_sampling = false;
opts.supressWarnings = false;
opts.minsp = minsp;
opts.udpooled = udp;

outputlist = {'I(A;B)'};

sbsmp = 5;
opts.subsamp = 5; % number of random sub-sampling

MIv = nan(3, tnpl, sbsmp); 
MInull = nan(3, tnpl, nshuff, sbsmp);
for pd = 1:tnpl
    for cs = 1:3
        vid = iind{cs}(:, pd)&ind1&ind2&zgind;
        % get stimulus
        S = cind(vid, pd);
        % get estimation
        E = estim(vid);
        for sb = 1:sbsmp
            % match number of samples for con/incon
            % get number of excess consistent samples
            if length(S) >= minsp
                did = randperm(length(S), minsp);
            else
                did = randperm(length(S));
            end
            did = sort(did);
            did = did';
            St = S(did);
            Et = E(did);
            
            inputs = {St', Et'};
            
            %% Compute MI (binned method)
            [va,~,nulld] = MI(inputs, outputlist, opts);
            MIv(cs, pd, sb) = va{1};
            MInull(cs, pd, :, sb) = nulld;
        end
    end
end
MIv = nanmean(MIv, 3);
MInull = nanmean(MInull,4);
sname = sprintf('ISE_%s_%s_pm%d_np%d_CC.mat', s, conditions{c}, fgmu, npl)
save(sname, 'MIv', 'MInull','opts')
fprintf('\n')
end


