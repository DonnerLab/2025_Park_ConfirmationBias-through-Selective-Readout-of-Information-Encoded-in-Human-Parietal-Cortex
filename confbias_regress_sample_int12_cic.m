function confbias_regress_sample_int12_cic(jobin)


%% Parameters

sn = jobin{1};
% convert subject index to string
s = sprintf('S%s', sn);

if strcmp(jobin{3}, 'attn')
    c = 2;
%     parind = 2;
elseif strcmp(jobin{3}, 'choice')
    c = 3;
%     parind = 1;
end
a = jobin{4};
areas = glasser_group(a);


fgmu = jobin{5}; % filter generative means 

conditions = {'all', 'attn', 'choice'};

sess_id = {[1 2 3 4], [1 2], [3 4]};

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

%% get behavioral data
fname = [s '_4decode_rtco.mat'];
beh = load(fname); % behavioral data
clear tind vind
tind = ismember(beh.sess, sess_id{c});

% % up/down indices up:1, down:0
% clear tmp vind
% tmp = logical(beh.ref_angles(tind)); % same length as size(comb_dict{1}.erfdata, 1)
% vind(:, 1) = tmp;
% vind(:, 2) = ~tmp;

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
tiind{3} = ~isnan(beh.samples(tind, :)); % all
for cs = 1:length(tiind)
    for pd = 1:tnpl
        tmp = tiind{cs}(:, 1+npl*(pd-1):npl*pd);
        iind{cs}(:, pd) = tmp(:);
    end
end

% valid trials based on interm button press
ind1 = beh.binary_choices(tind)~=99; % missed interm button press (includes both intervals)
ind1 = repmat(ind1, 1, npl);
ind1 = ind1(:);
ind2 = ~isnan(beh.estimations(tind)); % missed estimation17.9
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
for cs = 1:length(iind)
    minsp(cs, :) = sum(iind{cs}&ind1&ind2&zgind, 1);
end
minsp = min(minsp, [], 1);

%% Regress
Nsh = 100;
subsmp = 5;

t1 = tic;
for p = 1:length(areas)
RegAngStim = cell(length(iind), tnpl); RegAngEstim = cell(length(iind), tnpl);
RegAngStimsh = cell(length(iind), tnpl); RegAngEstimsh = cell(length(iind), tnpl);

CVSCORE = cell(length(iind), tnpl);
CVPRED = cell(length(iind), tnpl);

clear area
if a ~= 26
    area = sprintf('HCPMMP1_%s',areas{p});
else
    area = areas{p};
end

% load PCA compon
filename = sprintf('PCA_%s_%s_%s.mat', s, conditions{c}, area)
load(fullfile(sprintf('PCA/%s/',s), filename), 'score')

data1 = cellfun(@(x) permute(x, [3 1 2]), score, 'UniformOutput', false);
data1 = cell2mat(data1);
data1 = permute(data1, [2 3 1]); % ntrials x ncomp x ntimes
data = nan(size(data1, 1), size(data1, 2), length(etimep));
for tr = 1:size(data1, 1)
    [Y, Yx] = resample(permute(data1(tr, :, :), [3 2 1]), -0.5:1/200:4.5, 160);
    data(tr, :, :) = permute(Y, [2 1]);
end

% PCA comps are already downsampled to 200 Hz
% PCA comps are further downsampled to 160 Hz 
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
w = 1;
for pd = 1:tnpl
    tmp = tdata(:, 1+npl*(pd-1):npl*pd, :, :);
    data = reshape(tmp, size(tmp, 1)*size(tmp, 2), size(tmp, 3), size(tmp, 4));
    for cs = 1:length(iind)
        vid = iind{cs}(:, pd)&ind1&ind2&zgind; % up/down pooled
        R = data(vid, :, :);
        % get stimulus & estimation
        S = cind(vid, pd);
        E = estim(vid);

        [~, Npc, Nt] = size(R);
        Bs = nan(2,Npc, Nt, subsmp); STATSs = nan(4, Npc, Nt, subsmp);
        Be = nan(2,Npc, Nt, subsmp); STATSe = nan(4, Npc, Nt, subsmp);
        Bshs = nan(2,Npc, Nt, Nsh, subsmp);
        Bshe = nan(2,Npc, Nt, Nsh, subsmp);

        cvscoret = nan(10, Nt, Npc, 3, subsmp);

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
            Rt = R(did, :, :);
            St = S(did);
            Et = E(did);

            [Ntr, Npc, Nt] = size(Rt);
            % Regress
            t2 = tic;

            for t = 1 : Nt
                % shuffle
                for h = 1:Nsh
                    ni = randperm(Ntr);
                    % S, R, E: nDimsX x nTrials
                    for cp = 1:Npc
                        Bshs(:, cp, t, h, sb) = regress(St(ni), [10^12.*Rt(:, cp, t) ones(size(Rt, 1), 1)]);
                        Bshe(:, cp, t, h, sb) = regress(Et(ni), [10^12.*Rt(:, cp, t) ones(size(Rt, 1), 1)]);
                    end
                end

                % normal version
                % S, R, E: nDimsX x nTrials
                for cp = 1:Npc
                    [Bs(:, cp, t, sb),~,~,~,STATSs(:, cp, t, sb)] = regress(St, [10^12.*Rt(:, cp, t) ones(size(Rt, 1), 1)]);
                    [Be(:, cp, t, sb),~,~,~,STATSe(:, cp, t, sb)] = regress(Et, [10^12.*Rt(:, cp, t) ones(size(Rt, 1), 1)]);
                end
                w = w+1;

                if mod(w,10)==0
                    fprintf('*')
                    if mod(w,100)==0
                        fprintf('\n')
                    end
                end
            end

            CVO = cvpartition(Ntr, 'KFold', 10);
            cvpredts = []; cvpredte = [];
            v = 1;
            for i = 1:CVO.NumTestSets
                trIdx = CVO.training(i);
                teIdx = CVO.test(i);

                cvBs = []; cvBe = [];
                cvpreds = nan(CVO.TestSize(i), Nt, Npc); cvprede = nan(CVO.TestSize(i), Nt, Npc); 
                for t = 1 : Nt
                    % S, R, E: nDimsX x nTrials
                    for cp = 1:Npc
                        cvBs(:, cp) = regress(St(trIdx), [10^12.*Rt(trIdx, cp, t) ones(CVO.TrainSize(i), 1)]);
                        cvBe(:, cp) = regress(Et(trIdx), [10^12.*Rt(trIdx, cp, t) ones(CVO.TrainSize(i), 1)]);
    
                        cvpreds(:, t, cp) = Rt(teIdx, cp, t)*cvBs(1, cp)+cvBs(end, cp);
                        cvscoret(i, t, cp, 1, sb) = corr(St(teIdx), cvpreds(:, t, cp));
    
                        cvprede(:, t, cp) = Rt(teIdx, cp, t)*cvBe(1, cp)+cvBe(end, cp);
                        cvscoret(i, t, cp, 2, sb) = corr(Et(teIdx), cvprede(:, t, cp));
                        
                        % correlation between shat and E
                        cvscoret(i, t, cp, 3, sb) = corr(Et(teIdx), cvpreds(:, t, cp));
                    end
                    v = v+1;

                    if mod(v,10)==0
                        fprintf('.')
                        if mod(w,100)==0
                            fprintf('\n')
                        end
                    end
                end
                cvpredts = cat(1, cvpredts, cvpreds);
                cvpredte = cat(1, cvpredte, cvprede);
            end
        end
        RegAngStim{cs, pd}{1} = Bs;
        RegAngStim{cs, pd}{2} = STATSs;
        RegAngStimsh{cs, pd} = Bshs;

        RegAngEstim{cs, pd}{1} = Be;
        RegAngEstim{cs, pd}{2} = STATSe;
        RegAngEstimsh{cs, pd} = Bshe;
        CVSCORE{cs, pd} = cvscoret;
        CVPRED{cs, pd} = cat(3, cvpredts, cvpredte);
        fprintf('\n')
        tm = toc(t2);
        fprintf('one condition: %1.2f', tm)
        fprintf('\n')
    end
end
sname = sprintf('RegSamp_%s_%s_%s_int12_pm%d_rtco_udpool_np%d.mat', s, conditions{c}, areas{p}, fgmu, npl)
save(sname, 'RegAngStim', 'RegAngEstim', 'RegAngStimsh', 'RegAngEstimsh', 'CVSCORE', 'CVPRED')

end
ttt = toc(t1);
fprintf('total time: %1.2f hrs', ttt/3600)
fprintf('\n')
end
