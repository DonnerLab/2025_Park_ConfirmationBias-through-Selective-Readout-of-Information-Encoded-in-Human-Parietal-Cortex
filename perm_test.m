function [index,pvalcic] = perm_test(x,y,ARG, cfg)
%%%
% x, y should be in subj x time
if cfg.smooth
    tl = size(x, 2);
    % smooth with gaussian
    w = gausswin(7, 0.2); % make it taper at the end
    w = repmat(w', size(x,1), 1);
    % pad data
    tmp1 = cat(2, repmat(mean(x(:,1:3),2),1,3), x, repmat(mean(x(:,end-2:end),2),1,3)); x = [];
    tmp2 = cat(2, repmat(mean(y(:,1:3),2),1,3), y, repmat(mean(y(:,end-2:end),2),1,3)); y = [];

    ii = 1;
    for t = 1:tl
        x(:, ii) = mean(w.*tmp1(:, t:t+6), 2);
        y(:, ii) = mean(w.*tmp2(:, t:t+6), 2);
        ii = ii+1;
    end
end
%%%%%%%%%% do permutation test %%%%%%%%%%
tmp = x-y;
tval = sqrt(size(tmp, 1)-1)*(mean(tmp, 1)./std(tmp, [], 1));
tmpc = cat(1, x, y);
tval_boot = zeros(size(tmp,2),ARG.Nperm);
for b=1:ARG.Nperm
    % randomly switch labels between two groups
    rs = cat(1, ones(size(tmp,1), 1), 2*ones(size(tmp,1), 1));
    rs = logical(rs(randperm(length(rs)))-1);
    tmp1v = tmpc(rs, :);
    tmp2v = tmpc(~rs, :);
    tmpv = tmp1v-tmp2v;
    tval_boot(:,b) = sqrt(size(tmp,1)-1)*(mean(tmpv, 1)./std(tmpv, [], 1));
end
cfg.critval = abs(tinv(ARG.pcrit/2,size(tmp,1)-1)); % 0.05 two-sided. critical cutoff value for cluster members if parametric
[PosClus,NegClus] = eegck_clusterstats(cfg,tval',tval_boot);
[~,index] = eegck_stripclusters(PosClus,NegClus,[size(tval, 2) 1]);
clear tmp tmpc
%%%%%%%%%%% test 0.1-0.5 window %%%%%%%%%%%
if isfield(ARG, 'tw')
nperm = ARG.Nperm; d1 = ARG.tw(1); d2 = ARG.tw(2);
mx = mean(x(:, d1:d2),2); my = mean(y(:, d1:d2),2); 
if cfg.ttest
    [~, pvalcic] = ttest(mx, my);
elseif ~cfg.ttest
    if strcmp(cfg.perm, 'paired')
        diff_orig = mx-my;
        mu_orig = mean(diff_orig);
        dummy_null = diff_orig-mu_orig;
        tmpc = cat(1, diff_orig, dummy_null);
    elseif strcmp(cfg.perm, 'unpaired')
        diff_orig = mx-my;
        mu_orig = mean(diff_orig);
        tmpc = cat(1, mx, my);
    end
    mu_boot = zeros(size(diff_orig,2),nperm);
    for b=1:nperm
        % randomly switch labels between baseline and data
        rs = cat(1, ones(size(diff_orig,1), 1), 2*ones(size(diff_orig,1), 1));
        rs = logical(rs(randperm(length(rs)))-1);
        tmp1v = tmpc(rs, :);
        tmp2v = tmpc(~rs, :);
        tmpv = tmp1v-tmp2v;
        mu_boot(:,b) = mean(tmpv);
    end
    pvalcic = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
end
end
end

