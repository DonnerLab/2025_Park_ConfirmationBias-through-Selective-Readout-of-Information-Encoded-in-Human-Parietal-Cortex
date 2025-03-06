function [pval, mu_orig] = perm_test_behav(x,y,cfg)
%%%
% cfg.perm = 'paired', 'zero', 'unpaired', 
nperm = cfg.nperm; 

if cfg.ttest
    % paired t-ttest
    [~, pval] = ttest(x, y, 'tail', cfg.tail);
elseif strcmp(cfg.perm, 'paired')
    diff_orig = x-y;
    mu_orig = nanmean(diff_orig);
    dummy_null = diff_orig-mu_orig;
    tmpc = cat(1, diff_orig, dummy_null);
    mu_boot = zeros(size(diff_orig,2),nperm);
    for b=1:nperm
        % randomly switch labels between baseline and data
        rs = cat(1, ones(size(diff_orig,1), 1), 2*ones(size(diff_orig,1), 1));
        rs = logical(rs(randperm(length(rs)))-1);
        tmp1v = tmpc(rs, :);
        tmp2v = tmpc(~rs, :);
        tmpv = tmp1v-tmp2v;
        mu_boot(:,b) = nanmean(tmpv);
    end
    
elseif strcmp(cfg.perm, 'zero')
    mu_orig = nanmean(x);
    mu_boot = zeros(size(x, 2), nperm);
    for b=1:nperm
        % randomly assign signs
        rs = sign(rand(size(x, 1), 1)-0.5);
        mu_boot(:,b) = nanmean(rs.*x);
    end
    
elseif strcmp(cfg.perm, 'unpaired')
    diff_orig = x-y;
    mu_orig = nanmean(diff_orig);
    tmpc = cat(1, x, y);
    mu_boot = zeros(size(diff_orig, 2), nperm);
    for b=1:nperm
        % randomly switch labels between con and incon
        rs = cat(1, ones(size(diff_orig,1), 1), 2*ones(size(diff_orig,1), 1));
        rs = logical(rs(randperm(length(rs)))-1);
        tmp1v = tmpc(rs, :);
        tmp2v = tmpc(~rs, :);
        tmpv = tmp1v-tmp2v;
        mu_boot(:,b) = nanmean(tmpv);
    end
    
end
if ~cfg.ttest
    if strcmp(cfg.tail, 'both')
        pval = (sum(mu_boot>=abs(mu_orig)) + sum(mu_boot<=-abs(mu_orig)))/nperm;
    elseif strcmp(cfg.tail, 'right')
        pval = sum(mu_boot>=mu_orig)/nperm;
    elseif strcmp(cfg.tail, 'left')
        pval = sum(mu_boot<=mu_orig)/nperm;
    end
end
end