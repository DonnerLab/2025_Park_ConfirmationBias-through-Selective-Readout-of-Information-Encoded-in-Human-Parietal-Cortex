% Script to simulate the dependency of I(S;R) and II(S;R;E) on the number
% of trials

clear; close all;

MINT_path = '\MINT\'; % MINT toolbox path
addpath(genpath(MINT_path)); 

rng('default') % set random number generator to ensure reproducibility

% Set simulation parameters
params.num_trials = [20,50,100,200,500,750,1000];
params.repetitions = 500;
params.sigma_encoding = 2; % encoding gaussian noise sigma
params.sigma_readout = 2; % readout gaussian noise sigma
params.exp_trials = 200; % number of experimental trials

% Set info calculation parameters
opts.suppressWarnings = 'True';
opts.bin_method = {'eqpop'};
opts.n_bins = {3};
opts.bias = 'shuffSub'; % shuffle-subtracted bias correction
opts.shuff = 100; % number of shufflings for the shuffle subtraction bias corr

% info_val structures initialization, for memory efficiency
tmp = zeros(numel(params.num_trials),params.repetitions);
info_val.iISE_u = tmp; info_val.iIRS_u = tmp; info_val.iIRC_u = tmp; info_val.iII_u = tmp;
info_val.iISE_b = tmp; info_val.iIRS_b = tmp; info_val.iIRC_b = tmp; info_val.iII_b = tmp;
clear tmp 

% Run simulations
tic;
for rep = 1:params.repetitions % loop over params.repetitions
    if mod(rep,10)==0
        rep
    end
    for i = 1:numel(params.num_trials) % loop over number of trials
    
        params.numdata=params.num_trials(i);
        [tmp_info] = simulate_and_compute_info(params,opts);

        % unbiased quantities
        info_val.iISE_u(i,rep) = tmp_info.iISE_u;
        info_val.iIRS_u(i,rep) = tmp_info.iIRS_u;
        info_val.iIRC_u(i,rep) = tmp_info.iIRC_u;
        info_val.iII_u(i,rep) = tmp_info.iII_u;
    
        % biased quantities
        info_val.iISE_b(i,rep) = tmp_info.iISE_b;
        info_val.iIRS_b(i,rep) = tmp_info.iIRS_b;
        info_val.iIRC_b(i,rep) = tmp_info.iIRC_b;
        info_val.iII_b(i,rep) = tmp_info.iII_b;

    end
end
elapsed_time = toc;
disp(['Simulation ran in ',num2str(elapsed_time,2), ' sec'])

fname = ['sim_params.num_trials_II_',num2str(params.repetitions),'reps_sigma',num2str(params.sigma_encoding),'.m'];
save(fname,'info_val','params','opts')

%% Plot results
cols = [0.7,0,0;0.3,0.3,0.3]; % unbiased and biased lines colors

fig=figure('Position',[264,448,744,209]);

% I(S;R) plot
ax1 = subplot(1,3,1);
hold on
h(1)=linePlt(info_val.iIRS_u,params.num_trials,cols(1,:)); % unbiased line plot
h(2)=linePlt(info_val.iIRS_b,params.num_trials,cols(2,:)); % baised line plot
ax1=adjust_ax_scale(ax1,mean(info_val.iIRS_b,2)); % adjust y axis range
xline(params.exp_trials,'--') % add exp trials x axis line
yline(mean(info_val.iIRS_b(end,:),2),'b--') % add 'asymptotic' info y line
title('I(S;R)')
xlabel('# trials')
ylabel('info [bits]')

% I(R;E) plot
ax2=subplot(1,3,2);
hold on

% plot unbiased quantity
linePlt(info_val.iIRC_u,params.num_trials,cols(1,:)); % unbiased line plot
linePlt(info_val.iIRC_b,params.num_trials,cols(2,:)); % biased line plot
ax2.YLim = ax1.YLim; % adjust y axis range (match I(S;R) axis)
xline(params.exp_trials,'--') % add exp trials x axis line
yline(mean(info_val.iIRC_b(end,:),2),'b--') % add 'asymptotic' info y line
title('I(R;E)')
xlabel('# trials')
ylabel('info [bits]')

% II plot
ax=subplot(1,3,3);
hold on
linePlt(info_val.iII_u,params.num_trials,cols(1,:)); % unbiased line plot
linePlt(info_val.iII_b,params.num_trials,cols(2,:)); % baised line plot
ax=adjust_ax_scale(ax,mean(info_val.iII_b,2)); % adjust y axis range
xline(params.exp_trials,'--') % add exp trials x axis line
yline(mean(info_val.iII_b(end,:),2),'b--') % add 'asymptotic' info y line
title('II(S;R;E)')
xlabel('# trials')
ylabel('info [bits]')
legend([h(1),h(2)],{'shuffle-subtract','plugin'})

% Save figure
fname = ['info_wTrials_',num2str(params.repetitions),'reps_sigma',num2str(params.sigma_encoding)];
saveas(fig,[fname,'.png'])
saveas(fig,[fname,'.svg'])

%% helper functions

function [info_val] = simulate_and_compute_info(params,opts)
% Simulate data
numdata=params.numdata;
sigma_encoding=params.sigma_encoding;
sigma_readout=params.sigma_readout;

r=randn(1,numdata);
sample= r + sigma_encoding*randn(1,numdata);
est= r + sigma_readout*randn(1,numdata);

S=sample; 
C=est;
R=r;

% Compute info values
[iISE_u,iISE_b] = MI({S, C}, {'I(A;B)'}, opts);
[iII_u, iII_b] = II({S, R, C}, {'II(A,B,C)'}, opts);
[iIRS_u, iIRS_b] = MI({S, R}, {'I(A;B)'}, opts);
[iIRC_u, iIRC_b] = MI({C, R}, {'I(A;B)'}, opts);

% unbiased quantities
info_val.iISE_u = iISE_u{1};
info_val.iII_u = iII_u{1};
info_val.iIRS_u = iIRS_u{1};
info_val.iIRC_u = iIRC_u{1};

% biased quantities
info_val.iISE_b = iISE_b{1};
info_val.iII_b = iII_b{1};
info_val.iIRS_b = iIRS_b{1};
info_val.iIRC_b = iIRC_b{1};

end

% Plotting functions
function ax = adjust_ax_scale(ax,data)
% Adjsut Y axis range
yAxis_range = ax.YLim(2)-data(end);
ax.YLim(1) = data(end)-0.5*yAxis_range;
end

function h=linePlt(data_plt,x_vals,col)
mean_plt = mean(data_plt,2);
sem_plt = std(data_plt,[],2)/sqrt(size(data_plt,2));
h=plot(x_vals,mean_plt,'LineWidth',2,'Color',col);
shadedErrorBar(x_vals,mean_plt,sem_plt,'lineProps',{'color',col})
end
