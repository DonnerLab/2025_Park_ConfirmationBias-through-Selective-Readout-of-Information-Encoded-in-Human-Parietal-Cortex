
addpath(genpath('MINT'))


% Simulation 1
% Consistent sample - same encoded info and good readout

% Simulation 2
% Inconsistent sample - same encoded info and good readout

% Simulation 3
% Inconsistent sample -  bad encoded info and good readout

% Simulation 4
% Inconsistent sample - same encoded info and bad strength of readout with good alignment

% Simulation 5
% Inconsistent sample - same encoded info and readout with same strength but bad alignement 


numSimul = 5;
outAll = cell(numSimul, 1);

numdata = 10000;
eNoise = [1 1 2 1 1]; % 1: same encoded info; 2: bad encoded info
dNoise = [1 1 1 2 1]; % 1: good readout; 2: bad readout
encoding_weights = [1 3; -1 -3; -1 -3; -1 -3; -1 -3]; % negative sign is negative weighting 
decoding_weights = [1 3; 1 3; 1 3; 1 3; 3 1]; % flipped value is bad alignment for readout


ttxt1 = {'$$corr(S,E)$$','$$corr(S,\hat{S})$$','$$corr(E,\hat{E})$$','$$corr(\hat{S},E)$$'}; % ylabel(ttxt1, 'Interpreter', 'Latex')
ttxt2 = {'I(S;E)','I(S;R)','I(R;E)','II(S;R;E)'}; 

% create axes
hw = 18; vw = 22;

bax = 3;
sax = 1;

lmar = 0.05;
rmar = 0.02;
marh1 = 0.05;
marh2 = 0.03;
axw = (1-lmar-rmar-3*marh1-4*marh2)/4; baxw = 0.75*axw; saxw = 0.25*axw;
marv = 0.07;
tmar = (1-numSimul*baxw-(numSimul-1)*marv)/2; 
bmar = tmar;

% axh = (vw-tmar-bmar-(numSimul-1)*marv)/numSimul;

x_st = []; 
x_st(1) = lmar;
x_st(2) = x_st(1)+baxw+marh2;
x_st(3) = x_st(2)+saxw+marh1;
x_st(4) = x_st(3)+baxw+marh2;
x_st(5) = x_st(4)+saxw+marh1;
x_st(6) = x_st(5)+baxw+marh2;
x_st(7) = x_st(6)+saxw+marh1;
x_st(8) = x_st(7)+baxw+marh2;

x_w = [];
x_w([1 3 5 7]) = baxw;
x_w([2 4 6 8]) = saxw;

y_st = [];
y_st(1) = 1-tmar-baxw;
for s = 2:numSimul
    y_st(s) = y_st(1)-(s-1)*(baxw+marv);
end
y_h = []; 
y_h([1 3 5 7]) = baxw; 
y_h([2 4 6 8]) = saxw;

axpars.x_st = x_st;
axpars.y_st = y_st;
axpars.x_w = x_w;
axpars.y_h = y_h;

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 hw vw];

% simulate and plot
conv = cell(numSimul, 1);
for s = 1:numSimul
    params.numdata = numdata;
    params.level_encoding = eNoise(s);
    params.level_readout = dNoise(s);
    params.beta_encoding = encoding_weights(s,:);
    params.beta_readout = decoding_weights(s,:);
    outAll{s} = compute_corr_and_info(params);
    if s == 1
        conv{s}.dSC=0;
        conv{s}.dIRS=0;
        conv{s}.dIRC=0;
        conv{s}.dII=0;
        cols = [10 184 0; 0 0 0]./255; % consistent
        conv{s}.xtlb = {'C',''};
    else
        conv{s} = outAll{1};
        cols = [255 184 10; 0 0 0]./255; % inconsistent
        conv{s}.xtlb = {'I','C-I'};
    end
    plot_eight(outAll{s}, conv{s}, s, cols, ttxt1, ttxt2, axpars);
end

vis = figure;
vis.Units = 'centimeters';
vis.Position = [2 2 hw vw];
% plot existing out structure
load('Simul_outAll.mat', 'outAll')
conv = cell(numSimul, 1);
for s = 1:numSimul
    if s == 1
        conv{s}.dSC=0;
        conv{s}.dIRS=0;
        conv{s}.dIRC=0;
        conv{s}.dII=0;
        cols = [10 184 0; 0 0 0]./255; % consistent
        conv{s}.xtlb = {'C',''};
    else
        conv{s} = outAll{1};
        cols = [255 184 10; 0 0 0]./255; % inconsistent
        conv{s}.xtlb = {'I','C-I'};
    end
    plot_eight(outAll{s}, conv{s}, s, cols, ttxt1, ttxt2, axpars);
end

sname = sprintf('simulation')
% exportgraphics(vis, [sname '.eps'],'ContentType','vector',...
%        'BackgroundColor','none')
% saveas(vis, sname, 'png')
   
   

function [out] = compute_corr_and_info(params)
numdata=params.numdata;
level_encoding=params.level_encoding;
level_readout = params.level_readout;
beta_encoding=params.beta_encoding;
sigma_encoding=level_encoding*norm(beta_encoding);
beta_readout=params.beta_readout;

sigma_readout=level_readout*norm(beta_readout);
r=randn(2,numdata);
sample=beta_encoding*r + sigma_encoding*randn(1,numdata);
est=beta_readout*r + sigma_readout*randn(1,numdata);
sample_dec=beta_encoding*r + sigma_encoding*randn(1,numdata);
est_dec=beta_readout*r + sigma_readout*randn(1,numdata);
for i=1:2
  sample_dec_indiv(i,:)=beta_encoding(i)*r(i,:) + sigma_encoding*randn(1,numdata);
  est_dec_indiv(i,:)=beta_readout(i)*r(i,:) + sigma_readout*randn(1,numdata);
end    

[out.cISE]=corr(sample',est'); % ISE=abs(ISE)
%[cIS]=corr(sample',sample_dec');
[cIS1]=corr(sample',sample_dec_indiv(1,:)');
[cIS2]=corr(sample',sample_dec_indiv(2,:)');
out.cIS=(cIS1+cIS2)/2;
[cIC1]=corr(est',est_dec_indiv(1,:)');
[cIC2]=corr(est',est_dec_indiv(2,:)');
out.cIC=(cIC1+cIC2)/2;
%[cIC]=corr(est',est_dec');
%[cII]=corr(sample_dec',est'); %II=abs(II)
[cII1]=corr(sample_dec_indiv(1,:)',est'); %II=abs(II)
[cII2]=corr(sample_dec_indiv(2,:)',est');
out.cII=(cII1+cII2)/2;

out.sample=sample;
out.sample_dec=sample_dec;
out.est=est;
out.est_dec=est_dec;

S=sample;
R1=r(1,:);
R2=r(2,:);
C=est;
opts.suppressWarnings='True';
opts.bin_method = {'eqpop'};
opts.n_bins = {2};
[dII1, PID_SR_C, PID_RC_S] = II({S, R1, C}, {'II(A,B,C)'}, opts);
[dII2, PID_SR_C, PID_RC_S] = II({S, R2, C}, {'II(A,B,C)'}, opts);
[dISC] = MI({S, C}, {'I(A;B)'}, opts);
[dIRS1] = MI({S, R1}, {'I(A;B)'}, opts);
[dIRS2] = MI({S, R2}, {'I(A;B)'}, opts);
[dIRC1] = MI({C, R1}, {'I(A;B)'}, opts);
[dIRC2] = MI({C, R2}, {'I(A;B)'}, opts);

out.dII=(dII1{1}+dII2{1})/2;
out.dSC = dISC{1};
out.dIRS=(dIRS1{1}+dIRS2{1})/2;
out.dIRC=(dIRC1{1}+dIRC2{1})/2;

end