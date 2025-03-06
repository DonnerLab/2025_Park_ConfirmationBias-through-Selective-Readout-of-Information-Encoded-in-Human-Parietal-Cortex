function plot_eight(out, conv, row, cols, ttxt1, ttxt2, axpars)

xlm = [-30 30]; ylm = [-30 30]; xlmb = [0.2 2.8]; ylmb = [0 0.1];

ms = 3; mcol = [170 170 170]./255;
%%% corr(S,E), I(S;E)
axes(gcf, 'Position', [axpars.x_st(1), axpars.y_st(row), axpars.x_w(1), axpars.y_h(1)]); hold on;
% title([ttxt1{1} '=' num2str(out.cISE, '%1.3f')],'Interpreter', 'Latex')
% title(['r=' num2str(out.cISE, '%1.3f')])
plot(out.sample,out.est,'ok', 'markersize', ms, 'markerfacecolor', mcol, 'markeredgecolor', 'none')
% xlabel('S','Interpreter', 'Latex'); 
xlim(xlm)
set(gca, 'xticklabels', [], 'yticklabels', [],'tickdir', 'out')
% xlh.Position = [0 -0.1 -1];
% ylabel('E','Interpreter', 'Latex'); 
ylim(ylm)
axis square; box off;

axes('Position', [axpars.x_st(2), axpars.y_st(row), axpars.x_w(2), axpars.y_h(2)]); hold on;
% title([ttxt2{1} '=' num2str(out.dSC, '%1.3f')],'Interpreter', 'Latex')
% title(['I=' num2str(out.dSC, '%1.3f')])
% title([ttxt2{1}],'Interpreter', 'Latex')
bh = bar([1 2], [out.dSC conv.dSC-out.dSC]);
bh.FaceColor = 'flat';
bh.CData = cols;
bh.LineStyle = 'none';
bh.BarWidth = 0.8;
ax = gca;
set(ax, 'xtick', [1 2], 'xticklabels', conv.xtlb, 'ytick', 10^-1.*[0 1], 'yticklabels', string(10^-1.*[0 2]),'tickdir', 'out')
ax.YAxis.Exponent = -1;
xlim(xlmb);ylim(ylmb)
axis square; box off;


%%% corr(S,Shat), I(S;R)
axes('Position', [axpars.x_st(3), axpars.y_st(row), axpars.x_w(3), axpars.y_h(3)]); hold on;
% title([ttxt1{2} '=' num2str(out.cIS, '%1.3f')],'Interpreter', 'Latex')
% title(['r=' num2str(out.cIS, '%1.3f')])
plot(out.sample,out.sample_dec,'ok', 'markersize', ms, 'markerfacecolor', mcol, 'markeredgecolor', 'none')
set(gca, 'xticklabels', [], 'yticklabels', [],'tickdir', 'out')
% xlabel('S','Interpreter', 'Latex'); 
xlim(xlm)
% ylabel('$$\hat{S}$$', 'Interpreter', 'Latex'); 
ylim(ylm)
axis square; box off;

axes('Position', [axpars.x_st(4), axpars.y_st(row), axpars.x_w(4), axpars.y_h(4)]); hold on;
% title([ttxt2{2} '=' num2str(out.dIRS, '%1.3f')],'Interpreter', 'Latex')
% title(['I=' num2str(out.dIRS, '%1.3f')])
% title([ttxt2{2}],'Interpreter', 'Latex')
bh = bar([1 2], [out.dIRS conv.dIRS-out.dIRS]);
bh.FaceColor = 'flat';
bh.CData = cols;
bh.LineStyle = 'none';
bh.BarWidth = 0.8;
ax = gca;
set(ax, 'xtick', [1 2], 'xticklabels', conv.xtlb, 'ytick', 10^-1.*[0 1], 'yticklabels', string(10^-1.*[0 2]),'tickdir', 'out')
ax.YAxis.Exponent = -1;
xlim(xlmb);ylim(ylmb)
axis square; box off;


%%% corr(E,Ehat), I(R;E)
axes('Position', [axpars.x_st(5), axpars.y_st(row), axpars.x_w(5), axpars.y_h(5)]); hold on;
% title([ttxt1{3} '=' num2str(out.cIC, '%1.3f')],'Interpreter', 'Latex')
% title(['r=' num2str(out.cIC, '%1.3f')])
plot(out.est,out.est_dec,'ok', 'markersize', ms, 'markerfacecolor', mcol, 'markeredgecolor', 'none')
set(gca, 'xticklabels', [], 'yticklabels', [],'tickdir', 'out')
% xlabel('E','Interpreter', 'Latex'); 
xlim(xlm)
% ylabel('$$\hat{E}$$', 'Interpreter', 'Latex'); 
ylim(ylm)
axis square; box off;

axes('Position', [axpars.x_st(6), axpars.y_st(row), axpars.x_w(6), axpars.y_h(6)]); hold on;
% title([ttxt2{3} '=' num2str(out.dIRC, '%1.3f')],'Interpreter', 'Latex')
% title(['I=' num2str(out.dIRC, '%1.3f')])
% title([ttxt2{3}],'Interpreter', 'Latex')
bh = bar([1 2], [out.dIRC conv.dIRC-out.dIRC]);
bh.FaceColor = 'flat';
bh.CData = cols;
bh.LineStyle = 'none';
bh.BarWidth = 0.8;
ax = gca;
set(ax, 'xtick', [1 2], 'xticklabels', conv.xtlb, 'ytick', 10^-1.*[0 1], 'yticklabels', string(10^-1.*[0 2]),'tickdir', 'out')
ax.YAxis.Exponent = -1;
xlim(xlmb);ylim(ylmb)
axis square; box off;


%%% corr(Shat, E), II(S;R;E)
axes('Position', [axpars.x_st(7), axpars.y_st(row), axpars.x_w(7), axpars.y_h(7)]); hold on;
% title([ttxt1{4} '=' num2str(out.cII, '%1.3f')],'Interpreter', 'Latex')
% title(['r=' num2str(out.cII, '%1.3f')])
plot(out.sample_dec,out.est,'ok', 'markersize', ms, 'markerfacecolor', mcol, 'markeredgecolor', 'none')
set(gca, 'xticklabels', [], 'yticklabels', [],'tickdir', 'out')
% xlabel('$$\hat{S}$$', 'Interpreter', 'Latex'); 
xlim(xlm)
% ylabel('E','Interpreter', 'Latex'); 
ylim(ylm)
axis square; box off;

axes('Position', [axpars.x_st(8), axpars.y_st(row), axpars.x_w(8), axpars.y_h(8)]); hold on;
% title([ttxt2{4} '=' num2str(out.dII, '%1.3f')],'Interpreter', 'Latex')
% title(['I=' num2str(out.dII, '%1.3f')])
% title([ttxt2{4}],'Interpreter', 'Latex')
bh = bar([1 2], [out.dII conv.dII-out.dII]);
bh.FaceColor = 'flat';
bh.CData = cols;
bh.LineStyle = 'none';
bh.BarWidth = 0.8;
ax = gca;
set(ax, 'xtick', [1 2], 'xticklabels', conv.xtlb, 'ytick', 10^-2.*[0 5],'yticklabels', string(10^-2.*[0 5]),'tickdir', 'out')
ax.YAxis.Exponent = -2;
xlim(xlmb);ylim([0 0.05])
axis square; box off;
end