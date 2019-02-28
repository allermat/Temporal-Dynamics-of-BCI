% Load data
dataFname = 'bci_model_sel_MVPA_group_AV-c-av_5mdl_180619130318.mat';
load(fullfile(DEC_2_setupdir('final','anal_eeg_group_mvpa'),'sample-wise-sm-avg',...
dataFname));
% colors = {[44,123,182]./255,[253,174,97]./255,[171,217,233]./255,...
%           [215,25,28]./255,[254,224,144]./255};
% colors = {[69,117,180]./255,[254,224,144]./255,[145,191,219]./255,...
          % [215,48,39]./255,[224,243,248]./255};
% colors = {[33,102,172]./255,[253,219,199]./255,[103,169,207]./255,...
          % [178,24,43]./255,[239,138,98]./255};
% colors = {[33,102,172]./255,[239,138,98]./255,[103,169,207]./255,...
          % [178,24,43]./255,[253,219,199]./255};
% Colors I find best
% colors = {[69,117,180]./255,[252,141,89]./255,[145,191,219]./255,...
%           [215,48,39]./255,[254,224,144]./255};
% Colors recommended by Uta
colors = {[34,114,179]./255,[247,182,49]./255,[68,176,82]./255,...
          [191,38,38]./255,[145,191,219]./255};
colorM = cell2mat(colors');
modelNames = models(ismember(models,modelSel));

f = figure();
subplot(5,1,1);
h = area(timePoints,rfx.xp); 
ylim([0,1]);
xlim(timePoints([1,end]));
legend(modelNames,'Location','NorthEastOutside');
ylabel('p_c_u_m_u_l_a_t_i_v_e');
title('Exceedance probability');
subplot(5,1,2);
h = area(timePoints,rfx.pxp);
ylim([0,1]);
xlim(timePoints([1,end]));
legend(modelNames,'Location','NorthEastOutside');
ylabel('p_c_u_m_u_l_a_t_i_v_e');
title('Protected exceedance probability');
% This below should have been rfx.exp_r for posterior probability,
% but I messed it up in the initial plotting, this should be
% re-plotted
subplot(5,1,3);
h = area(timePoints,rfx.exp_r);
ylim([0,1]);
xlim(timePoints([1,end]));
legend(modelNames,'Location','NorthEastOutside');
ylabel('p_c_u_m_u_l_a_t_i_v_e');
title('Posterior probability');
colormap(colorM(ismember(models,modelSel),:));

% subplot(6,1,4);
% for i = 1:numel(modelNames)
%     idx = find(ismember(models,modelSel));
%     herr(i) = shadedErrorBar(timePoints,rfx.logLikeMeans(:,i),rfx.logLikeSEMs(:,i),...
%                              {'color',colorM(idx(i),:),'lineWidth',1.5},true);%#ok
%     if i == 1, hold on; end
% end
% legend([herr.mainLine],modelNames,'Location','NorthEastOutside');
% ylabel('log-Likelihood');
% title('log-Likelihood');

subplot(5,1,4);
for i = 1:numel(modelNames)
    idx = find(ismember(models,modelSel));
    herr(i) = shadedErrorBar(timePoints,rfx.R2Means(:,i),rfx.R2SEMs(:,i),...
                             {'color',colorM(idx(i),:),'lineWidth',1.5},true);
    if i == 1, hold on; end
        
end
legend([herr.mainLine],modelNames,'Location','NorthEastOutside');
ylabel('R^2');
title('Coefficient of determination');

s = subplot(5,1,5);
herr = shadedErrorBar(timePoints,rfx.plVarMean,rfx.plVarSEM,...
                         {'color','k','lineWidth',1.5},true);
ylabel('Variance');
title('Variance of decoded labels');