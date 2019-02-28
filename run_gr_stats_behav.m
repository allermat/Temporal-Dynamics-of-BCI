expStage = 'final';
dataFileName = 'BEHAV_ANAL_group.mat';
load(fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),dataFileName),...
    'grWav','grVisLocCorr','grAudLocCorr');

% Statistics for Wav
indivData = table2struct(grWav.data(grWav.data.groupings == 'all',3:end))';
grStatsWav = mvpa.compAVmodelWeightStats(indivData,1,'mode','behav');

% Statistics for visual localization performance (comparing performance
% between reliability levels)
rhoVlow = atanh(grVisLocCorr.data.rho(grVisLocCorr.data.relV == 12));
rhoVhigh = atanh(grVisLocCorr.data.rho(grVisLocCorr.data.relV == 2));
rhoA = atanh(grAudLocCorr.data.rho);

rhoDiff_Vhigh_A = rhoVhigh-rhoA;
rhoDiff_Vlow_A = rhoVlow-rhoA;

grStatsVisCorr = permttest(rhoDiff_Vlow_A,0);
grStatsVisCorr = cat(1,grStatsVisCorr,permttest(rhoDiff_Vhigh_A,0));
grStatsVisCorr.Properties.RowNames = {'Vlow_vs_A','Vhigh_vs_A'};
save(fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),dataFileName),...
     '-append','grStatsWav','grStatsVisCorr');
 
% % Statistics for psalpha high vs low
% varNames = grWav.data.Properties.VariableNames;
% varNames = varNames(~cellfun(@isempty,regexp(varNames,'_wav$','once')));
% % Getting psalpha low and high data separately
% tempTable = grWav.data(~cellfun(@isempty,regexp(cellstr(grWav.data.groupings),...
%                                                 'psalpha_low','once')),:);
% temp_low = varfun(@(x) x,tempTable,'InputVariables',varNames,...
%                   'GroupingVariables',{'subID'},'OutputFormat','cell');
% tempTable = grWav.data(~cellfun(@isempty,regexp(cellstr(grWav.data.groupings),...
%                                                 'psalpha_high','once')),:);
% temp_high = varfun(@(x) x,tempTable,'InputVariables',varNames,...
%                   'GroupingVariables',{'subID'},'OutputFormat','cell');
% nSubj = size(temp_low,1);
% nTimePoints = size(temp_low{1,1},1);
% % Averaging over the conditions of the original 2x2x2 design
% temp_low = mat2cell(mean(cell2mat(temp_low),2),ones(nSubj,1)*nTimePoints,1);
% temp_high = mat2cell(mean(cell2mat(temp_high),2),ones(nSubj,1)*nTimePoints,1);

% indivData = cell2struct(cat(2,temp_low,temp_high),{'psalow_wav','psahigh_wav'},2)';

% grStatsPsalpha = mvpa.compAVmodelWeightStats(indivData,1:nTimePoints,...
%                                              'mode','behav',...
%                                              'expDesign','1way');

% save(fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),dataFileName),...
%      '-append','grStatsPsalpha');

% indivWavs = grWav.data(grWav.data.psalpha == 'all',:);
% indivWavsMain = indivWavs(:,1:2);
% indivWavsMain.('r_wav') = circ_mean([indivWavs.('r_d_a_wav'),indivWavs.('r_d_v_wav'),...
%                                      indivWavs.('r_D_a_wav'),indivWavs.('r_D_v_wav')],[],2);
% indivWavsMain.('R_wav') = circ_mean([indivWavs.('R_d_a_wav'),indivWavs.('R_d_v_wav'),...
%                                      indivWavs.('R_D_a_wav'),indivWavs.('R_D_v_wav')],[],2);
% indivWavsMain.('a_wav') = circ_mean([indivWavs.('r_d_a_wav'),indivWavs.('R_d_a_wav'),...
%                                      indivWavs.('r_D_a_wav'),indivWavs.('R_D_a_wav')],[],2);
% indivWavsMain.('v_wav') = circ_mean([indivWavs.('r_d_v_wav'),indivWavs.('R_d_v_wav'),...
%                                      indivWavs.('r_D_v_wav'),indivWavs.('R_D_v_wav')],[],2);
% indivWavsMain.('d_wav') = circ_mean([indivWavs.('r_d_a_wav'),indivWavs.('r_d_v_wav'),...
%                                      indivWavs.('R_D_a_wav'),indivWavs.('R_d_v_wav')],[],2);
% indivWavsMain.('D_wav') = circ_mean([indivWavs.('r_D_a_wav'),indivWavs.('r_D_v_wav'),...
%                                      indivWavs.('R_D_a_wav'),indivWavs.('R_D_v_wav')],[],2);

% indivWavs2way = indivWavs(:,1:2);
% indivWavs2way.('r_d_wav') = circ_mean([indivWavs.('r_d_a_wav'),indivWavs.('r_d_v_wav')],[],2);
% indivWavs2way.('r_D_wav') = circ_mean([indivWavs.('r_D_a_wav'),indivWavs.('r_D_v_wav')],[],2);
% indivWavs2way.('R_d_wav') = circ_mean([indivWavs.('R_d_a_wav'),indivWavs.('R_d_v_wav')],[],2);
% indivWavs2way.('R_D_wav') = circ_mean([indivWavs.('R_D_a_wav'),indivWavs.('R_D_v_wav')],[],2);
% indivWavs2way.('r_a_wav') = circ_mean([indivWavs.('r_d_a_wav'),indivWavs.('r_D_a_wav')],[],2);
% indivWavs2way.('r_v_wav') = circ_mean([indivWavs.('r_d_v_wav'),indivWavs.('r_D_v_wav')],[],2);
% indivWavs2way.('R_a_wav') = circ_mean([indivWavs.('R_d_a_wav'),indivWavs.('R_D_a_wav')],[],2);
% indivWavs2way.('R_v_wav') = circ_mean([indivWavs.('R_d_v_wav'),indivWavs.('R_D_v_wav')],[],2);
% indivWavs2way.('d_a_wav') = circ_mean([indivWavs.('r_d_a_wav'),indivWavs.('R_d_a_wav')],[],2);
% indivWavs2way.('d_v_wav') = circ_mean([indivWavs.('r_d_v_wav'),indivWavs.('R_d_v_wav')],[],2);
% indivWavs2way.('D_a_wav') = circ_mean([indivWavs.('r_D_a_wav'),indivWavs.('R_D_a_wav')],[],2);
% indivWavs2way.('D_v_wav') = circ_mean([indivWavs.('r_D_v_wav'),indivWavs.('R_D_v_wav')],[],2);
% nSubj = size(indivWavs2way,1);

% % Main effects
% analNames = {'VR','Task','Disp'};
% analVarNames = {{'r_wav','R_wav'},{'a_wav','v_wav'},{'d_wav','D_wav'}};
% for i = 1:numel(analNames)
    
%     % Arrays for data collection
%     grStats.(analNames{i}) = table;
%     varNames = analVarNames{i};
    
%     % Cycling through bins
%     alpha = table2array(indivWavsMain(:,varNames));
%     idx = repmat([1,2],nSubj,1);
%     idx = idx(:);

%     [p,t] = circ_wwtest(alpha,idx);
    
%     temp = table(p(1),{t},'VariableNames',[strcat({'p_'},analNames{i}),'ANOVA_table']);
%     grStats.(analNames{i}) = cat(1,grStats.(analNames{i}),temp);
% end

% % Two-way interactions
% analNames = {'VR_x_Disp','VR_x_Task','Disp_x_Task'};
% analFactorNames = {{'VR','Disp'},{'VR','Task'},{'Disp','Task'}};
% analVarNames = {{'r_d_wav','r_D_wav','R_d_wav','R_D_wav'},{'r_a_wav','r_v_wav','R_a_wav','R_v_wav'},...
%     {'d_a_wav','d_v_wav','D_a_wav','D_v_wav'}};

% for i = 1:numel(analNames)
    
%     % Arrays for data collection
%     grStats.(analNames{i}) = table;
%     factorNames = analFactorNames{i};
%     varNames = analVarNames{i};
    
%     % Cycling through bins
%     alpha = table2array(indivWavs2way(:,varNames));
%     % alpha = radians(alpha(:));
%     idp = repmat([1,1,2,2],nSubj,1);
%     idp = idp(:);
%     idq = repmat([1,2,1,2],nSubj,1);
%     idq = idq(:);
%     [p,t] = circ_hktest(alpha,idp,idq,1,factorNames);
    
%     temp = table(p(1),p(2),p(3),{t},'VariableNames',[strcat({'p_','p_'},factorNames),'p_Int','ANOVA_table']);
%     grStats.(analNames{i}) = cat(1,grStats.(analNames{i}),temp);
% end