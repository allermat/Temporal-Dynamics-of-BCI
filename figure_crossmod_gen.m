function s = figure_crossmod_gen(subID,trainMethod)

expStage = 'final';

s = struct();

if strcmp(subID,'group')
    dir_analysis = fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),trainMethod);
    matchStr = 'gr_er_tr-[AV]_stim_gen-[AV].mat';
else
    dir_analysis = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),trainMethod);
    matchStr = ['er_tr-','[AV]','_s_-100-5-1000_gen-.*'];
end

saveDf = cd(dir_analysis);

fileList = cellstr(ls('*.mat'));
matchIdx = ~cellfun(@isempty,regexp(fileList,matchStr));
fileList = fileList(matchIdx);

in = struct;
in.f = figure();
h = {};
figOrder = [3,4,1,2];
for i = 1:size(fileList,1)
    
    switch i
      case 1
        condFname = 'trA_genA';
        title = 'Train A, generalize A';
      case 2
        condFname = 'trA_genV';
        title = 'Train A, generalize V';
      case 3
        condFname = 'trV_genA';
        title = 'Train V, generalize A';
      case 4
        condFname = 'trV_genV';
        title = 'Train V, generalize V';
    end
    
    s.(condFname) = mvpares(fileList{i});
    
    h{i} = subplot(2,2,figOrder(i));
    in.parentAxes = h{i};
    in.title = title;
    s.(condFname).setPlotTimeWin([-0.1,0.7]);
    showGenPerf(s.(condFname),'genTime','tr_x_tr','perfIdx','rf',...
        'cvEval','poolFolds','in',in);
    % axis('square');
    c(i,:) = caxis;
end

% Adding contours for significant bidirectional cross-modal
% generalization
stats = s.trA_genV.getStats('genPerf','genTime','tr_x_tr');
h_AV = stats.h_rf_poolFolds;
stats = s.trV_genA.getStats('genPerf','genTime','tr_x_tr');
h_VA = stats.h_rf_poolFolds;
contour(h{2},h_AV & h_VA,[1,1],'m');
contour(h{3},h_AV & h_VA,[1,1],'m');

if strcmp(subID,'group')
    suplabel(sprintf('Group level'),'t');
else
   suplabel(sprintf('Subject %s',subID),'t'); 
end

% cMaxAbs = max(abs(c(:)));
cMax = max(c(:));
cMin = min(c(:));
for i = 1:numel(h), caxis(h{i},[cMin,cMax]), end;

% Setting figure size
set(in.f,'Units','normalized','Position',[1/6,1/6,2/3,2/3],...
         'PaperPositionMode','auto');

cd(saveDf);