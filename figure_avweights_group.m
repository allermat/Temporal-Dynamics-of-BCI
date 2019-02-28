function figure_avweights_group(grWav)
% All
% Converting radians to degrees
grWav.meanWav{:,3:end} = arrayfun(@degrees,grWav.meanWav{:,3:end});

figure();
data = grWav.meanWav(grWav.meanWav.groupings == 'all',:);
errorbar([1,3]',[data.R_d_v_wav,data.R_D_v_wav]',[data.R_d_v_ci,data.R_D_v_ci]',...
    'Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'LineWidth',1.5,'LineStyle','--'); hold('on');
errorbar([1,3]',[data.r_d_v_wav,data.r_D_v_wav]',[data.r_d_v_ci,data.r_D_v_ci]',...
    'Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'LineWidth',1.5,'LineStyle','-');
errorbar([1,3]',[data.R_d_a_wav,data.R_D_a_wav]',[data.R_d_a_ci,data.R_D_a_ci]',...
    'Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','--');
errorbar([1,3]',[data.r_d_a_wav,data.r_D_a_wav]',[data.r_d_a_ci,data.r_D_a_ci]',...
    'Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','-');
line([[0,4]',[0,4]'],[[90,90]',[0,0]'],'Color',[0,0,0]);
xlim([0,4]);
ylim([-10,100]);
legend('V, VR+','V, VR-','A, VR+','A, VR-','Location','NorthEastOutside');
title('Group Level Relative audio-visual weight (W_A_V)');
ylabel('W_A_V (degree)');
set(gca,'XTick',[1,3],'XTickLabel',{'D-','D+'});

% High vs low alpha
figure();
dataLow = grWav.meanWav(grWav.meanWav.groupings == 'psalpha_low',:);
dataHigh = grWav.meanWav(grWav.meanWav.groupings == 'psalpha_high',:);
xLoc = [1,3];
offset = [0.1,0.1];
errorbar((xLoc - offset)',[dataLow.R_d_v_wav,dataLow.R_D_v_wav]',[dataLow.R_d_v_ci,dataLow.R_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'LineWidth',1.5,'LineStyle','--'); hold('on');
errorbar((xLoc + offset)',[dataHigh.R_d_v_wav,dataHigh.R_D_v_wav]',[dataHigh.R_d_v_ci,dataHigh.R_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'LineWidth',1.5,'LineStyle','--');
errorbar((xLoc - offset)',[dataLow.r_d_v_wav,dataLow.r_D_v_wav]',[dataLow.r_d_v_ci,dataLow.r_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'LineWidth',1.5,'LineStyle','-');
errorbar((xLoc + offset)',[dataHigh.r_d_v_wav,dataHigh.r_D_v_wav]',[dataHigh.r_d_v_ci,dataHigh.r_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'LineWidth',1.5,'LineStyle','-');
errorbar((xLoc - offset)',[dataLow.R_d_a_wav,dataLow.R_D_a_wav]',[dataLow.R_d_a_ci,dataLow.R_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','--');
errorbar((xLoc + offset)',[dataHigh.R_d_a_wav,dataHigh.R_D_a_wav]',[dataHigh.R_d_a_ci,dataHigh.R_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','--');
errorbar((xLoc - offset)',[dataLow.r_d_a_wav,dataLow.r_D_a_wav]',[dataLow.r_d_a_ci,dataLow.r_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','-');
errorbar((xLoc + offset)',[dataHigh.r_d_a_wav,dataHigh.r_D_a_wav]',[dataHigh.r_d_a_ci,dataHigh.r_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','-');

line([[0,4]',[0,4]'],[[90,90]',[0,0]'],'Color',[0,0,0]);
xlim([0,4]);
% ylim([-10,100]);
legend('V, VR+, alpha-','V, VR+, alpha+','V, VR-, alpha-','V, VR-, alpha+',...
       'A, VR+, alpha-','A, VR+, alpha+','A, VR-, alpha-','A, VR-, alpha+','Location','NorthEastOutside');
title(sprintf(['Group Level Relative audio-visual weight (W_A_V)',...
      '\nhigh vs low prestim alpha power']));
ylabel('W_A_V (degree)');
set(gca,'XTick',[1,3],'XTickLabel',{'D-','D+'});

% High vs low cmb
figure();
dataLow = grWav.meanWav(grWav.meanWav.groupings == 'cmb_low',:);
dataHigh = grWav.meanWav(grWav.meanWav.groupings == 'cmb_high',:);
xLoc = [1,3];
offset = [0.1,0.1];
errorbar((xLoc - offset)',[dataLow.R_d_v_wav,dataLow.R_D_v_wav]',[dataLow.R_d_v_ci,dataLow.R_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'LineWidth',1.5,'LineStyle','--'); hold('on');
errorbar((xLoc + offset)',[dataHigh.R_d_v_wav,dataHigh.R_D_v_wav]',[dataHigh.R_d_v_ci,dataHigh.R_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'LineWidth',1.5,'LineStyle','--');
errorbar((xLoc - offset)',[dataLow.r_d_v_wav,dataLow.r_D_v_wav]',[dataLow.r_d_v_ci,dataLow.r_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'LineWidth',1.5,'LineStyle','-');
errorbar((xLoc + offset)',[dataHigh.r_d_v_wav,dataHigh.r_D_v_wav]',[dataHigh.r_d_v_ci,dataHigh.r_D_v_ci]',...
    'o','Color',[0,0,0],'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'LineWidth',1.5,'LineStyle','-');
errorbar((xLoc - offset)',[dataLow.R_d_a_wav,dataLow.R_D_a_wav]',[dataLow.R_d_a_ci,dataLow.R_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','--');
errorbar((xLoc + offset)',[dataHigh.R_d_a_wav,dataHigh.R_D_a_wav]',[dataHigh.R_d_a_ci,dataHigh.R_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','--');
errorbar((xLoc - offset)',[dataLow.r_d_a_wav,dataLow.r_D_a_wav]',[dataLow.r_d_a_ci,dataLow.r_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','-');
errorbar((xLoc + offset)',[dataHigh.r_d_a_wav,dataHigh.r_D_a_wav]',[dataHigh.r_d_a_ci,dataHigh.r_D_a_ci]',...
    'o','Color',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','-');

line([[0,4]',[0,4]'],[[90,90]',[0,0]'],'Color',[0,0,0]);
xlim([0,4]);
% ylim([-10,100]);
legend('V, VR+, cmb-','V, VR+, cmb+','V, VR-, cmb-','V, VR-, cmb+',...
       'A, VR+, cmb-','A, VR+, cmb+','A, VR-, cmb-','A, VR-, cmb+','Location','NorthEastOutside');
title(sprintf(['Group Level Relative audio-visual weight (W_A_V)',...
      '\nhigh vs low CMB']));
ylabel('W_A_V (degree)');
set(gca,'XTick',[1,3],'XTickLabel',{'D-','D+'});


end

