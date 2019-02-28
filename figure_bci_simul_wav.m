% Load data
load(fullfile(DEC_2_setupdir('final','anal_behav_group'),'bci_simul_wav_BEHAV_group.mat'));
% Converting radians to degrees
grWav.meanWav{:,3:end} = arrayfun(@degrees,grWav.meanWav{:,3:end});

for i = 1:size(grWav.meanWav,1)
    figure();
    data = grWav.meanWav(i,:);
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
    title(sprintf('Audio-visual weights predicted by the\n%s model',...
                  char(grWav.meanWav.model(i))));
    ylabel('W_A_V (degree)');
    set(gca,'XTick',[1,3],'XTickLabel',{'D-','D+'});
end