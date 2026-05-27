function [] = bar_diff_exp_old_20251118(fig_save_path)
x_vec = [1,2,3];
y_vec_dur = [-0.3260,0.2303,-0.2921];
y_vec_total = [-0.4089,0.3194,-0.0839];

figure(1)
paperpos=[0,0,100,100]*0.6;
papersize=[100 100]*0.6;
draw_pos=[10,10,80,80]*0.6;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
%scatter(x_vec,y_vec_dur,15,'filled');hold on
bar(y_vec_dur,'BaseValue',-0.0);hold on
set(gca,'XTick',1:3,'XTickLabels',{},'YTick',-0.5:0.1:0.5,'YTickLabels',{})% '-0.5','0','0.5','1.0','1.5'

YL = ylim;
ylim([-0.5,0.5])
xlim([0,4])
saveas(gcf,strcat(fig_save_path,'bar_duration_diff_exp_old_20251118'),'pdf')
close()

figure(2)
paperpos=[0,0,100,100]*0.6;
papersize=[100 100]*0.6;
draw_pos=[10,10,80,80]*0.6;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
% scatter(x_vec,y_vec_total,15,'filled');hold on
% set(gca,'XTick',1:3,'XTickLabels',{},'YTickLabels',{})
bar(y_vec_total,'BaseValue',-0.0);hold on
set(gca,'XTick',1:3,'XTickLabels',{},'YTick',-0.5:0.1:0.5,'YTickLabels',{}) %'-0.5','0','0.5','1.0','1.5'

YL = ylim;
ylim([-0.5,0.5])
xlim([0,4])
saveas(gcf,strcat(fig_save_path,'bar_total_diff_exp_old_20251118'),'pdf')
close()

end
