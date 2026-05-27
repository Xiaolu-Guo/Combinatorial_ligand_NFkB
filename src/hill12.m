
save_fig_path_sub = '../subfigures2025/revision/figureS3/'; % should be update to figureS4
Kd = 0.001116;
Vmax = 18.79;
x = 0:0.000001: 0.005;

hill_coef = 1;
y1 = Vmax *x.^hill_coef./(Kd^hill_coef+x.^hill_coef);
hill_coef = 2;
y2 = Vmax *x.^hill_coef./(Kd^hill_coef+x.^hill_coef);

figure(1)
fold_zoom_1 = 0.8;
paperpos=[0,0,180,100];
papersize=[180 100];
draw_pos=[10,10,160,80];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)
plot(x,y1);hold on
plot(x,y2);hold on
xlim([0,0.004])

saveas(gcf,strcat(save_fig_path_sub,'hill12','_withlabel'),'pdf')
set(gca,'XTickLabel',{},'YTickLabel',{})
saveas(gcf,strcat(save_fig_path_sub,'hill12'),'pdf')
close