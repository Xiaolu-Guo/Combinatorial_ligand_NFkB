figure_save_name = 'nocompete';


load(strcat(data_save_file_path,'W_diff_mat_',figure_save_name,'.mat'))

w_diff_mat_non_compete = W_diff_mat(1:15,1:15);
CpG_wids_no_compete = w_diff_mat_non_compete(13,[1:12,14:15]);
CpG_wids_no_compete_31cond = W_diff_mat(13,[1:12,14:31]);
% idx = triu(true(size(w_diff_mat_non_compete)), 1);    % upper triangle, k=1 removes diagonal
% upper_no_diag = w_diff_mat_non_compete(idx);

figure_save_name = 'compete';

load(strcat(data_save_file_path,'W_diff_mat_',figure_save_name,'.mat'))

w_diff_mat_compete = W_diff_mat(1:15,1:15);
CpG_wids_compete = w_diff_mat_compete(13,[1:12,14:15]);
CpG_wids_compete_31cond = W_diff_mat(13,[1:12,14:31]);


%% plot overlaping 
if 1 % plot hist for 15 conditions (only single & dual)
figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.5;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

histogram(CpG_wids_no_compete,0:0.25:1);hold on
histogram(CpG_wids_compete,0:0.25:1);hold on

if 1 % run Wilcoxon signed-rank tests
    [p, ~, ~] = signrank(CpG_wids_no_compete, CpG_wids_compete,'tail','left')
end

legend('non-compete','compete')

YL = ylim();
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_vs_non_compete_withlabel'),'pdf')

set(gca,'XTickLabel',{},'YTickLabel',{})
legend off
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_vs_non_compete'),'pdf')
close()

end

%% plot histogram seperately 

if 1 % plot hist for 15 conditions (only single & dual)
figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.5;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

histogram(CpG_wids_no_compete,0:0.25:1);hold on
% histogram(CpG_wids_compete,0:0.25:1);hold on

ylim(YL)
saveas(gcf,strcat(fig_save_path,'hist_Wdis_non_compete_withlabel'),'pdf')
set(gca,'XTickLabel',{},'YTickLabel',{})
saveas(gcf,strcat(fig_save_path,'hist_Wdis_non_compete'),'pdf')
close()


figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.5;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

% histogram(CpG_wids_no_compete,0:0.25:1);hold on
histogram(CpG_wids_compete,0:0.25:1);hold on
ylim(YL)
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_withlabel'),'pdf')

set(gca,'XTickLabel',{},'YTickLabel',{})
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete'),'pdf')
close()

end

%% plot overlapping

if 1 % plot hist for 31 conditions, all comb

figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.5;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

histogram(CpG_wids_no_compete_31cond,0:0.25:1);hold on
histogram(CpG_wids_compete_31cond,0:0.25:1);hold on
legend('non-compete','compete')
YL2 = ylim();
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_vs_non_compete_31cond_withlabel'),'pdf')

set(gca,'XTickLabel',{},'YTickLabel',{})
legend off
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_vs_non_31cond_compete'),'pdf')
close()
end


%% plot seperately

if 1 % plot hist for 31 conditions, all comb

figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.5;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

histogram(CpG_wids_no_compete_31cond,0:0.25:1);hold on
% histogram(CpG_wids_compete_31cond,0:0.25:1);hold on

ylim(YL2)
saveas(gcf,strcat(fig_save_path,'hist_Wdis_non_compete_31cond_withlabel'),'pdf')

set(gca,'XTickLabel',{},'YTickLabel',{})
legend off
saveas(gcf,strcat(fig_save_path,'hist_Wdis_non_compete_31cond'),'pdf')
close()


figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.5;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

% histogram(CpG_wids_no_compete_31cond,0:0.25:1);hold on
histogram(CpG_wids_compete_31cond,0:0.25:1);hold on

ylim(YL2)
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_31cond_withlabel'),'pdf')

set(gca,'XTickLabel',{},'YTickLabel',{})
legend off
saveas(gcf,strcat(fig_save_path,'hist_Wdis_compete_31cond'),'pdf')
close()
end

