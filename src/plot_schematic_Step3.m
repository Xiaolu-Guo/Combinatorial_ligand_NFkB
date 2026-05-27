i_row = i_row_plot_cond+1;
i_col = i_col_plot_cond+1;
total_plot_column = 3;
i_ind = 1;

for i_species = 1:total_species

data_traj = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(i_species:12:end,:);
ymax = max([max(data_traj(plot_cell_index_cond(i_ind),:),[],'all'),...
    max(data_traj_D1(plot_cell_index_cond(i_ind),:),[],'all'),...
    max(data_traj_D2(plot_cell_index_cond(i_ind),:),[],'all')]);


subplot(total_species,total_plot_column,(i_species - 1) * total_plot_column + 1)
plot(0:5:480,data_traj_D1(plot_cell_index_cond(i_ind),:), 'r','LineWidth',2); hold on
xlim([0,480])
ylim([0,ymax*1.05])

subplot(total_species,total_plot_column,(i_species - 1) * total_plot_column + 2)
plot(0:5:480,data_traj_D2(plot_cell_index_cond(i_ind),:),'m','LineWidth',2); hold on
xlim([0,480])
ylim([0,ymax*1.05])

subplot(total_species,total_plot_column,(i_species - 1) * total_plot_column + 3)
plot(0:5:480,data_traj(plot_cell_index_cond(i_ind),:),'b','LineWidth',2);hold on

xlim([0,480])
ylim([0,ymax*1.05])
end