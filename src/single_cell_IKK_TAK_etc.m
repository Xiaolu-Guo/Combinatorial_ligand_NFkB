if 0
%% synergy condition: under 5 by 5 condtion
switch metric_name
    case'peak'
        i_row_plot_cond = 2;
        i_col_plot_cond = 1;
    case 'total'
        i_row_plot_cond = 2;
        i_col_plot_cond = 2;
end
%synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

% score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
synergy_score_nfkb_cond = Synergy_score(:,i_row_plot_cond,i_col_plot_cond);
plot_cell_index_cond = find(synergy_score_nfkb_cond>1.25 & responder_cells_nfkb_cond);

length(plot_cell_index_cond)
sum(responder_cells_nfkb_cond)

% total 28

switch metric_name
    case'peak'
        total_plot_row = 4;
        total_plot_column = 7;

        fold_zoom = 5;
        paperpos=[0,0,140,140];
        papersize=[140 140];
        draw_pos=[10,10,120,120];

    case 'total'
        total_plot_row = 1;
        total_plot_column = 5;

        fold_zoom = 5;
        paperpos=[0,0,100,30];
        papersize=[100 30];
        draw_pos=[10,10,80,10];
end


% fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
fig_handle = plot_single_cell_IKK_TAK_etc(data,plot_cell_index_cond([1,2,3,7,17,28]),i_row_plot_cond+1,i_col_plot_cond+1);
fig_handle 
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

saveas(fig_handle,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_plot_','IKK_TAK_etc_',metric_name),'pdf')
close()
end

if 1
%% antag condition: under 5 by 5 condtion
switch metric_name
    case'peak'
        i_row_plot_cond = 4;
i_col_plot_cond = 5;
    case 'total'
        i_row_plot_cond = 5;
i_col_plot_cond = 5;
end

%synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

% score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
antag_score_nfkb_cond = Antag_score(:,i_row_plot_cond,i_col_plot_cond);
plot_cell_index_cond = find(antag_score_nfkb_cond<0.9 & responder_cells_nfkb_cond);


length(plot_cell_index_cond)
sum(responder_cells_nfkb_cond)

switch metric_name
    case'peak'
        total_plot_row = 6;
        total_plot_column = 6;
    case 'total'
        total_plot_row = 6;
        total_plot_column = 6;
end

fold_zoom = 5;
paperpos=[0,0,140,140];
papersize=[140 140];
draw_pos=[10,10,120,120];

paperpos=[0,0,200,140];
papersize=[200 140];
draw_pos=[10,10,180,120];

[~,ind_sort] = sort(Antag_score(plot_cell_index_cond,i_row_plot_cond,i_col_plot_cond),'descend');
plot_cell_index_cond = plot_cell_index_cond(ind_sort(1:36));

%fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
fig_handle = plot_single_cell_IKK_TAK_etc(data,plot_cell_index_cond([2,12,10,26,35,25]),i_row_plot_cond+1,i_col_plot_cond+1);
fig_handle 
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

saveas(fig_handle ,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_plot_','IKK_TAK_etc_',metric_name),'pdf')
close()
end 

%% subfunction

function [fig_handle] = plot_single_cell_IKK_TAK_etc(data,plot_cell_index_cond,i_row,i_col)

total_plot_row = length(plot_cell_index_cond);

species_vec = [1:6,8,9,12];%1:6
total_plot_column = length(species_vec);
for i_species = 1:total_plot_column
    species_index = species_vec(i_species);

    if ~isempty(plot_cell_index_cond)

        figure(1)
        data_traj = data.model_sim{(i_row - 1)*6 + i_col}(species_index:12:end,:);
        data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(species_index:12:end,:);
        data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(species_index:12:end,:);

        for i_ind = 1:total_plot_row
            i_plot_row = i_ind;
            i_plot_col = i_species;

            subplot(total_plot_row,total_plot_column,(i_plot_row - 1) * total_plot_column + i_plot_col)
            plot(0:5:480,data_traj(plot_cell_index_cond(i_ind),:),'b');hold on
            plot(0:5:480,data_traj_D1(plot_cell_index_cond(i_ind),:), 'r--','LineWidth',2); hold on
            plot(0:5:480,data_traj_D2(plot_cell_index_cond(i_ind),:),'m--','LineWidth',2); hold on
            
            xlim([0,480])
            if i_ind > length(plot_cell_index_cond)
                break
            end


            if i_ind > length(plot_cell_index_cond)
                break
            end
        end
    end
end

fig_handle = gcf;

end