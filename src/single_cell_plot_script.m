% update name, move to /src
% relpace/ delete / or move to debug: the original codes single_cell_plot_script.m 
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
        paperpos=[0,0,140,80];
        papersize=[140 80];
        draw_pos=[10,10,120,60];

    case 'total'
        total_plot_row = 1;
        total_plot_column = 5;

        fold_zoom = 5;
        paperpos=[0,0,100,30];
        papersize=[100 30];
        draw_pos=[10,10,80,10];
end


fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
fig_handle
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

saveas(fig_handle,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_single_cell_plot_','_',metric_name),'pdf')
close()

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

[~,ind_sort] = sort(Antag_score(plot_cell_index_cond,i_row_plot_cond,i_col_plot_cond),'descend');
plot_cell_index_cond = plot_cell_index_cond(ind_sort(1:36));

fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

saveas(fig_handle ,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_single_cell_plot_','_',metric_name),'pdf')
close()

%% subfunction

function [fig_handle] = plot_single_cell_traj(data,plot_cell_index_cond,i_row,i_col, total_plot_row,total_plot_column, i_species)

% if nargin < 5
% 
% total_plot_row = 4;
% total_plot_column = 7;
% end


if nargin < 7
    i_species = 1;
end


if ~isempty(plot_cell_index_cond)
    i_ind = 1;
    figure(1)
    data_traj = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
    data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
    data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(i_species:12:end,:);

    for i_plot_row = 1:total_plot_row
        for i_plot_col = 1:total_plot_column

            subplot(total_plot_row,total_plot_column,(i_plot_row - 1) * total_plot_column + i_plot_col)
            plot(0:5:480,data_traj(plot_cell_index_cond(i_ind),:),'b');hold on
            plot(0:5:480,data_traj_D1(plot_cell_index_cond(i_ind),:), 'r--','LineWidth',2); hold on
            plot(0:5:480,data_traj_D2(plot_cell_index_cond(i_ind),:),'m--','LineWidth',2); hold on
            i_ind = i_ind + 1;
            xlim([0,480])
            if i_ind > length(plot_cell_index_cond)
                break
            end

        end
        if i_ind > length(plot_cell_index_cond)
            break
        end
    end
end

fig_handle = gcf;

end