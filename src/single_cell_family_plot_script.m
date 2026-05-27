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
    plot_cell_index_back = find(synergy_score_nfkb_cond<1.0 & responder_cells_nfkb_cond);

    length(plot_cell_index_cond)
    sum(responder_cells_nfkb_cond)

    % total 28
    for i_species = 1:12
        % fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
        fig_handle = plot_single_cell_family_traj(data,plot_cell_index_cond,plot_cell_index_back,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column,i_species,i_species);%3,3
        set(gcf, 'PaperUnits','points')
        fold_zoom = 2.5;
        paperpos=[0,0,50,150];
        papersize=[50,150];
        draw_pos=[10,10,30,130];
        set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

        saveas(fig_handle,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_plot_',data_info.species_outputname{i_species},'_',metric_name),'pdf')
        close()
    end

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
    responder_cells_nfkb_cond_all_cond = Valid_SNR_all_cond(:,i_row+1,i_col+1)...
                                |Valid_SNR_all_cond(:,i_row+1,1)|Valid_SNR_all_cond(:,1,i_col+1);

    % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
    antag_score_nfkb_cond = Antag_score(:,i_row_plot_cond,i_col_plot_cond);
    plot_cell_index_cond = find(antag_score_nfkb_cond<0.9 & responder_cells_nfkb_cond);
    plot_cell_index_back = find(antag_score_nfkb_cond>1.05 & responder_cells_nfkb_cond);

    plot_cell_index_cond = find(antag_score_nfkb_cond<0.9 & responder_cells_nfkb_cond_all_cond);
    plot_cell_index_back = find(antag_score_nfkb_cond>1.05 & responder_cells_nfkb_cond_all_cond);


    length(plot_cell_index_cond)
    sum(responder_cells_nfkb_cond)


    [~,ind_sort] = sort(Antag_score(plot_cell_index_cond,i_row_plot_cond,i_col_plot_cond),'descend');
    plot_cell_index_cond = plot_cell_index_cond(ind_sort(1:36));

    total_plot_row =1;total_plot_column = 1;
    if 0
    for i_species = 1:12
        % fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
        fig_handle = plot_single_cell_family_traj(data,plot_cell_index_cond,plot_cell_index_back,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column,i_species,i_species);%,'switch_order',3,3
        set(gcf, 'PaperUnits','points')
        fold_zoom = 2.5;
        paperpos=[0,0,50,150];
        papersize=[50,150];
        draw_pos=[10,10,30,130];
        set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

        saveas(fig_handle,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_vs_nonatg11_single_cell_plot_',data_info.species_outputname{i_species},'_',metric_name),'pdf')
        close()
    end
    end

    if 1 % plot -5 min to 120 min
    for i_species = 1:12
        save_fig_path_sub_short_period = strcat(save_fig_path_sub,'short_period/');
        % fig_handle = plot_single_cell_traj(data,plot_cell_index_cond,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column);
        fig_handle = plot_single_cell_family_traj_short_period(data,plot_cell_index_cond,plot_cell_index_back,i_row_plot_cond+1,i_col_plot_cond+1,total_plot_row,total_plot_column,i_species,i_species);%,'switch_order',3,3
        set(gcf, 'PaperUnits','points')
        fold_zoom = 2.5;
        paperpos=[0,0,50,150];
        papersize=[50,150];
        draw_pos=[15,10,20,130];
        set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

        saveas(fig_handle,strcat(save_fig_path_sub_short_period,ligand1,ligand2,'_antag_vs_nonatg11_single_cell_plot_',data_info.species_outputname{i_species},'_',metric_name),'pdf')
        close()
    end

    end

end
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

% debug:
% max_vec = max(data_traj_D1,[],2);
% a = find(max_vec == max(max_vec));
% Ikk_a = data_traj_D1(a,:);
% 
% i_species = 1;
% data_traj_D1_nfkB = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
% NFkb_a = data_traj_D1_nfkB(a,:);
% 
% 
% data_traj_nfkB = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
% NFkb_comb = data_traj_nfkB(a,:);
% plot(1:97, NFkb_comb)


%% subfunction
function [fig_handle] = plot_single_cell_family_traj(data,plot_cell_index_cond,plot_cell_index_back,i_row,i_col, total_plot_row,total_plot_column, i_species_start,i_species_end,switch_order)

% if nargin < 5
%
total_plot_row = 3;
total_plot_column = (i_species_end-i_species_start)+1;
% end
Color_vec = {[1,0,0,0.3],[0.7,0.7,0.7,0.3]};
index_plt = {plot_cell_index_cond,plot_cell_index_back};

if nargin == 10
    if strcmp(switch_order,'switch_order')
        Color_vec = {[0.7,0.7,0.7,0.15],[1,0,0,0.15]};
        index_plt = {plot_cell_index_back,plot_cell_index_cond};
    end
end

% if nargin < 7
%     i_species = 1;
% end


ylim_max = -ones(1,(i_species_end-i_species_start)+1);
for i_class_plot = 1:length(index_plt)
    figure(1)
    plot_cell_index = index_plt{i_class_plot};
    for i_col_plot = 1:(i_species_end-i_species_start)+1

        i_species = i_species_start-1 + i_col_plot;
        data_traj = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
        data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
        data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(i_species:12:end,:);
       
        ylim_max(i_class_plot) = max([max(data_traj(plot_cell_index,:),[],'all'),...
            max(data_traj_D1(plot_cell_index,:),[],'all'),...
            max(data_traj_D2(plot_cell_index,:),[],'all')]);
        YL_max = max(ylim_max,[],'all');
        subplot(total_plot_row,total_plot_column,0*total_plot_column + i_col_plot)
        for i_ind = 1:length(plot_cell_index)
            plot(0:5:480,data_traj_D1(plot_cell_index(i_ind),:), 'Color',Color_vec{i_class_plot},'LineWidth',1); hold on
        end
        xlim([0,480])
        ylim([0,YL_max])

        subplot(total_plot_row,total_plot_column,1*total_plot_column + i_col_plot)
        for i_ind = 1:length(plot_cell_index)
            plot(0:5:480,data_traj_D2(plot_cell_index(i_ind),:), 'Color',Color_vec{i_class_plot},'LineWidth',1); hold on
        end
        xlim([0,480])
        ylim([0,YL_max])

        subplot(total_plot_row,total_plot_column,2*total_plot_column + i_col_plot)
        for i_ind = 1:length(plot_cell_index)
            plot(0:5:480,data_traj(plot_cell_index(i_ind),:), 'Color',Color_vec{i_class_plot},'LineWidth',1);hold on
        end
        xlim([0,480])
        ylim([0,YL_max])
    end
end

% if ~isempty(plot_cell_index_cond)
%     i_ind = 1;
%     figure(1)
%     data_traj = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
%     data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
%     data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(i_species:12:end,:);
%
%     for i_plot_row = 1:total_plot_row
%         for i_plot_col = 1:total_plot_column
%
%             subplot(total_plot_row,total_plot_column,(i_plot_row - 1) * total_plot_column + i_plot_col)
%             plot(0:5:480,data_traj(plot_cell_index_cond(i_ind),:),'b');hold on
%             plot(0:5:480,data_traj_D1(plot_cell_index_cond(i_ind),:), 'r--','LineWidth',2); hold on
%             plot(0:5:480,data_traj_D2(plot_cell_index_cond(i_ind),:),'m--','LineWidth',2); hold on
%             i_ind = i_ind + 1;
%             xlim([0,480])
%             if i_ind > length(plot_cell_index_cond)
%                 break
%             end
%
%         end
%         if i_ind > length(plot_cell_index_cond)
%             break
%         end
%     end
% end

fig_handle = gcf;

end

%% subfunction for short period clusters trajectory plot
function [fig_handle] = plot_single_cell_family_traj_short_period(data,plot_cell_index_cond,plot_cell_index_back,i_row,i_col, total_plot_row,total_plot_column, i_species_start,i_species_end,switch_order)

total_plot_row = 3;
total_plot_column = (i_species_end-i_species_start)+1;
% end
Color_vec = {[1,0,0,0.3],[0.7,0.7,0.7,0.15]};
index_plt = {plot_cell_index_cond,plot_cell_index_back};

if nargin == 10
    if strcmp(switch_order,'switch_order')
        Color_vec = {[0.7,0.7,0.7,0.15],[1,0,0,0.15]};
        index_plt = {plot_cell_index_back,plot_cell_index_cond};
    end
end


ylim_max = -ones(1,(i_species_end-i_species_start)+1);
for i_class_plot = 1:length(index_plt)
    figure(1)
    plot_cell_index = index_plt{i_class_plot};
    for i_col_plot = 1:(i_species_end-i_species_start)+1

        i_species = i_species_start-1 + i_col_plot;
        data_traj = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
        data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
        data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(i_species:12:end,:);
       
        ylim_max(i_class_plot) = max([max(data_traj(plot_cell_index,:),[],'all'),...
            max(data_traj_D1(plot_cell_index,:),[],'all'),...
            max(data_traj_D2(plot_cell_index,:),[],'all')]);
        YL_max = max(ylim_max,[],'all');
        subplot(total_plot_row,total_plot_column,0*total_plot_column + i_col_plot)
        for i_ind = 1:length(plot_cell_index)
            plot(-5:5:480,[data_traj_D1(plot_cell_index(i_ind),1),data_traj_D1(plot_cell_index(i_ind),:)], 'Color',Color_vec{i_class_plot},'LineWidth',1); hold on
        end
        xlim([-5,120])
        ylim([0,YL_max])
        set(gca,'XTick',[0:60:120])

        subplot(total_plot_row,total_plot_column,1*total_plot_column + i_col_plot)
        for i_ind = 1:length(plot_cell_index)
            plot(-5:5:480,[data_traj_D2(plot_cell_index(i_ind),1),data_traj_D2(plot_cell_index(i_ind),:)], 'Color',Color_vec{i_class_plot},'LineWidth',1); hold on
        end
        xlim([-5,120])
        ylim([0,YL_max])
        set(gca,'XTick',[0:60:120])

        subplot(total_plot_row,total_plot_column,2*total_plot_column + i_col_plot)
        for i_ind = 1:length(plot_cell_index)
            plot(-5:5:480,[data_traj(plot_cell_index(i_ind),1),data_traj(plot_cell_index(i_ind),:)], 'Color',Color_vec{i_class_plot},'LineWidth',1);hold on
        end
        xlim([-5,120])
        ylim([0,YL_max])
        set(gca,'XTick',[0:60:120])
    end
end

% if ~isempty(plot_cell_index_cond)
%     i_ind = 1;
%     figure(1)
%     data_traj = data.model_sim{(i_row - 1)*6 + i_col}(i_species:12:end,:);
%     data_traj_D1 = data.model_sim{(i_row - 1)*6 + 1}(i_species:12:end,:);
%     data_traj_D2 = data.model_sim{(1 - 1)*6 + i_col}(i_species:12:end,:);
%
%     for i_plot_row = 1:total_plot_row
%         for i_plot_col = 1:total_plot_column
%
%             subplot(total_plot_row,total_plot_column,(i_plot_row - 1) * total_plot_column + i_plot_col)
%             plot(0:5:480,data_traj(plot_cell_index_cond(i_ind),:),'b');hold on
%             plot(0:5:480,data_traj_D1(plot_cell_index_cond(i_ind),:), 'r--','LineWidth',2); hold on
%             plot(0:5:480,data_traj_D2(plot_cell_index_cond(i_ind),:),'m--','LineWidth',2); hold on
%             i_ind = i_ind + 1;
%             xlim([0,480])
%             if i_ind > length(plot_cell_index_cond)
%                 break
%             end
%
%         end
%         if i_ind > length(plot_cell_index_cond)
%             break
%         end
%     end
% end

fig_handle = gcf;

end