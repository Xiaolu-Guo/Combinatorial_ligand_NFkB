% to move to src, if decided to use for publication
fig_save_path_sub = fig_save_path;

% data_save_file_path,data_file_list,fig_save_path
data_file_list = {'Sim1008_TNF_LPS_sc_dose2_IKK_hill1.mat'
    'Sim1008_TNF_CpG_sc_dose2_IKK_hill1.mat'
    'Sim1008_TNF_PolyIC_sc_dose2_IKK_hill1.mat'
    'Sim1008_TNF_Pam3CSK_sc_dose2_IKK_hill1.mat'
    'Sim1008_LPS_CpG_sc_dose2_IKK_hill1.mat'
    'Sim1008_LPS_PolyIC_sc_dose2_IKK_hill1.mat'
    'Sim1008_LPS_Pam3CSK_sc_dose2_IKK_hill1.mat'
    'Sim1008_CpG_PolyIC_sc_dose2_IKK_hill1.mat'
    'Sim1008_CpG_Pam3CSK_sc_dose2_IKK_hill1.mat'
    'Sim1008_Pam3CSK_PolyIC_sc_dose2_IKK_hill1.mat'
    'Sim1008_TNF_LPS_sc_dose2_IKK_hill2.mat'
    'Sim1008_TNF_CpG_sc_dose2_IKK_hill2.mat'
    'Sim1008_TNF_PolyIC_sc_dose2_IKK_hill2.mat'
    'Sim1008_TNF_Pam3CSK_sc_dose2_IKK_hill2.mat'
    'Sim1008_LPS_CpG_sc_dose2_IKK_hill2.mat'
    'Sim1008_LPS_PolyIC_sc_dose2_IKK_hill2.mat'
    'Sim1008_LPS_Pam3CSK_sc_dose2_IKK_hill2.mat'
    'Sim1008_CpG_PolyIC_sc_dose2_IKK_hill2.mat'
    'Sim1008_CpG_Pam3CSK_sc_dose2_IKK_hill2.mat'
    'Sim1008_Pam3CSK_PolyIC_sc_dose2_IKK_hill2.mat'};

species_num = 12;

single_index_hill(:,:,1) = [1,1;1,2;2,2;3,2;4,2];
single_index_hill(:,:,2) =  [11,1;11,2;12,2;13,2;14,2];

for hill_coeff = 1:2


    clear metric_cal data_info_cal total_activity_single total_activity_comb
    total_activity_single = [];
    total_activity_comb = [];
    %%
    i_cond = 1;
    for i_single_ligand = 1:size(single_index_hill(:,:,hill_coeff),1)

        load(strcat(data_save_file_path,data_file_list{single_index_hill(i_single_ligand,1,hill_coeff)}))
        index_to_cal = single_index_hill(i_single_ligand,2,hill_coeff)+1;

        if i_cond == 1

            field_list = fieldnames(metrics{1});
            data_info_cal.info_ligand = data.info_ligand(index_to_cal)';
            data_info_cal.info_dose_str = data.info_dose_str(index_to_cal)';
            metric_cal = cell(1,2);

        else
            data_info_cal.info_ligand = [data_info_cal.info_ligand;data.info_ligand(index_to_cal)'];
            data_info_cal.info_dose_str = [data_info_cal.info_dose_str;data.info_dose_str(index_to_cal)'];
        end

        species_num_metric = 1;
        for i_id = 1:length(index_to_cal) %1:length(metrics)
            i_id_met = index_to_cal(i_id);
            data_info_cal.data_label{i_cond} = 'sim_WT';
            for i_field = 1:length(field_list)
                metric_cal{i_cond}.(field_list{i_field}) = metrics{i_id_met}.(field_list{i_field})(1:species_num_metric:end,:);
            end
            total_activity_single = [total_activity_single;metrics{i_id_met}.max_pos_integral(1:species_num_metric:end,:)];
            i_cond = i_cond + 1;
        end
    end


    for i_data_file_list = ((hill_coeff-1) * 10 + 1):((hill_coeff-1) * 10 +length(data_file_list)/2)
        load(strcat(data_save_file_path,data_file_list{i_data_file_list}))
        % switch i_data_file_list
        %     case 1
        %         index_to_cal = 2:4;
        %     case 3
        %
        %         index_to_cal = 2:4;
        %     case 6
        %         index_to_cal = [2,4];
        %     case 10
        %         index_to_cal = 2:4;
        %     case 12
        %
        %         index_to_cal = 2:4;
        %     case 15
        %         index_to_cal = [2,4];
        %     otherwise
        %         index_to_cal = 4;
        % end

        index_to_cal = 4;


        data_info_cal.info_ligand = [data_info_cal.info_ligand;data.info_ligand(index_to_cal)'];
        data_info_cal.info_dose_str = [data_info_cal.info_dose_str;data.info_dose_str(index_to_cal)'];


        species_num_metric = 1;
        for i_id = 1:length(index_to_cal) %1:length(metrics)
            i_id_met = index_to_cal(i_id);
            data_info_cal.data_label{i_cond} = 'sim_WT';
            for i_field = 1:length(field_list)
                metric_cal{i_cond}.(field_list{i_field}) = metrics{i_id_met}.(field_list{i_field})(1:species_num_metric:end,:);
            end
            total_activity_comb = [total_activity_comb;metrics{i_id_met}.max_pos_integral(1:species_num_metric:end,:)];
            i_cond = i_cond + 1;
        end
    end

    collect_feature_vects_integral_single{hill_coeff} = total_activity_single;
    collect_feature_vects_integral_comb{hill_coeff} = total_activity_comb;
collect_feature_vects_integral_all{hill_coeff} = [total_activity_single;total_activity_comb];

end




figure (1)
set(gcf, 'PaperUnits','points')

zoom_folder = 0.7;
paper_pos = [0,0,200,150];
paper_size = [200,150];
set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

hill_coeff = 1;
histogram(collect_feature_vects_integral_all{hill_coeff},0:0.1:4 );hold on
hill_coeff = 2;
histogram(collect_feature_vects_integral_all{hill_coeff} ,0:0.1:4)
%legend('Hill = 1','Hill = 2')
saveas(gcf,strcat(plot_opt.fig_save_path,'hist_integral_hill12_with_label'),'pdf')
set(gca,'XTickLabel',{},'YTickLabel',{})
saveas(gcf,strcat(plot_opt.fig_save_path,'hist_integral_hill12_no_label'),'pdf')
close()

