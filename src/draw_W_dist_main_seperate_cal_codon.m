% to move to src, if decided to use for publication

% file list
%
% Sim1008_CpG_PolyIC_sc_dose2_IKK_hill2.mat
% Sim1008_LPS_CpG_sc_dose2_IKK_hill2.mat
% Sim1008_LPS_Pam3CSK_sc_dose2_IKK_hill2.mat
% Sim1008_LPS_PolyIC_sc_dose2_IKK_hill2.mat
% Sim1008_Pam3CSK_PolyIC_sc_dose2_IKK_hill2.mat
% Sim1008_TNF_CpG_sc_dose2_IKK_hill2.mat
% Sim1008_TNF_LPS_sc_dose2_IKK_hill2.mat
% Sim1008_TNF_Pam3CSK_sc_dose2_IKK_hill2.mat
% Sim1008_TNF_PolyIC_sc_dose2_IKK_hill2.mat
%
% Sim1008_CpG_PolyIC_sc_signergy_different_dose.mat
% Sim1008_CpG_PolyIC_sc_synergy_different_dose.mat
% Sim1008_LPS_CpG_sc_synergy_different_dose.mat
% Sim1008_LPS_Pam3CSK_sc_synergy_different_dose.mat
% Sim1008_LPS_PolyIC_sc_synergy_different_dose.mat
% Sim1008_Pam3CSK_PolyIC_sc_synergy_different_dose.mat
% Sim1008_TNF_CpG_sc_synergy_different_dose.mat
% Sim1008_TNF_LPS_sc_synergy_different_dose.mat
% Sim1008_TNF_Pam3CSK_sc_synergy_different_dose.mat
% Sim1008_TNF_PolyIC_sc_synergy_different_dose.mat
%
% Sim1008_LPS_Pam3CSK_sc_antag_Diff_CD14.mat
plot_wdist_heatmap = 0;
plot_wdist_hist = 1;
plot_wdist_hist_seperately = 1;
plot_single_ligand_wdist_hist = 0;

if 1
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


        clear metric_cal data_info_cal
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
                i_cond = i_cond + 1;
            end
        end

        % delete here
        %
        % if hill_coeff == 2
        %     %%
        % i_cond = 1;
        % for i_single_ligand = 1:size(single_index_hill2,1)
        %
        %     load(strcat(data_save_file_path,data_file_list{single_index_hill2(i_single_ligand,1)}))
        %     index_to_cal = single_index_hill2(i_single_ligand,2)+1;
        %
        %     if i_cond == 1
        %
        %         field_list = fieldnames(metrics{1});
        %         data_info_cal.info_ligand = data.info_ligand(index_to_cal)';
        %         data_info_cal.info_dose_str = data.info_dose_str(index_to_cal)';
        %         metric_cal = cell(1,2);
        %
        %     else
        %         data_info_cal.info_ligand = [data_info_cal.info_ligand;data.info_ligand(index_to_cal)'];
        %         data_info_cal.info_dose_str = [data_info_cal.info_dose_str;data.info_dose_str(index_to_cal)'];
        %     end
        %
        %
        %     species_num_metric = 1;
        %     for i_id = 1:length(index_to_cal) %1:length(metrics)
        %         i_id_met = index_to_cal(i_id);
        %         data_info_cal.data_label{i_cond} = 'sim_WT';
        %         for i_field = 1:length(field_list)
        %             metric_cal{i_cond}.(field_list{i_field}) = metrics{i_id_met}.(field_list{i_field})(1:species_num_metric:end,:);
        %         end
        %         i_cond = i_cond + 1;
        %     end
        % end
        %
        % for i_data_file_list = 11:length(data_file_list)
        %     load(strcat(data_save_file_path,data_file_list{i_data_file_list}))
        %     % switch i_data_file_list
        %     %     case 1
        %     %         index_to_cal = 2:4;
        %     %     case 3
        %     %
        %     %         index_to_cal = 2:4;
        %     %     case 6
        %     %         index_to_cal = [2,4];
        %     %     case 10
        %     %         index_to_cal = 2:4;
        %     %     case 12
        %     %
        %     %         index_to_cal = 2:4;
        %     %     case 15
        %     %         index_to_cal = [2,4];
        %     %     otherwise
        %     %         index_to_cal = 4;
        %     % end
        %
        %     index_to_cal = 4;
        %
        %
        %     data_info_cal.info_ligand = [data_info_cal.info_ligand;data.info_ligand(index_to_cal)'];
        %     data_info_cal.info_dose_str = [data_info_cal.info_dose_str;data.info_dose_str(index_to_cal)'];
        %
        %
        %     species_num_metric = 1;
        %     for i_id = 1:length(index_to_cal) %1:length(metrics)
        %         i_id_met = index_to_cal(i_id);
        %         data_info_cal.data_label{i_cond} = 'sim_WT';
        %         for i_field = 1:length(field_list)
        %             metric_cal{i_cond}.(field_list{i_field}) = metrics{i_id_met}.(field_list{i_field})(1:species_num_metric:end,:);
        %         end
        %         i_cond = i_cond + 1;
        %     end
        % end
        %
        % end

        collect_feature_vects = calculate_codon_from_metric20251118(data_info_cal,metric_cal);

        if 0
            %% plot codon distribution
            codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

            for i_codon = 1:length(codon_list)
                figure(1)
                paperpos=[0,0,800,100];
                papersize=[800 100];
                draw_pos=[10,10,780,80];
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

                clear y

                y = collect_feature_vects.(codon_list{i_codon});

                % z = collect_feature_vects.(codon_list{i_codon}); % 1000 by 10
                std_cal_y = [];
                for i_y = 1:length(y)
                    std_cal_y = [std_cal_y;y{i_y}];
                end
                std_y = std(std_cal_y(:))/4000;

                % subplot(1,length(vis_data_field),i_data_field)

                al_goodplot_pair_RMSD_diff_size(y,[1:15],0.5,ones(size(y,2),1)*[255 0 0]/255 ,'bilateral',[],std_y); %left
                % al_goodplot_pair_RMSD_diff_size(z,[2:3:5],0.5,ones(size(z,2),1)*[0 0 0]/255,'bilateral',[],std_y);


                xlim([0 16]);

                xticks([0:16]);
                xticklabels({});
                %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

                ylim([-2.1,2.1]);
                % for i_x = 1:10
                %     plot([i_x,i_x],[0,5],'--','Color','k');hold on
                % end
                set(gca,'fontsize',14,'fontname','Arial');
                %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

                saveas(gcf,strcat(fig_save_path_sub,'codon_distrib_hill',num2str(hill_coeff),'/Codon_',codon_list{i_codon}),'pdf');%,save_name_tag
                close

            end
        end

        if 1

            %% plot W-dist

            plot_opt.cal_Wdis_mat = 1;
            plot_opt.save_Wdis_mat = 1;
            plot_opt.data_save_file_path = data_save_file_path;
            plot_opt.plot_wdis_mat = 1;
            plot_opt.plot_hierach_tree = 0;
            plot_opt.fig_save_path = fig_save_path;
            plot_opt.fig_save_name_wdis_mat = strcat('IKK_hill',num2str(hill_coeff),'_10cond_comb_1119');
            wdis_mat_name = strcat('IKK_hill',num2str(hill_coeff),'_10cond_Wdist_comb.mat');
            w3 = draw_Wdist_mat(collect_feature_vects,wdis_mat_name,plot_opt);


        end

        collect_feature_vects_cell{hill_coeff} = collect_feature_vects;
        w_dis_cell{hill_coeff} = w3;
    end



    if 1

        %% bar plot of W-dist
        %[1,2,5,6,10]
        % w3_1 = w3([3,4,7:9,11:14],[3,4,7:9,11:14]);
        % w3_2 = w3((14+[3,4,7:9,11:14]),(14+[3,4,7:9,11:14]));
        % histogram(w3_1);hold on
        % histogram(w3_2)


        if plot_wdist_hist
            w3_1 = w_dis_cell{1}(1:15,1:15);
            w3_2 = w_dis_cell{2}(1:15,1:15);

            w3_1_vec = w3_1(tril(true(size(w3_1))));
            w3_2_vec = w3_2(tril(true(size(w3_2))));

            figure (1)
            set(gcf, 'PaperUnits','points')

            zoom_folder = 0.7;
            paper_pos = [0,0,200,150];
            paper_size = [200,150];
            set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

            histogram(w3_1_vec,0:0.2:1.0);hold on
            histogram(w3_2_vec,0:0.2:1.0)

            if 0 % run Wilcoxon signed-rank tests
                [p, ~, ~] = signrank(w3_1_vec, w3_2_vec)
            end
            %legend('Hill = 1','Hill = 2')
            YL = ylim();
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_Wdis_hill12_with_label'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_Wdis_hill12_no_label'),'pdf')
            close()
        end

        if plot_wdist_hist_seperately

                        w3_1 = w_dis_cell{1}(1:15,1:15);
            w3_2 = w_dis_cell{2}(1:15,1:15);

            w3_1_vec = w3_1(tril(true(size(w3_1))));
            w3_2_vec = w3_2(tril(true(size(w3_2))));

            figure (1)
            set(gcf, 'PaperUnits','points')

            zoom_folder = 0.7;
            paper_pos = [0,0,200,150];
            paper_size = [200,150];
            set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

            histogram(w3_1_vec,0:0.2:1.0);hold on
            % histogram(w3_2_vec,0:0.2:1.0)
            %legend('Hill = 1','Hill = 2')
            ylim(YL)
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_Wdis_hill1_with_label'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_Wdis_hill1_no_label'),'pdf')
            close()

            figure (1)
            set(gcf, 'PaperUnits','points')

            zoom_folder = 0.7;
            paper_pos = [0,0,200,150];
            paper_size = [200,150];
            set(gcf, 'PaperPosition', zoom_folder* paper_pos,'PaperSize',zoom_folder*paper_size,'Position',zoom_folder*[10,10,180,130])

            %histogram(w3_1_vec,0:0.2:1.0);hold on
            histogram(w3_2_vec,0:0.2:1.0)
            ylim(YL)
            %legend('Hill = 1','Hill = 2')
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_Wdis_hill2_with_label'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_Wdis_hill2_no_label'),'pdf')
            close()

        end


        if plot_single_ligand_wdist_hist
            w3_1_single = w_dis_cell{1}(1:5,1:5);
            w3_2_single = w_dis_cell{2}(1:5,1:5);

            w3_1_single_vec = w3_1_single(tril(true(size(w3_1_single))));
            w3_2_single_vec = w3_2_single(tril(true(size(w3_2_single))));

            figure (1)
            set(gcf, 'PaperUnits','points')

            paper_pos = [0,0,200,150];
            paper_size = [200,150];
            set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size,'Position',[10,10,180,130])

            histogram(w3_1_single_vec,0:0.2:1.0);hold on
            histogram(w3_2_single_vec,0:0.2:1.0)
            % legend('Hill = 1','Hill = 2')
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_single_Wdis_hill12_with_label'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_single_Wdis_hill12_no_label'),'pdf')
            close()
        end


        if 0 % to delete? 12/01/2025

            w3_1_dual = w_dis_cell{1}(6:15,6:15);
            w3_2_dual = w_dis_cell{2}(6:15,6:15);

            w3_1_dual_vec = w3_1_dual(tril(true(size(w3_1_dual))));
            w3_2_dual_vec = w3_2_dual(tril(true(size(w3_2_dual))));

            figure (1)
            set(gcf, 'PaperUnits','points')

            paper_pos = [0,0,200,150];
            paper_size = [200,150];
            set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size,'Position',[10,10,180,130])

            histogram(w3_1_dual_vec,0:0.2:1.0);hold on
            histogram(w3_2_dual_vec,0:0.2:1.0)
            legend('Hill = 1','Hill = 2')
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_dual_Wdis_hill12_with_label'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_dual_Wdis_hill12_no_label'),'pdf')
            close()

        end


        if 0 % to delete? 12/01/2025

            dual_to_single_index = [6,1,2;
                7,1,3;
                8,1,4;
                9,1,5;
                10,2,3;
                11,2,4;
                12,2,5;
                13,3,4;
                14,3,5;
                15,4,5];
            figure (1)
            set(gcf, 'PaperUnits','points')

            paper_pos = [0,0,1000,150];
            paper_size = [1000,150];
            set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size,'Position',[10,10,980,130])

            for i_pair = 1:length(dual_to_single_index)
                w3_1_dual = w_dis_cell{1}(dual_to_single_index(i_pair,:),dual_to_single_index(i_pair,:));
                w3_2_dual = w_dis_cell{2}(dual_to_single_index(i_pair,:),dual_to_single_index(i_pair,:));

                w3_1_dual_vec = w3_1_dual(tril(true(size(w3_1_dual))));
                w3_2_dual_vec = w3_2_dual(tril(true(size(w3_2_dual))));

                subplot(1,10,i_pair)

                histogram(w3_1_dual_vec,0:0.2:1.0);hold on
                histogram(w3_2_dual_vec,0:0.2:1.0)
                legend('Hill = 1','Hill = 2')

            end
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_dual_Wdis_hill12_with_label'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(plot_opt.fig_save_path,'hist_dual_Wdis_hill12_no_label'),'pdf')
            close()

        end

    end

end
