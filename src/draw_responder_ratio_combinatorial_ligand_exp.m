% version 202502

%% initializing do not copy this
%debug_or_not = 0;

addpath('../raw_data2025/');%_fay_parameter/';
%fig_save_path = fig_save_path_diff;

% addpath('./lib/')
% addpath('./src/')
% addpath('./bin/')

Colors_vecs = get(groot, 'defaultAxesColorOrder');
% load('Supriya_data_metrics.mat')
% supriya_index = [4;27;22;6;15];
% supriya_dual_ligand_index = [3;5;8;9;26;14;25;18;19];

% % sampled data
% data_filename = 'Sim5_codon_all5dose_metric.mat';
% load(strcat(data_save_file_path_1,data_filename))
% sample_dual_ligand_index = [1;2;3;4;10;9;6;8;5;7];
% sample_index = [13;12;11;14;15];
        dot_size = 12;

%% draw exp. multi ligand
if 1 % compare with exp. 3 ligand
    load('Exp_data_metrics_2024.mat')
    %     data_index = [19;8;18];
    %     data_index = [4;22;6;15]; % 11
    cal_responder_ratio = 1;
    index_exp_ligand_vects = {[4];% 1-TNF
        [6,24]; %2- lps
        [1,10,22]; %3 cpg
        [7,12,15,17]; %4
        [2,11,23,27]; %5
        [8];%6
        [5];%7
        [9];%8
        [3];%9
        [25];%10
        [19];%11
        [26];%12
        [13,16,18];%13
        [];%14
        [14];%15
        [];%16
        [];%17
        [];%18
        [20];%19
        [];%20
        [21];%21
        [];%22
        [];%23
        [];%24
        [];%25
        [];%26
        [];%27
        [];%28
        [];%29
        [];%30
        [28]};%31
    if cal_responder_ratio

        %% plot lines

        index_exp_selected = {[4];% 1-TNF
            [6]; %2- lps
            [22]; %3 cpg
            [15]; %4
            [11]; %5
            [8];%6
            [5];%7
            [9];%8
            [3];%9
            [25];%10
            [19];%11
            [26];%12
            [18];%13
            [];%14
            [14];%15
            [];%16
            [];%17
            [];%18
            [20];%19
            [];%20
            [21];%21
            [];%22
            [];%23
            [];%24
            [];%25
            [];%26
            [];%27
            [];%28
            [];%29
            [];%30
            [28]};%31
        plot_lines = 1;
        if plot_lines
            figure(1)
            paperpos = [0,0,220,50]*1.8;
            papersize = [220,50]*1.8;
            draw_pos = [10,10,200,30]*1.8;

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

            figure(2)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)



            integral_pos_threshold_vec = [0.2,0.3,0.40];
            peak_threshold_vec = [0.14,0.18,0.22];

            % Use arrayfun for element-wise operation
            integral_pos_threshold_legend = arrayfun(@num2str, integral_pos_threshold_vec, 'UniformOutput', false);
            peak_threshold_legend = arrayfun(@num2str, peak_threshold_vec, 'UniformOutput', false);

            for i_thresh = 1:length(integral_pos_threshold_vec)
                integral_pos_threshold = integral_pos_threshold_vec(i_thresh);
                peak_threshold = peak_threshold_vec(i_thresh);

                clear responder_ratio responder_ratio_scatter_y responder_ratio_scatter_x
                clear responder_ratio_by_peak_scatter_y responder_ratio_by_peak_scatter_x responder_ratio_by_peak
                index_scatter_y = 1;
                for i_cond = 1:length(index_exp_selected)
                    index_exp_ligand = index_exp_selected{i_cond};
                    responder_ratio(i_cond) = 0;
                    responder_ratio_by_peak(i_cond) = 0;
                    for i_rpc = 1:length(index_exp_ligand)
                        responder_ratio_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end) > integral_pos_threshold)/length(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end));
                        responder_ratio_scatter_x(index_scatter_y) = i_cond;
                        responder_ratio(i_cond) =responder_ratio(i_cond)+ responder_ratio_scatter_y(index_scatter_y);

                        responder_ratio_by_peak_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end) > peak_threshold)/length(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end));
                        responder_ratio_by_peak_scatter_x(index_scatter_y) = i_cond;
                        responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)+responder_ratio_by_peak_scatter_y(index_scatter_y);
                        index_scatter_y = index_scatter_y+1;
                    end

                    if responder_ratio(i_cond)
                        responder_ratio(i_cond) = responder_ratio(i_cond)/length(index_exp_ligand);
                        responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)/length(index_exp_ligand);
                    end
                end

                responder_ratio_exp{i_thresh} = responder_ratio;
                responder_ratio_scatter_y_exp{i_thresh} = responder_ratio_scatter_y;
                responder_ratio_scatter_x_exp{i_thresh} = responder_ratio_scatter_x;
                responder_ratio_by_peak_exp{i_thresh} = responder_ratio_by_peak;


                one_index = 1:5;
                two_index = [6:13,15];
                three_index = [19,21];
                four_index = [];
                five_more_index = 31;
                respondera_ratio_bar_exp(i_thresh,:) = [mean(responder_ratio(one_index)),...
                    mean(responder_ratio(two_index)),...
                    mean(responder_ratio(three_index)),...
                    mean(responder_ratio(four_index)),...
                    mean(responder_ratio(five_more_index))];
                respondera_ratio_bar_by_peak_exp(i_thresh,:) = [mean(responder_ratio_by_peak(one_index)),...
                    mean(responder_ratio_by_peak(two_index)),...
                    mean(responder_ratio_by_peak(three_index)),...
                    mean(responder_ratio_by_peak(four_index)),...
                    mean(responder_ratio_by_peak(five_more_index))];


                figure(1)
                responder_ratio_nonzero_index = find(responder_ratio > 0);
                %plot(responder_ratio_nonzero_index,responder_ratio(responder_ratio_nonzero_index),'LineWidth',2);hold on
                % scatter(responder_ratio_scatter_x,responder_ratio_scatter_y,'filled');hold
                % figure(1)
                % plot(1:length(responder_ratio),responder_ratio,'LineWidth',2);hold on
                scatter(responder_ratio_nonzero_index,responder_ratio(responder_ratio_nonzero_index),dot_size,Colors_vecs(i_thresh,:),'filled') ;hold on

                figure(2)
                responder_ratio_by_peak_nonzero_index = find(responder_ratio_by_peak > 0);
                %plot(1:length(responder_ratio_by_peak),responder_ratio_by_peak,'LineWidth',2);hold on
                scatter(responder_ratio_by_peak_nonzero_index,responder_ratio_by_peak(responder_ratio_by_peak_nonzero_index),dot_size,Colors_vecs(i_thresh,:),'filled') ;hold on

            end

            length_responder_ratio_total = 31;
            Xtick_labels ={};
            for i_Xtick_labels = 1:length_responder_ratio_total
                Xtick_labels(i_Xtick_labels) = {''};
            end

            figure(1)
            %legend(integral_pos_threshold_legend)

            set(gca,'XTick',1:length_responder_ratio_total,'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
                'FontSize',8)
            ylim([0,1])
            xlim([0.5,31.5])
            %saveas(gcf,strcat(fig_save_path,'Responder_ratio_scatter_by_integral_exp_with7sti'),'pdf')
            close()

            figure(2)
            %legend(integral_pos_threshold_legend)

            set(gca,'XTick',1:length_responder_ratio_total,'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
                'FontSize',8)
            ylim([0,1])
            xlim([0.5,31.5])
            %saveas(gcf,strcat(fig_save_path,'Responder_ratio_scatter_by_peak_exp_with7sti'),'pdf')
            close()

            paperpos = [0,0,70,50]*1.8;
            papersize = [70,50]*1.8;
            draw_pos = [10,10,50,30]*1.8;

            figure(1)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

            for i_bar = 1:size(respondera_ratio_bar_exp,1)
                plot([1:3,5], respondera_ratio_bar_exp(i_bar,[1:3,5]),'LineWidth',2);hold on
            end

            set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
                'FontSize',8)

            ylim([0,1])
            xlim([1,5])
            %saveas(gcf,strcat(fig_save_path,'Responder_ratio_ligand_number_curve_by_integral_exp_match'),'pdf')
            close()

            figure(2)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

            for i_bar = 1:size(respondera_ratio_bar_by_peak_exp,1)
                plot([1:3,5], respondera_ratio_bar_by_peak_exp(i_bar,[1:3,5]),'LineWidth',2);hold on
            end

            set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
                'FontSize',8)

            ylim([0,1])
            xlim([1,5])
            %saveas(gcf,strcat(fig_save_path,'Responder_ratio_ligand_number_curve_by_peak_exp_match'),'pdf')
            close()

            responder_ratio_exp_by_integral = responder_ratio_exp{2};
            responder_ratio_exp_by_peak = responder_ratio_by_peak_exp{2};
            %save(strcat('../raw_data2025/','exp_responder_ratio_med'),"responder_ratio_exp_by_integral","responder_ratio_exp_by_peak");
            %save(strcat('../raw_data2025/','sim_responder_ratio_med'),"responder_ratio_sim_by_integral","responder_ratio_sim_by_peak");

        end

    end

end
