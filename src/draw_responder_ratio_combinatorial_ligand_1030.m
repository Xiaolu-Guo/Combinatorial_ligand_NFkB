% version 20251030; for updated model simulation with pIC-TLR3 deg increase
% requires variables: save_metric_name, fig_save_path ,data_save_file_path 

Colors_vecs = get(groot, 'defaultAxesColorOrder');

%% draw multi ligand stim match
if 1 % draw all comb ligand stim for weighted

    cal_responder_ratio = 1;
    integral_pos_threshold_vec = [0.4,0.5,0.6];
    peak_threshold_vec = [0.12,0.14,0.16];
        dot_size = 12;
    % Use arrayfun for element-wise operation
    integral_pos_threshold_legend = arrayfun(@num2str, integral_pos_threshold_vec, 'UniformOutput', false);
    peak_threshold_legend = arrayfun(@num2str, peak_threshold_vec, 'UniformOutput', false);

    if cal_responder_ratio

        paperpos = [0,0,220,50]*1.8;
        papersize = [220,50]*1.8;
        draw_pos = [10,10,200,30]*1.8;

        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        figure(2)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        load(strcat(data_save_file_path,save_metric_name))

        clear respondera_ratio_bar_sim respondera_ratio_bar_by_peak_sim
        for i_thresh = 1:length(integral_pos_threshold_vec)
            integral_pos_threshold = integral_pos_threshold_vec(i_thresh);
            peak_threshold = peak_threshold_vec(i_thresh);


            for i_cond = 1:length(metrics)
                responder_ratio(i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));

            end

            responder_ratio_sim{i_thresh} = responder_ratio;
            responder_ratio_by_peak_sim{i_thresh} = responder_ratio_by_peak;

            respondera_ratio_bar_sim(i_thresh,:) = [mean(responder_ratio(1:5)),...
                mean(responder_ratio(6:15)),...
                mean(responder_ratio(16:25)),...
                mean(responder_ratio(26:30)),...
                mean(responder_ratio(31))];
            respondera_ratio_bar_by_peak_sim(i_thresh,:) = [mean(responder_ratio_by_peak(1:5)),...
                mean(responder_ratio_by_peak(6:15)),...
                mean(responder_ratio_by_peak(16:25)),...
                mean(responder_ratio_by_peak(26:30)),...
                mean(responder_ratio_by_peak(31))];

            figure(1)
            % plot(1:length(responder_ratio),responder_ratio,'LineWidth',2);hold on
            scatter(1:length(responder_ratio),responder_ratio,dot_size,Colors_vecs(i_thresh,:),'filled') ;hold on

            figure(2)
            %plot(1:length(responder_ratio_by_peak),responder_ratio_by_peak,'LineWidth',2);hold on
            scatter(1:length(responder_ratio_by_peak),responder_ratio_by_peak,dot_size,Colors_vecs(i_thresh,:),'filled') ;hold on

        end

        Xtick_labels ={};
        for i_Xtick_labels = 1:length(responder_ratio)
            Xtick_labels(i_Xtick_labels) = {''};
        end

        figure(1)

        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)
        ylim([0,1])
        xlim([0.5,31.5])
        saveas(gcf,strcat(fig_save_path,'Responder_ratio_scatter_by_integral_sim_match'),'pdf')
        close()

        figure(2)
        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)

        ylim([0,1])
        xlim([0.5,31.5])

        saveas(gcf,strcat(fig_save_path,'Responder_ratio_scatter_by_peak_sim_match'),'pdf')
        close()


        paperpos = [0,0,70,50]*1.8;
        papersize = [70,50]*1.8;
        draw_pos = [10,10,50,30]*1.8;

        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        for i_bar = 1:size(respondera_ratio_bar_sim,1)
            plot(1:size(respondera_ratio_bar_sim,2), respondera_ratio_bar_sim(i_bar,:),'LineWidth',2);hold on
        end

        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)

        ylim([0,1])
        xlim([1,5])
        saveas(gcf,strcat(fig_save_path,'Responder_ratio_ligand_number_curve_by_integral_sim_match'),'pdf')
        close()

        figure(2)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        for i_bar = 1:size(respondera_ratio_bar_by_peak_sim,1)
            plot(1:size(respondera_ratio_bar_by_peak_sim,2), respondera_ratio_bar_by_peak_sim(i_bar,:),'LineWidth',2);hold on
        end

        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)

        ylim([0,1])
        %xlim([0,6])
        xlim([1,5])
        saveas(gcf,strcat(fig_save_path,'Responder_ratio_ligand_number_curve_by_peak_sim_match'),'pdf')
        close()


        responder_ratio_sim_by_integral = responder_ratio_sim{2};
        responder_ratio_sim_by_peak = responder_ratio_by_peak_sim{2};
        save(strcat(data_save_file_path,'sim_responder_ratio_med'),"responder_ratio_sim_by_integral","responder_ratio_sim_by_peak");

    end

end



