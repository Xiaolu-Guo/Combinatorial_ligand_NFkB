if 0 % updated & tested : 12/01
    %% synergy condition: under 5 by 5 condtion

    fold_zoom = 13;
    paperpos=[0,0,110,100];
    papersize=[110 100];
    draw_pos=[10,10,100,90];



    switch metric_name
        case'peak'
            i_row_plot_cond = 2;
            i_col_plot_cond = 1;

        case 'total'
            i_row_plot_cond = 2;
            i_col_plot_cond = 2;
    end
    i_cond = i_row_plot_cond * 6 + i_col_plot_cond + 1;

    params_cond = gene_info_all{i_cond}.gene_parameter_value_vec_genotype{1}{1}';
    params_names = gene_info_all{i_cond}.parameter_name_vec{1};

    % remove extra useless parameters: first p52n2, overwrite by second p52n2
    params_cond = params_cond(:,[1:2,4:end]);
    params_names = params_names([1:2,4:end]);

    params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);
    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
    responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

    % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
    synergy_score_nfkb_cond = Synergy_score(:,i_row_plot_cond,i_col_plot_cond);
    plot_cell_index_cond = synergy_score_nfkb_cond>1.25 & responder_cells_nfkb_cond;
    plot_cell_index_back = synergy_score_nfkb_cond<1.0 & responder_cells_nfkb_cond;

    if 0 % plot scatter plots and pdf

        %%
        len_para = length(params_names);
        W_diff = zeros(1,len_para);

        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

        for i_para = 1:len_para
            for j_para = i_para:len_para

                subplot(len_para,len_para,(i_para-1)* len_para+ j_para)
                if i_para == j_para
                    nbins = 10;
                    para_all_responder = params_cond(plot_cell_index_back,i_para);
                    para_synergy_antag = params_cond(plot_cell_index_cond,i_para);
                    %edges = linspace(min(para_all_responder),
                    %max(min(para_all_responder)), nbins+1); ,edges
                    % Estimate and plot the distribution for para_all_responder

                    [f1, x1] = ksdensity(log10(para_all_responder));
                    plot(x1, f1, 'Color',[0.7,0.7,0.7],'LineWidth', 2); hold on;

                    % Estimate and plot the distribution for para_synergy_antag
                    f2 = ksdensity(log10(para_synergy_antag),x1);
                    plot(x1, f2,'Color',[1,0,0], 'LineWidth', 2); hold on;
                    xlabel(strcat('log10(',params_names_label{i_para},')'))
                    ylabel('pdf')
                    % % Add labels and legend
                    % xlabel('Parameter Value');
                    % ylabel('Probability Density');
                    % legend('All Responders', 'Synergy/Antagonism');
                    % title('Smoothed Distribution Curves');

                    % histogram(para_all_responder,'Normalization', 'pdf');hold on
                    % histogram(para_synergy_antag,'Normalization', 'pdf');hold on
                    % set(gca,'XScale','log')
                    %% synergistic cells distribution vs default distribution

                    W_diff(i_para) = w_distance(log10(para_all_responder),log10(para_synergy_antag), 2);

                else
                    figure(1)
                    scatter(log10(params_cond(plot_cell_index_back,j_para)),log10(params_cond(plot_cell_index_back,i_para)),16,[0.7,0.7,0.7],'filled');hold on
                    scatter(log10(params_cond(plot_cell_index_cond,j_para)),log10(params_cond(plot_cell_index_cond,i_para)),16,[1,0,0],'filled')
                    %set(gca,'XScale','log','YScale','log')
                    %set(gca,'XTickLabel',{},'YTickLabel',{})
                    xlabel(strcat('log10(',params_names_label{j_para},')'))
                    ylabel(strcat('log10(',params_names_label{i_para},')'))

                end

            end
        end

        figure(1)
        saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_distrib_',metric_name),'pdf')
        close()

        figure(2)
        set(gcf, 'PaperUnits','points')
        fold_zoom_fig2 = 1.5;
        set(gcf, 'PaperPosition', paperpos*fold_zoom_fig2,'PaperSize', papersize*fold_zoom_fig2,'Position',draw_pos*fold_zoom_fig2)

        bar(W_diff)
        set(gca,'XTickLabel',params_names_label)
        ylabel('W-Distance')
        %ylim([0,0.7])
        saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_WDist_',metric_name),'pdf')
        close()

    end

    if 1 % plot single parameter distribution
%%
        plot_para_index = [1,6,7,8];

        for i_para = 1:length(plot_para_index)

            index_para = plot_para_index(i_para);

            para_all_responder = params_cond(plot_cell_index_back,index_para);
            para_synergy_antag = params_cond(plot_cell_index_cond,index_para);


            [f1, x1] = ksdensity(log10(para_all_responder));

            f2 = ksdensity(log10(para_synergy_antag),x1);

            figure

            fold_zoom_1 = 0.8;
            paperpos=[0,0,110,100];
            papersize=[110 100];
            draw_pos=[10,10,90,80];
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)

            plot(x1, f1, 'Color',[0.7,0.7,0.7],'LineWidth', 1); hold on;

            plot(x1, f2,'Color',[1,0,0], 'LineWidth', 1); hold on;
            %xlabel(strcat('-log10(',params_names_label{1},'*',params_names_label{6},')'))
            %ylabel('pdf')
            saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_',params_names_label{index_para},metric_name,'_withlabel'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_',params_names_label{index_para},metric_name),'pdf')
            

            close()
        end
    end

    if 0 % to delete not useful % plot hist of k99 * k101,
%%
        W_diff = 0;
        if 1 % k99 * k101

            para_all_responder = params_cond(plot_cell_index_back,1).*...
                params_cond(plot_cell_index_back,6);

            para_synergy_antag = params_cond(plot_cell_index_cond,1).*...
                params_cond(plot_cell_index_cond,6);


            [f1, x1] = ksdensity(-log10(para_all_responder));

            f2 = ksdensity(-log10(para_synergy_antag),x1);

            figure

            fold_zoom_1 = 1.25;
            paperpos=[0,0,110,100];
            papersize=[110 100];
            draw_pos=[10,10,90,80];
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)


            plot(x1, f1, 'Color',[0.7,0.7,0.7],'LineWidth', 1); hold on;

            plot(x1, f2,'Color',[1,0,0], 'LineWidth', 1); hold on;
            xlabel(strcat('-log10(',params_names_label{1},'*',params_names_label{6},')'))
            ylabel('pdf')
            saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_',params_names_label{1},params_names_label{6},metric_name),'pdf')
            close()

            W_diff(1) = w_distance(-log10(para_all_responder),-log10(para_synergy_antag), 2)

        end

        if 1 % k68/k75
            para_all_responder = params_cond(plot_cell_index_back,7)./...
                params_cond(plot_cell_index_back,8);

            para_synergy_antag = params_cond(plot_cell_index_cond,7)./...
                params_cond(plot_cell_index_cond,8);

            [f1, x1] = ksdensity(log10(para_all_responder));

            f2 = ksdensity(log10(para_synergy_antag),x1);

            figure
            fold_zoom_1 = 1.25;
            paperpos=[0,0,110,100];
            papersize=[110 100];
            draw_pos=[10,10,90,80];
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)

            plot(x1, f1, 'Color',[0.7,0.7,0.7],'LineWidth', 1); hold on;
            plot(x1, f2,'Color',[1,0,0], 'LineWidth', 1); hold on;
            xlabel(strcat('log10(',params_names_label{7},'/',params_names_label{8},')'))
            ylabel('pdf')

            %saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_try',metric_name),'pdf')

            saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_vs_nonsyn_single_cell_para_',params_names_label{7},params_names_label{8},metric_name),'pdf')
            close()

            W_diff(2) = w_distance(log10(para_all_responder),log10(para_synergy_antag), 2)

        end


    end
end

if 1 % updated & tested : 12/01
    %% antag condition: under 5 by 5 condtion
    % add NFkB

    switch metric_name
        case'peak'
            i_row_plot_cond = 4;
            i_col_plot_cond = 5;
        case 'total'
            i_row_plot_cond = 5;
            i_col_plot_cond = 5;
    end

    i_cond = i_row_plot_cond * 6 + i_col_plot_cond + 1;
    params_cond = gene_info_all{i_cond}.gene_parameter_value_vec_genotype{1}{1}';
    params_names = gene_info_all{i_cond}.parameter_name_vec{1};

    % remove extra useless parameters: first p52n2, overwrite by second p52n2
    params_cond = params_cond(:,[1:2,4:end]);
    params_names = params_names([1:2,4:end]);

    % estimate NFkB for each cell
    max_cell_by_cond = -ones(1000,36);

    for i_cond_max = 1:36
        max_cell_by_cond(:,i_cond_max) = max(data.model_sim{i_cond_max}(1:12:end,:),[],2);
    end

    est_max_nuc_NFkB = max(max_cell_by_cond,[],2);



    params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);

    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
    responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

    % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
    antag_score_nfkb_cond = Antag_score(:,i_row_plot_cond,i_col_plot_cond);
    plot_cell_index_cond = find(antag_score_nfkb_cond<0.9 & responder_cells_nfkb_cond);
    plot_cell_index_back = find(antag_score_nfkb_cond>1.05 & responder_cells_nfkb_cond);


    if 0 % plot scatter plots and pdf
        %%
        fold_zoom = 13;
        paperpos=[0,0,110,100];
        papersize=[110 100];
        draw_pos=[10,10,100,90];

        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

        len_para = length(params_names);
        W_diff = zeros(1,len_para);
        for i_para = 1:len_para
            for j_para = i_para:len_para
                subplot(len_para,len_para,(i_para-1)* len_para+ j_para)
                if i_para == j_para
                    nbins = 10;
                    para_all_responder = params_cond(plot_cell_index_back,i_para);
                    para_synergy_antag = params_cond(plot_cell_index_cond,i_para);
                    %edges = linspace(min(para_all_responder),
                    %max(min(para_all_responder)), nbins+1); ,edges
                    % Estimate and plot the distribution for para_all_responder

                    [f1, x1] = ksdensity(log10(para_all_responder));
                    plot(x1, f1, 'Color',[0.7,0.7,0.7],'LineWidth', 2); hold on;

                    % Estimate and plot the distribution for para_synergy_antag
                    f2 = ksdensity(log10(para_synergy_antag),x1);
                    plot(x1, f2,'Color',[1,0,0], 'LineWidth', 2); hold on;
                    xlabel(strcat('log10(',params_names_label{i_para},')'))
                    ylabel('pdf')
                    % % Add labels and legend
                    % xlabel('Parameter Value');
                    % ylabel('Probability Density');
                    % legend('All Responders', 'Synergy/Antagonism');
                    % title('Smoothed Distribution Curves');

                    % histogram(para_all_responder,'Normalization', 'pdf');hold on
                    % histogram(para_synergy_antag,'Normalization', 'pdf');hold on
                    % set(gca,'XScale','log')
                    %% synergistic cells distribution vs default distribution

                    W_diff(i_para) = w_distance(log10(para_all_responder),log10(para_synergy_antag), 2);

                else
                    scatter(log10(params_cond(plot_cell_index_back,j_para)),log10(params_cond(plot_cell_index_back,i_para)),16,[0.7,0.7,0.7],'filled');hold on
                    scatter(log10(params_cond(plot_cell_index_cond,j_para)),log10(params_cond(plot_cell_index_cond,i_para)),16,[1,0,0],'filled')
                    % set(gca,'XScale','log','YScale','log')
                    % set(gca,'XTickLabel',{},'YTickLabel',{})
                    xlabel(strcat('log10(',params_names_label{j_para},')'))
                    ylabel(strcat('log10(',params_names_label{i_para},')'))
                end

            end
        end

        saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_vs_nonant_single_cell_para_distrib_',metric_name),'pdf')
        close()

        figure(2)
        set(gcf, 'PaperUnits','points')
        fold_zoom_fig2 = 1.5;
        set(gcf, 'PaperPosition', paperpos*fold_zoom_fig2,'PaperSize', papersize*fold_zoom_fig2,'Position',draw_pos*fold_zoom_fig2)

        bar(W_diff)
        set(gca,'XTickLabel',params_names_label)
        ylabel('W-Distance')
        %ylim([0,0.7])
        saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_vs_nonant_single_cell_para_WDist_',metric_name),'pdf')
        close()

    end

    if 1 % plot single parameter distribution
%%
        plot_para_index = [5,7,8];

        for i_para = 1:length(plot_para_index)

            index_para = plot_para_index(i_para);

            para_all_responder = params_cond(plot_cell_index_back,index_para);
            para_synergy_antag = params_cond(plot_cell_index_cond,index_para);

            if 0 %  One-sided Mann–Whitney U tests
                params_names{index_para}
                'left'
                [p1, ~, ~] = ranksum(para_all_responder, para_synergy_antag, 'tail','left')
                'right'
                [p1, ~, ~] = ranksum(para_all_responder, para_synergy_antag, 'tail','right')

            end

            [f1, x1] = ksdensity(log10(para_all_responder));

            f2 = ksdensity(log10(para_synergy_antag),x1);

            figure

            fold_zoom_1 = 0.8;
            paperpos=[0,0,110,100];
            papersize=[110 100];
            draw_pos=[10,10,90,80];
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)

            plot(x1, f1, 'Color',[0.7,0.7,0.7],'LineWidth', 1); hold on;

            plot(x1, f2,'Color',[1,0,0], 'LineWidth', 1); hold on;
            %xlabel(strcat('-log10(',params_names_label{1},'*',params_names_label{6},')'))
            %ylabel('pdf')
            saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_vs_nonantag_single_cell_para_',params_names_label{index_para},metric_name,'_withlabel'),'pdf')
            set(gca,'XTickLabel',{},'YTickLabel',{})
            saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_vs_nonantag_single_cell_para_',params_names_label{index_para},metric_name),'pdf')
            

            close()
        end
    end
end