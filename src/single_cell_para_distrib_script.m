if 1
%% synergy condition: under 5 by 5 condtion

fold_zoom = 13;
paperpos=[0,0,110,100];
papersize=[110 100];
draw_pos=[10,10,100,90];

figure(1)
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)


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
plot_cell_index_cond = find(synergy_score_nfkb_cond>1.25 & responder_cells_nfkb_cond);

len_para = length(params_names);
W_diff = zeros(1,len_para);

figure(1)
for i_para = 1:len_para
    for j_para = i_para:len_para

        subplot(len_para,len_para,(i_para-1)* len_para+ j_para)
        if i_para == j_para
            nbins = 10;
            para_all_responder = params_cond(responder_cells_nfkb_cond,i_para);
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
            scatter(log10(params_cond(responder_cells_nfkb_cond,j_para)),log10(params_cond(responder_cells_nfkb_cond,i_para)),16,[0.7,0.7,0.7],'filled');hold on
            scatter(log10(params_cond(plot_cell_index_cond,j_para)),log10(params_cond(plot_cell_index_cond,i_para)),16,[1,0,0],'filled')
            %set(gca,'XScale','log','YScale','log')
            %set(gca,'XTickLabel',{},'YTickLabel',{})
            xlabel(strcat('log10(',params_names_label{j_para},')'))
            ylabel(strcat('log10(',params_names_label{i_para},')'))

        end

    end
end

figure(1)
saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_single_cell_para_distrib_',metric_name),'pdf')
close()

figure(2)
set(gcf, 'PaperUnits','points')
fold_zoom_fig2 = 1.5;
set(gcf, 'PaperPosition', paperpos*fold_zoom_fig2,'PaperSize', papersize*fold_zoom_fig2,'Position',draw_pos*fold_zoom_fig2)

bar(W_diff)
set(gca,'XTickLabel',params_names_label)
ylabel('W-Distance')
%ylim([0,0.7])
saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_single_cell_para_WDist_',metric_name),'pdf')
close()
end

if 1
%% antag condition: under 5 by 5 condtion
fold_zoom = 13;
paperpos=[0,0,110,100];
papersize=[110 100];
draw_pos=[10,10,100,90];

figure(1)
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

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

params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);

%synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

% score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
antag_score_nfkb_cond = Antag_score(:,i_row_plot_cond,i_col_plot_cond);
plot_cell_index_cond = find(antag_score_nfkb_cond<0.9 & responder_cells_nfkb_cond);


len_para = length(params_names);
W_diff = zeros(1,len_para);
for i_para = 1:len_para
    for j_para = i_para:len_para
        subplot(len_para,len_para,(i_para-1)* len_para+ j_para)
        if i_para == j_para
            nbins = 10;
            para_all_responder = params_cond(responder_cells_nfkb_cond,i_para);
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
            scatter(log10(params_cond(responder_cells_nfkb_cond,j_para)),log10(params_cond(responder_cells_nfkb_cond,i_para)),16,[0.7,0.7,0.7],'filled');hold on
            scatter(log10(params_cond(plot_cell_index_cond,j_para)),log10(params_cond(plot_cell_index_cond,i_para)),16,[1,0,0],'filled')
            % set(gca,'XScale','log','YScale','log')
            % set(gca,'XTickLabel',{},'YTickLabel',{})
            xlabel(strcat('log10(',params_names_label{j_para},')'))
            ylabel(strcat('log10(',params_names_label{i_para},')'))
        end

    end
end

saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_single_cell_para_distrib_',metric_name),'pdf')
close()

figure(2)
set(gcf, 'PaperUnits','points')
fold_zoom_fig2 = 1.5;
set(gcf, 'PaperPosition', paperpos*fold_zoom_fig2,'PaperSize', papersize*fold_zoom_fig2,'Position',draw_pos*fold_zoom_fig2)

bar(W_diff)
set(gca,'XTickLabel',params_names_label)
ylabel('W-Distance')
%ylim([0,0.7])
saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_single_cell_para_WDist_',metric_name),'pdf')
close()
end