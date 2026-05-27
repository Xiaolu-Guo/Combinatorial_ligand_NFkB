function [] = codon_distrib_comb(data_save_file_path,data_file,fig_save_path,compete_or_not)%,save_name_tag

if compete_or_not
    figure_save_name = 'compete';
else
    figure_save_name = 'nocompete';
end

load(strcat(data_save_file_path,data_file))
codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

field_list = fieldnames(metrics{1});
data_info.info_ligand = data.info_ligand;
data_info.info_dose_str = data.info_dose_str;

metric_cal = cell(1,2);
for i_id = 1:length(metrics)
    data_info.data_label{i_id} = 'sim_WT';
    for i_field = 1:length(field_list)

        metric_cal{i_id}.(field_list{i_field}) = metrics{i_id}.(field_list{i_field})(1:9:end,:);
    end
end
collect_feature_vects = calculate_codon_from_metric20251118(data_info,metric_cal);

if 1
    for i_codon = 1:length(codon_list)
        figure(1)
        paperpos=[0,0,1500,100];
        papersize=[1500 100];
        draw_pos=[10,10,1480,90];
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

        al_goodplot_pair_RMSD_diff_size(y,[1:31],0.5,ones(size(y,2),1)*[255 0 0]/255 ,'bilateral',[],std_y); %left
        % al_goodplot_pair_RMSD_diff_size(z,[2:3:5],0.5,ones(size(z,2),1)*[0 0 0]/255,'bilateral',[],std_y);


        xlim([0 32]);

        xticks([0:32]);
        xticklabels({});
        %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})

        ylim([-2.1,2.1]);
        % for i_x = 1:10
        %     plot([i_x,i_x],[0,5],'--','Color','k');hold on
        % end
        set(gca,'fontsize',14,'fontname','Arial');
        %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');

        saveas(gcf,strcat(fig_save_path,'Codon_',codon_list{i_codon}),'pdf');%,save_name_tag
        close

    end
end

if 1 % plot heatmap differences across different conditions
    %
    if 1
        w_dis = [];

        for i_codon = 1:length(codon_list)
            for i_sti_1 = 1:length(collect_feature_vects.Duration)
                for i_sti_2 = 1:length(collect_feature_vects.Duration)

                    W_diff(i_sti_1,i_sti_2,i_codon) =... % for diagnal, differences between exp and simulaiton conditions
                        w_distance(collect_feature_vects.(codon_list{i_codon}){i_sti_1},...
                        collect_feature_vects.(codon_list{i_codon}){i_sti_2}, 2);


                end
            end
        end

        W_diff_mat = mean(W_diff,3);

        index_reorder = 1:length(collect_feature_vects.Duration);

        W_diff_mat = W_diff_mat(index_reorder,index_reorder);

        save(strcat(data_save_file_path,'W_diff_mat_',figure_save_name,'.mat'),"W_diff_mat");

    else

        load(strcat(data_save_file_path,'W_diff_mat_',figure_save_name,'.mat'));

    end



    % if 1
    %     codon_list = {'Speed','PeakAmplitude','TotalActivity','EarlyVsLate','OscVsNonOsc'};
    %     w_dis = [];
    %
    %     for i_codon = 1:length(codon_list)
    %         for i_sti_1 = 1:length(collect_feature_vects.Duration)
    %             for i_sti_2 = 1:length(collect_feature_vects.Duration)
    %
    %                 W_diff(i_sti_1,i_sti_2,i_codon) =... % for diagnal, differences between exp and simulaiton conditions
    %                     w_distance(collect_feature_vects.(codon_list{i_codon}){i_sti_1},...
    %                     collect_feature_vects.(codon_list{i_codon}){i_sti_2}, 2);
    %
    %
    %             end
    %         end
    %     end
    %
    %     W_diff_mat = mean(W_diff,3);
    %
    %     index_reorder = 1:length(collect_feature_vects.Duration);
    %
    %     W_diff_mat = W_diff_mat(index_reorder,index_reorder);
    %     save(strcat(data_save_file_path,'W_diff_mat_5codon.mat'),"W_diff_mat");
    % else
    %     load(strcat(data_save_file_path,'W_diff_mat_5codon.mat'))
    % end

    %% Figure xx:
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

    if 1
        figure (1)
        set(gcf, 'PaperUnits','points')

        paper_pos = [0,0,450,450];
        paper_size = [450,450];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]


        heatmap(W_diff_mat,'colormap',mymap,'ColorbarVisible','off')%,'ColorbarVisible','off'
        caxis([0,1])

        % Getting handle of the heatmap
        h = gca;

        % Removing tick labels from x and y axes
        h.XDisplayLabels = repmat({''}, size(h.XDisplayData));
        h.YDisplayLabels = repmat({''}, size(h.YDisplayData));

        saveas(gcf,strcat(fig_save_path,'Wdis_31_cond_',figure_save_name),'pdf')
        close()

    end

    if 0 % bar plot of average distance
        rank_dist = sum(W_diff_mat,1)/30;

        figure (1)
        set(gcf, 'PaperUnits','points')

        paper_pos = [0,0,80,80];
        paper_size = [80,80];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]


        bar(rank_dist(1:5))
        %heatmap(W_diff_mat,'colormap',mymap,'ColorbarVisible','off')%,'ColorbarVisible','off'
        ylim([0,1.5])

        saveas(gcf,strcat(fig_save_path,'Rank_average_Wdis_single_ligand'),'pdf')
        close()

    end


end


if 1 % hierarchicle cluster
    %

    load(strcat(data_save_file_path,'W_diff_mat_',figure_save_name,'.mat'))
    y = [];
    for i_row = 1:size(W_diff_mat,1)
        y = [y,W_diff_mat(i_row,(i_row+1) : size(W_diff_mat,1))];
    end
    % y = y';
    % Z = linkage(y,'ward');%'single'
    % labels_all = cellfun(@replace, ...
    %     data_info.info_ligand, ...
    %     repmat({'_'}, size(data_info.info_ligand)), ...
    %     repmat({'-'}, size(data_info.info_ligand)), 'UniformOutput', false);
    % %Z(30,3) = 3;
    % dendrogram(Z,0,'Labels',labels_all,'Orientation','left')
    %dendrogram(Z,0,'Orientation','left')
    %dendrogram(Z,'Orientation','left');

    Z = linkage(y,'average');

% 'average'	
% 'centroid'	
% 'complete'	
% 'median'	
% 'single'	
% 'ward'	
% 'weighted'

    % Get default order
    [H, T, perm] = dendrogram(Z, 0);
    % Swap first two branches(they are under the same cluster)
    perm_new = perm;
    perm_new([1, 2]) = perm_new([2, 1]);

    labels_all = cellfun(@replace, ...
        data_info.info_ligand, ...
        repmat({'_'}, size(data_info.info_ligand)), ...
        repmat({'-'}, size(data_info.info_ligand)), 'UniformOutput', false);

    labels_all{1} = labels_all{4};
    labels_all{4} = 'TNF';
    dendrogram(Z,0,'Reorder', perm_new,'Labels',labels_all,'Orientation','left')

    set(gca,'XTick',[0,0.2,0.4,0.6,0.8,1.0,1.2])
    saveas(gcf,strcat(fig_save_path,'cond_hiearchical_cluster_tree_',figure_save_name),'pdf')
    close()

end

end

% function zscores = zscore_codon(signaling_codons,interval_cell)
% if nargin <2
%     interval_cell = 1;
% end
% all_codons = [];
% for i = 1:length(signaling_codons)
%     all_codons = [all_codons;signaling_codons{i}(1:interval_cell:end,:)];
% end
% mean_zscore = mean(all_codons);
% std_zscore = std(all_codons);
% zscores = cell(length(signaling_codons),1);
% for i = 1:length(signaling_codons)
%     zscores{i} = (signaling_codons{i}(1:interval_cell:end,:)-mean_zscore)/std_zscore;
% end
% end