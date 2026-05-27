function [W_diff_mat] = draw_Wdist_mat(collect_feature_vects,Wdiff_savefile_name,plot_opt)
% to move to lib, if decided to use for publication
if nargin<2
    Wdiff_savefile_name = '';
    % 'W_diff_mat.mat'
end

if nargin< 3
    plot_opt.cal_Wdis_mat = 1;
    plot_opt.save_Wdis_mat = 0;
    plot_opt.plot_wdis_mat = 0;
    plot_opt.plot_hierach_tree = 0;
    plot_opt.fig_save_path = './';
    plot_opt.fig_save_name_wdis_mat = 'Wdis_31_cond';
end


if plot_opt.cal_Wdis_mat
    w_dis = [];

    codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

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

    if plot_opt.save_Wdis_mat
        save(strcat(plot_opt.data_save_file_path,Wdiff_savefile_name),"W_diff_mat");
    end
else
    load(strcat(plot_opt.data_save_file_path,Wdiff_savefile_name))
end

%% Figure xx:
mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
    ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

if plot_opt.plot_wdis_mat
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

    saveas(gcf,strcat(plot_opt.fig_save_path,plot_opt.fig_save_name_wdis_mat),'pdf')
    close()

end

if plot_opt.plot_hierach_tree % hierarchicle cluster
    %
    y = [];
    for i_row = 1:size(W_diff_mat,1)
        y = [y,W_diff_mat(i_row,(i_row+1) : size(W_diff_mat,1))];
    end
    % y = y';
    Z = linkage(y,'ward');%'single'
    labels_all = cellfun(@replace, ...
        data_info.info_ligand, ...
        repmat({'_'}, size(data_info.info_ligand)), ...
        repmat({'-'}, size(data_info.info_ligand)), 'UniformOutput', false);
    %Z(30,3) = 3;
    dendrogram(Z,0,'Labels',labels_all,'Orientation','left')
    %dendrogram(Z,0,'Orientation','left')
    %dendrogram(Z,'Orientation','left');
    saveas(gcf,strcat(plot_opt.fig_save_path,'cond_hiearchical_cluster_tree'),'pdf')
    close()

end