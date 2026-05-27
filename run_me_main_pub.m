addpath('./lib/')
addpath('./src/')
addpath('./bin/')
data_save_file_path = '/Users/xiaoluguo/Documents/Postdoc_projects/Projects/Ongoing_combinatorial_ligand/revision/data_upload_cowork/';
fig_save_path = '../subfigures/';
addpath(data_save_file_path)


%% Figure 1 & S1: responder ratio,  run on 11/18/2025
if 0
    if 0 % draw sim [tested 05/24/2026]
    % updated using PIC-TLR3 deg increased model
    % save_metric_name = 'Sim2025_all_comb_SR.mat';
    save_metric_name = 'Sim2025_all_comb_CpG_PIC_no_compete.mat';

    draw_responder_ratio_combinatorial_ligand_1030

    end

    if 0 % draw exp [tested 05/24/2026]
 
        draw_responder_ratio_combinatorial_ligand_exp
    end

end

%% Figure S1 & 2AB & S4: plot heatmaps of trajectories, for updated PIC-TLR3 deg increase: run on 11/18/2025
if 0
    %% plot 31 conditions for SR
    if 0 %Figure 2E
        save_metric_name = 'Sim2025_all_comb_SR.mat';
        % fig_save_path = '/Users/xiaoluguo/Documents/Postdoc_projects/Projects/Ongoing_combinatorial_ligand/subfigures2025/draw_heatmaps_revision/';

        load(strcat(data_save_file_path,save_metric_name))

        save_name_tag = 'draw'
        for i_traj = 1:length(data.model_sim)
            data_traj = data.model_sim{i_traj}(1:9:end,:);
            fig_gcf = plot_heatmap(data_traj);
            saveas(gcf,strcat(fig_save_path,save_name_tag,'_',data.info_ligand{i_traj}),'pdf');
            close
        end
    end

    %% plot 31 conditions for non-CpG-pIC conditions

    if 0 % Figure 2A, updated: 11/18/2025
        save_metric_name = 'Sim2025_all_comb_CpG_PIC_no_compete.mat';

        load(strcat(data_save_file_path,save_metric_name))

        save_name_tag = 'draw'
        for i_traj = 1:length(data.model_sim)
            data_traj = data.model_sim{i_traj}(1:9:end,:);
            fig_gcf = plot_heatmap(data_traj);
            saveas(gcf,strcat(fig_save_path,save_name_tag,'_',data.info_ligand{i_traj}),'pdf');
            close
        end
    end

    %% Figure 2B: plot exp data
    if 0
        % draw_heatmaps_exp_supriya
        %fig_save_path = '../subfigures2025/revision/figureS1/heatmaps/';
        %data_save_file_path_sim = '../raw_data2025_revision/';
        load(strcat(data_save_file_path,'All_exp_sim_codon_dual_ligand_202511_nonminmaxscaled.mat'))
        save_name_tag = 'draw_exp_';

        index_exp_cond = [1,3,4,5,13,17,18]; % condtions: TNF,CpG,LPS,pIC TNF_LPS, CpG_LPS, CpG_pIC
        for i_cond = 1:length(index_exp_cond)
            data_traj = metrics{index_exp_cond(i_cond)}.time_series;
            fig_gcf = plot_heatmap(data_traj);
            saveas(gcf,strcat(fig_save_path,save_name_tag,'_',data_info.info_ligand{index_exp_cond(i_cond)}),'pdf');
            close
        end

end
end



%% Figure 2CF &S3 : tested 05/26/2026
if 0   
    compare_sample_supriya_202605
end
%% Figure 2G: draw duration total diff: successfully run on 05/24/2026
if 0
    % scatter_diff_exp_old_new(fig_save_path_diff)
    bar_diff_exp_old_20251118(fig_save_path)
end




%% Figure S1 & S4 & 3A & 3B & 4E & 5E: codon distribution for all comb; W-dist & hierachicle clusters run on 11/18/2025
if 0

    if 0 % run CpG-pIC non-compete model %Figure 4E, Figure EV1
        save_metric_name = 'Sim2025_all_comb_CpG_PIC_no_compete.mat'; % no CpG-pIC compete, with pIC-TLR deg increased
        compete_or_not = 0;
        codon_distrib_comb(data_save_file_path,save_metric_name,fig_save_path,compete_or_not)

    end

    if 0 % run CpG-pIC compete model %Figure 3A, Figure 4E, Figure EV2
        save_metric_name = 'Sim2025_all_comb_SR.mat'; % for CpG-pIC compete, with pIC-TLR deg increased
        compete_or_not = 1;
        codon_distrib_comb(data_save_file_path,save_metric_name,fig_save_path,compete_or_not)

    end

    if 0 % Figure 5E: plot the histgram of Wdis bewteen CpG-pIC and other conditions 
        % for compete vs non-compete model
        plot_hist_Wdis_compete_vs_non_compete

    end

end

%% Figure 3 related: save the experimental/simulation dataset to ML format for training the classifier
if 0 % exp data
    save_metric_name = 'Exp_data_metrics_2024.mat';
    codon_save_path = data_save_file_path;
    codon_save_to_ML_format_202506(data_save_file_path,save_metric_name,codon_save_path);
end

if 0 % simulation data
    % todo:
    % need to rerun, using the updated cal_codon_metric_20251118
    % need to recalculte all the classifier results
    % check what's to_check_codon_save_to_ML_format_sim_202510.m for 
    save_metric_name = 'Sim2025_all_comb_SR.mat'; % for CpG-pIC compete, with pIC-TLR deg increased

    codon_save_path = data_save_file_path;

    index_to_save = mat2cell(1:31,[1],ones(1,31));
    % codon_save_to_ML_format_202506(data_save_file_path,save_metric_name,codon_save_path,index_to_save);%  to delete

    codon_save_to_ML_format_202604(data_save_file_path,save_metric_name,codon_save_path,index_to_save);
    % need to rename the file (to delete _exp) and move to the python codes
    % folder
    a = 1;
    % traj_save_to_ML_format(data_save_file_path,save_metric_name);
end



%% Figure4 & Figure 5 & supp figs: synergy & antagnoism: run on 12/01;
if 0 % new version 0526/2026
    single_cell_synergy

end

%% plot heatmaps Figure S8
if 0
    save_fig_path = fig_save_path;

    all_pairs = {{'CpG','PolyIC'},... 1
        {'LPS','PolyIC'},... 2
        {'TNF','PolyIC'},... 3
        {'TNF','LPS'},... 4
        {'LPS','CpG'},... 5
        {'CpG','Pam3CSK'},... 6
        {'Pam3CSK','PolyIC'},... 7
        {'LPS','Pam3CSK'},... 8
        {'TNF','CpG'},... 9
        {'TNF','Pam3CSK'}}; %10

    total_pair = 10;
    for i_pair = 8%1:total_pair

        %i_pair = 7;
        ligand1 = all_pairs{i_pair}{1};
        ligand2 = all_pairs{i_pair}{2};

        % Sim0620_CpG_PolyIC_sc_signergy_different_dose
        save_metric_name = strcat('Sim1008_',ligand1,'_',ligand2,'_sc_synergy_different_dose.mat');
        load(strcat(data_save_file_path ,save_metric_name))


        for i_species = 1:6
            %i_species = 5

            % calculate y limits for each speicies
            i_cond = 1;
            data_traj_all = [];
            for dose_ligand1 = 1:6
                for dose_ligand2 = 1:6
                    data_traj_all = [data_traj_all;data.model_sim{i_cond}(i_species:12:end,:)];
                    i_cond = i_cond+1;

                    

                end
            end
            ylim_species = [0,0];
            ylim_species(1) = prctile(data_traj_all(:),0.25);
            ylim_species(2) = prctile(data_traj_all(:),99.75);
            i_species
            ylim_species


            % draw the heatmap
            i_cond = 1;
            for dose_ligand1 = 1:6
                for dose_ligand2 = 1:6
                    species_name = data_info.species_outputname{i_species};
                    save_name_file = strcat('heatmap_',ligand1,'D',num2str(dose_ligand1-1),'_',ligand2,'D',num2str(dose_ligand2-1),'_',species_name);
                    data_traj = data.model_sim{i_cond}(i_species:12:end,:);
                    fig_gcf = plot_heatmap(data_traj,ylim_species);
                    saveas(gcf,strcat(save_fig_path,save_name_file),'pdf');
                    close
                    i_cond = i_cond+1;
                end
            end
        end

    end
end





%% Figure4DEF plot wdist for Hill =1, Hill =2, etc
if 0 % Figure 4D
    Fig_compare_synergy_hill1_hill2
    
end

if 0
    if 0 %Figure S11A
        draw_heatmaps_hill12
    end

    if 0 %Figure 4EF
        
        fig_save_path_sub = fig_save_path;
        draw_W_dist_main_seperate_cal_codon
    end

    if 0 %Figure S11E
        % plot the integral distribution for synergy test
        plot_integral_hist_hill1_hill2
    end
end




%% [run sim]
% run all codes under ./run_on_server/ , to get all related simulation
% results


%% subfuncitons for plot heatmaps
function fig_gcf = plot_heatmap(data_traj,ylim_vec)
% be default plot by integral

% data_traj = metrics{1, 1};
if sum(isnan(data_traj  ), 'all')
    data_traj_to_sort=  data_traj;
    data_traj_to_sort(isnan(data_traj_to_sort)) = 0;
    [~,data_order_index] = sort(sum(data_traj_to_sort,2),"descend");
else
    [~,data_order_index] = sort(sum(data_traj,2),"descend");
end


figure()

if nargin ==1
    color_limits = [-0.001,0.3];
elseif nargin ==2
    color_limits = ylim_vec;
else
    error('too many/too few inputs')
end
paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

cell_num=size(data_traj,1);
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

% subplot(1,length(vis_data_field),i_data_field)
h=heatmap(data_traj(data_order_index,:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limits);%[-0.001,0.2] for TNF

XLabels = 0:5:((size(data_traj,2)-1)*5);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/60);
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:cell_num;
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.YDisplayLabels = YCustomXLabels;

% xlabel('Time (hours)');
% ylabel(vis_data_field{i_data_field});
% clb=colorbar;
% clb.Label.String = 'NFkB(A.U.)';
colorbar('off')

set(gca,'fontsize',14,'fontname','Arial');
fig_gcf = gcf;

end