
% delete / or move to debug: the original codes in run_me_main.m


%% user define:
% decide whether run synergy score or antagnism score:
synergy_1_antag_0 = 1;

% decide whether synergy cell ratio devided by responder cells only
devided_responder_or_all_cells = 1;

% decide whether plot the heatmaps for scores
plot_score_heatmap = 0; % set i_pair to 1:total_pair

% decide whether plot the bar for synergy/antag ratio under different hill
plot_bar_synergy_ratio_diff_hill = 1;

% decide whether plot the histogram for the score distribution for given
% pair
plot_histogram = 0; % set i_pair to 8

% decide whether plot single-cell trajectories
plot_single_cell = 0; % set i_pair to 8

% decide whether plot parameter distribution
plot_para_distrib = 0;

% decide whether run and plot the logistric regression classifier using the
% parameter to predict synergy/antagnism
plot_para_score_logistic_reg_classifier = 0;

% decide whether run global parameter sensitivity analysis
global_sens_analysis = 0;

% decide whether hill = 1 or hill = 2
hill_coeff = 1;

% fig save paths and data save paths
fig_save_path_sub = fig_save_path;

%% initialize


all_pairs = {{'TNF','LPS'},... 1
    {'TNF','CpG'},... 2
    {'TNF','PolyIC'},... 3
    {'TNF','Pam3CSK'}, ... 4
    {'LPS','CpG'},... 5
    {'LPS','PolyIC'},... 6
    {'LPS','Pam3CSK'},... 7
    {'CpG','PolyIC'},... 8
    {'CpG','Pam3CSK'},... 9
    {'Pam3CSK','PolyIC'}};% 10

%%  run and plot

% different orders from figure 3, remember to update
% all_pairs = {{'CpG','PolyIC'},... 1
%     {'LPS','PolyIC'},... 2
%     {'TNF','PolyIC'},... 3
%     {'TNF','LPS'},... 4
%     {'LPS','CpG'},... 5
%     {'CpG','Pam3CSK'},... 6
%     {'Pam3CSK','PolyIC'},... 7
%     {'LPS','Pam3CSK'},... 8
%     {'TNF','CpG'},... 9
%     {'TNF','Pam3CSK'}}; %10
total_pair = 10;
total_species = 6;
ligand1_dose_num = 1;
ligand2_dose_num = 1;
total_dose_ligand1 = 1;
total_dose_ligand2 = 1;

metric_name_all = {'peak','total'};

for hill_coeff = 1:2

    if hill_coeff == 1
        data_file_list = {'Sim1008_TNF_LPS_sc_dose2_IKK_hill1.mat'
            'Sim1008_TNF_CpG_sc_dose2_IKK_hill1.mat'
            'Sim1008_TNF_PolyIC_sc_dose2_IKK_hill1.mat'
            'Sim1008_TNF_Pam3CSK_sc_dose2_IKK_hill1.mat'
            'Sim1008_LPS_CpG_sc_dose2_IKK_hill1.mat'
            'Sim1008_LPS_PolyIC_sc_dose2_IKK_hill1.mat'
            'Sim1008_LPS_Pam3CSK_sc_dose2_IKK_hill1.mat'
            'Sim1008_CpG_PolyIC_sc_dose2_IKK_hill1.mat'
            'Sim1008_CpG_Pam3CSK_sc_dose2_IKK_hill1.mat'
            'Sim1008_Pam3CSK_PolyIC_sc_dose2_IKK_hill1.mat'};
    end

    if hill_coeff == 2
        data_file_list = {'Sim1008_TNF_LPS_sc_dose2_IKK_hill2.mat'
            'Sim1008_TNF_CpG_sc_dose2_IKK_hill2.mat'
            'Sim1008_TNF_PolyIC_sc_dose2_IKK_hill2.mat'
            'Sim1008_TNF_Pam3CSK_sc_dose2_IKK_hill2.mat'
            'Sim1008_LPS_CpG_sc_dose2_IKK_hill2.mat'
            'Sim1008_LPS_PolyIC_sc_dose2_IKK_hill2.mat'
            'Sim1008_LPS_Pam3CSK_sc_dose2_IKK_hill2.mat'
            'Sim1008_CpG_PolyIC_sc_dose2_IKK_hill2.mat'
            'Sim1008_CpG_Pam3CSK_sc_dose2_IKK_hill2.mat'
            'Sim1008_Pam3CSK_PolyIC_sc_dose2_IKK_hill2.mat'};
    end

    for i_metric_name = 1% 1:length(metric_name_all)
        metric_name = metric_name_all{i_metric_name};

        for i_pair = 1:total_pair % 8 %
            %i_pair = 7;
            ligand1 = all_pairs{i_pair}{1};
            ligand2 = all_pairs{i_pair}{2};

            load(strcat(data_save_file_path,data_file_list{i_pair}))
            %
            % % Sim0620_CpG_PolyIC_sc_signergy_different_dose
            % save_metric_name = strcat('Sim1008_',ligand1,'_',ligand2,'_sc_synergy_different_dose.mat');
            % load(strcat(data_save_file_path ,save_metric_name))

            synergy_ratio = -ones(total_species,ligand1_dose_num,ligand2_dose_num);
            antag_ratio = -ones(total_species,ligand1_dose_num,ligand2_dose_num);
            for i_species = 1% :total_species
                %i_species = 5
                species_name = data_info.species_outputname{i_species};

                % calculate the metric (peak or integral) for each single cell
                [tensor_sc_metric,vis_index_ligand1,vis_index_ligand2,vis_sti] = calculate_sc_metric(data_info,data,species_name,metric_name);

                % initialize tensors for saving values
                Synergy_score = nan(1000,ligand1_dose_num,ligand2_dose_num);
                Antag_score = nan(1000,ligand1_dose_num,ligand2_dose_num);
                Valid_SNR = false(1000,ligand1_dose_num,ligand2_dose_num);
                Valid_SNR_all_cond = false(1000,ligand1_dose_num+1,ligand2_dose_num+1);

                % calculate synergy score and antagonism score for each cell
                for i_cell = 1:size(tensor_sc_metric,1)
                    cell_mat = squeeze(tensor_sc_metric(i_cell,:,:));

                    % responder definition should be revised for integral
                    % metric: 5 fold of integral compared to basal is too
                    % strong for responder def
                    % this works for peak metric; thus focus on Peak
                    Valid_SNR(i_cell,:,:) = (cell_mat(2:end,2:end) ) > cell_mat(1,1)*5;
                    Valid_SNR_all_cond(i_cell,:,:) = (cell_mat(1:end,1:end) ) > cell_mat(1,1)*5;

                    cell_additive_baseline = ((cell_mat(2:end,1) * ones(1,ligand2_dose_num) - cell_mat(1,1)) + ...
                        (ones(ligand1_dose_num,1) * cell_mat(1,2:end)  - cell_mat(1,1)));
                    Synergy_score(i_cell,:,:) = (cell_mat(2:end,2:end) - cell_mat(1,1)) ./ cell_additive_baseline;

                    cell_max_baseline = max((cell_mat(2:end,1) * ones(1,ligand2_dose_num) - cell_mat(1,1)),(ones(ligand1_dose_num,1) * cell_mat(1,2:end)  - cell_mat(1,1)));
                    Antag_score(i_cell,:,:) = (cell_mat(2:end,2:end) - cell_mat(1,1)) ./ cell_max_baseline;

                end

                if 1
                    for i_row = 1:total_dose_ligand1
                        for i_col = 1:total_dose_ligand2

                            % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
                            score_to_plot = Synergy_score(Valid_SNR_all_cond(:,i_row+1,i_col+1)...
                                |Valid_SNR_all_cond(:,i_row+1,1)|Valid_SNR_all_cond(:,1,i_col+1),i_row,i_col);

                            if devided_responder_or_all_cells

                                % responder_cells = sum(Valid_SNR(:,i_row,i_col));
                                responder_cells = sum(Valid_SNR_all_cond(:,i_row+1,i_col+1)...
                                    |Valid_SNR_all_cond(:,i_row+1,1)|Valid_SNR_all_cond(:,1,i_col+1));
                                if responder_cells == 0
                                    synergy_ratio(i_species,i_row,i_col) = -0.3;
                                else
                                    synergy_ratio(i_species,i_row,i_col) = sum(score_to_plot>1.25)/responder_cells;
                                end
                            else
                                synergy_ratio(i_species,i_row,i_col) = sum(score_to_plot>1.25)/1000;
                            end

                            if plot_histogram
                                figure(2)
                                subplot(5,5,(i_row-1) * 5 + i_col)
                                h1 = histogram(score_to_plot,[0.5:0.1:2]);
                                xlim([0.5,2]);
                                xline(1.25, '--r', 'LineWidth', 2);

                            end

                            % score_to_plot = Antag_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
                            score_to_plot = Antag_score(Valid_SNR_all_cond(:,i_row+1,i_col+1)...
                                |Valid_SNR_all_cond(:,i_row+1,1)|Valid_SNR_all_cond(:,1,i_col+1),i_row,i_col);

                            if devided_responder_or_all_cells
                                % responder_cells = sum(Valid_SNR(:,i_row,i_col));
                                responder_cells = sum(Valid_SNR_all_cond(:,i_row+1,i_col+1)...
                                    |Valid_SNR_all_cond(:,i_row+1,1)|Valid_SNR_all_cond(:,1,i_col+1));
                                if responder_cells == 0
                                    antag_ratio(i_species,i_row,i_col) = -0.3;
                                else
                                    antag_ratio(i_species,i_row,i_col) = sum(score_to_plot<0.9)/responder_cells;
                                end
                            else
                                antag_ratio(i_species,i_row,i_col) = sum(score_to_plot<0.9)/1000;
                            end
                            if plot_histogram
                                figure(3)
                                subplot(5,5,(i_row-1) * 5 + i_col)
                                h2 = histogram(score_to_plot,[0.7:0.025:1.2]);
                                xlim([0.7,1.2]);
                                xline(0.9, '--r', 'LineWidth', 2);
                            end

                        end
                    end
                end

                if plot_histogram

                    figure(2)
                    saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_synergy_hist_plot_','_',metric_name),'pdf')
                    close()
                    figure(3)
                    saveas(gcf,strcat(save_fig_path_sub,ligand1,ligand2,'_antag_hist_plot_','_',metric_name),'pdf')
                    close()
                end

                if plot_single_cell
                    % when run, pause on save, and adjuct the figure size
                    % before saving the figures.
                    if 1
                        % update name, move to /src
                        single_cell_plot_script
                    end

                    if 1
                        single_cell_family_plot_script
                    end
                    %plot_schematic_Step3
                end

                if plot_para_distrib
                    % can run till this point. then stop, go to the script,
                    % and change the parameters or figure setting as you
                    % want, and can be run directely without rerun this
                    % main function.
                    if 0 % plot synergy vs all responder (same for antag)
                        single_cell_para_distrib_script

                    end

                    if 1 % plot synergy vs non-synergy responder  (same for antag)
                        %single_cell_para_distrib_vs_non_syn_antag_script
                        single_cell_para_distrib_vs_non_syn_antag_script_0806
                    end
                end

                if global_sens_analysis

                    score_para_PAWN_analysis
                end

                if plot_para_score_logistic_reg_classifier
                    % can run till this point. then stop, go to the script,
                    % and change the parameters or figure setting as you
                    % want, and can be run directely without rerun this
                    % main function.
                    single_cell_para_score_logistic_reg_classifier_script
                end



                if synergy_1_antag_0
                    ratio_plot = squeeze(synergy_ratio(i_species,:,:));
                    synergy_ratio_plot_all(hill_coeff,i_pair) = ratio_plot;
                else
                    ratio_plot = squeeze(antag_ratio(i_species,:,:));
                    antag_ratio_plot_all(hill_coeff,i_pair) = ratio_plot;
                end



            end

        end

        if plot_score_heatmap

            figure(i_metric_name)

            paperpos = [0,0,500,300];
            papersize = [500,300];

            draw_pos=[10,10,480,280];

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

            %% setting colormap
            n = 13;
            xq = linspace(-0.3, 0.3, n);
            % define key colors at –0.3, 0, +0.3
            keyX  = [-0.3 , 0 , 0.3];
            keyRGB = [1  0  0    % blue
                1  1  1    % white
                0  0  1 ]; % red
            % interpolate each channel
            cmap = interp1(keyX, keyRGB, xq, 'linear');

            %% plot the heatmap
            if synergy_1_antag_0
                h = heatmap(synergy_ratio_plot_all(hill_coeff,:), ...
                    'Colormap', cmap, ...
                    'ColorLimits', [-0.3 0.3]);
            else
                h = heatmap(antag_ratio_plot_all(hill_coeff,:), ...
                    'Colormap', cmap, ...
                    'ColorLimits', [-0.3 0.3]);
            end

            colorbar off
            h.XDisplayLabels = repmat("", 1, 10);
            h.YDisplayLabels = repmat("", 1, 1);

            %% save score heatmap
            if synergy_1_antag_0
                saveas(gcf,strcat(fig_save_path,'synergy_ratio_hill',num2str(hill_coeff),'_',metric_name),'pdf')
            else
                saveas(gcf,strcat(fig_save_path,'antag_ratio_hill',num2str(hill_coeff),'_',metric_name),'pdf')
            end

            close()
        end
    end
end

if plot_bar_synergy_ratio_diff_hill

    figure

    fold_zoom = 0.8;
    paperpos = [0,0,200,100];
    papersize = [200,100];

    draw_pos=[10,10,180,80];

    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos*fold_zoom,'PaperSize', papersize*fold_zoom,'Position',draw_pos*fold_zoom)

    b2 = bar(synergy_ratio_plot_all');

    for k = 1:numel(b2)
        b2(k).BaseValue = -0.01;
    end
    ylim([-0.02,0.12])
    %% save score heatmap
    if synergy_1_antag_0
        saveas(gcf,strcat(fig_save_path,'synergy_ratio_bar_hill',num2str(hill_coeff),'_',metric_name,'_withlabel'),'pdf')
    else
        saveas(gcf,strcat(fig_save_path,'antag_ratio_bar_hill',num2str(hill_coeff),'_',metric_name,'_withlabel'),'pdf')
    end

    set(gca,'XTickLabel',{},'YTickLabel',{})
    %% save score heatmap
    if synergy_1_antag_0
        saveas(gcf,strcat(fig_save_path,'synergy_ratio_bar_hill',num2str(hill_coeff),'_',metric_name),'pdf')
    else
        saveas(gcf,strcat(fig_save_path,'antag_ratio_bar_hill',num2str(hill_coeff),'_',metric_name),'pdf')
    end

    close()

end



%% function plot heatmap from tensor
function fig_handle = plot_synergy_heatmap(tensor_A)
cell_num = 1000;
% wrong code, to delete
% A = tensor_A;
%
% synergy_mat = A(:,2:end,2:end) ...                % [M x N    x P]
%     - A(:,2:end,1)     ...                % [M x N    x 1]  → broadcasts over 3rd dim
%     - permute(A(:,1,2:end),[1,3,2]);      % [M x 1    x P]  → broadcasts over 2nd dim
%
% synergy_ratio = squeeze(sum(synergy_mat>0,1)/cell_num);

if 1 % after debugging
    A = tensor_A;
    cell_num = size(A,1);

    synergy_ratio = squeeze( ...
        sum( ...
        A(:,2:end,2:end) ...
        > ( A(:,2:end,1) + A(:,1,2:end) ), ...
        1 ) ...
        ) / cell_num;


end
% vis_mat_total(:,2:end,2:end) - vis_mat_total(:,2:end,1) * ones(1,5) - ones(5,1) * vis_mat_total(:,1,2:end);
%figure
heatmap(synergy_ratio)
fig_handle = gcf;

end


function fig_handle = plot_antagonism_heatmap(tensor_A)
cell_num = 1000;
A = tensor_A;

antagonism_mat = A(:,2:end,2:end) ...                % [M x N    x P]
    - max(A(:,2:end,1), permute(A(:,1,2:end),[1,3,2]));  % [M x 1    x P]  → broadcasts over 2nd dim

antagonism_ratio = squeeze(sum(antagonism_mat<0,1)/cell_num);
% vis_mat_total(:,2:end,2:end) - vis_mat_total(:,2:end,1) * ones(1,5) - ones(5,1) * vis_mat_total(:,1,2:end);
%figure
heatmap(antagonism_ratio)
fig_handle = gcf;

end

%% function for calculating synergy
function [tensor_sc_metric,vis_index_ligand1,vis_index_ligand2,vis_index_sti] = calculate_sc_metric(data_info,data,species_name,metric_name)

[tf, species_idx] = ismember(species_name, data_info.species_outputname);
if tf
else
    error('"%s" not in species_outputname!', species_name);
end

tensor_sc_metric = zeros(1000,1,1);

i_sti = 1;
for i_ligand1_dose = 1:2
    for i_ligand2_dose = 1:2
        vis_index_ligand1(i_ligand1_dose,i_ligand2_dose) = i_ligand1_dose;
        vis_index_ligand2(i_ligand1_dose,i_ligand2_dose) = i_ligand2_dose;
        vis_index_sti(i_ligand1_dose,i_ligand2_dose) = i_sti;
        switch metric_name
            case 'peak'
                tensor_sc_metric(:,i_ligand1_dose,i_ligand2_dose) = max(data.model_sim{i_sti}(species_idx:12:end,:),[],2);%mean(
            case 'total'
                tensor_sc_metric(:,i_ligand1_dose,i_ligand2_dose) = sum(data.model_sim{i_sti}(species_idx:12:end,:),2);
            otherwise
                error('no such metrics allowed')
        end
        i_sti = i_sti + 1;
    end
end


end

%% function for plotting synergy

function [heatmap_handle,fig_handle] = plot_synergy_heatmap_robust(tensor_A,min_factor)
cell_num = 1000;
A = tensor_A;

synergy_ratio = zeros(size(tensor_A,2)-1,size(tensor_A,3) -1 );

for j_tensor = 2:size(tensor_A,2) % different time range
    %j_tensor = 5;
    for k_tensor = 2:size(tensor_A,3) % different metrics
        %k_tensor = 6;
        x_vec = zeros(1000,1);
        y_vec = zeros(1000,1);
        z_vec = x_vec;
        base_vec = x_vec;
        x_vec(:) =  A(:,j_tensor,1)  ;
        y_vec(:) =  A(:,1,k_tensor)  ;
        z_vec(:) = A(:,j_tensor,k_tensor);
        base_vec(:) = A(:,1,1);

        synergy_ratio(j_tensor-1,k_tensor-1) = sum((z_vec > base_vec) & (z_vec > x_vec) & (z_vec > y_vec) & (z_vec  > ((x_vec + y_vec ) *1.000001))) /cell_num;%.000001
        % synergy_ratio(j_tensor-1,k_tensor-1) = sum((z_vec > min_factor) & (z_vec > x_vec) & (z_vec > y_vec) & (z_vec  > ((x_vec + y_vec ) *1.000001))) /cell_num;%.000001

        % if (j_tensor == 6) && (k_tensor == 3)
        %     idx = find((z_vec > base_vec) & (z_vec > x_vec) & (z_vec > y_vec) & (z_vec  > ((x_vec + y_vec ) *1.001)));
        %
        % end
    end
end

% vis_mat_total(:,2:end,2:end) - vis_mat_total(:,2:end,1) * ones(1,5) - ones(5,1) * vis_mat_total(:,1,2:end);
%figure
heatmap_handle = heatmap(synergy_ratio,'CellLabelColor','none');
fig_handle = gcf;

end

function [heatmap_handle,fig_handle] = plot_antagnism_heatmap_robust(tensor_A)
cell_num = 1000;
A = tensor_A;

antagnism_ratio = zeros(size(tensor_A,2)-1,size(tensor_A,3) -1 );

for j_tensor = 2:size(tensor_A,2)
    %j_tensor = 5;
    for k_tensor = 2:size(tensor_A,3)
        %k_tensor = 6;
        x_vec = zeros(1000,1);
        y_vec = zeros(1000,1);
        z_vec = x_vec;
        x_vec(:) =  A(:,j_tensor,1)  ;
        y_vec(:) =  A(:,1,k_tensor)  ;
        z_vec(:) = A(:,j_tensor,k_tensor);

        antagnism_ratio(j_tensor-1,k_tensor-1) = sum(z_vec < (max(x_vec, y_vec) -1e-4)) /cell_num;
    end
end
% vis_mat_total(:,2:end,2:end) - vis_mat_total(:,2:end,1) * ones(1,5) - ones(5,1) * vis_mat_total(:,1,2:end);
%figure
heatmap_handle = heatmap(antagnism_ratio,'CellLabelColor','none');
fig_handle = gcf;

end