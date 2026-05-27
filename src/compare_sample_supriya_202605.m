% clear all
% load the codon for all sampled data and supriya and ade's data for
% comparison
% comapring ade's data and supriya's data: the batch effects.
% comparing surpiya's data and model prediction
% check each chunk to set 1 or 0 to run or not run the corresponding part
%
cal_codon_all = 0;
cal_codon_dual_ligand = 1;
cal_codon = 1;
compete_or_not = 1;

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

if compete_or_not
    save_metric_name_sim = 'Sim2025_all_comb_SR.mat'; % for CpG-pIC compete, with pIC-TLR deg increased
else
    save_metric_name_sim = 'Sim2025_all_comb_CpG_PIC_no_compete.mat';
end

data_save_file_path_exp = data_save_file_path;
data_save_file_path_sim = data_save_file_path;


if cal_codon_all % not yet debug: 11/18/2025; maybe needed for Figure S2 update
    %%  calculate/load the codon for all data
    if cal_codon
                
        % Supriya's data
        load(strcat(data_save_file_path_exp,'Exp_data_metrics_2024.mat'))
        data_supriya = data;
        metrics_supriya = metrics;
        
        % TNF 10ng/mL       supriya: 03 index: 4;       ade index: 3    sample: 13
        % Pam3CSK 100ng/mL  supriya: 24 index: 27;      ade index: 18   sample: 12
        % CpG 100nM         supriya: 23 index: 22;      ade index: 11   sample: 11
        % LPS 10ng/mL       supriya: 04 index: 6;       ade index: 6    sample: 14
        % PolyIC 100ug/mL   supriya: 20 index: 15;      ade index: 16   sample: 15
        
        % goal: supriya's data, ade's data, fitting data, and sampled data for
        % benchamarking
        % compare the prediction of smapled data with supriya's data
        
        % TNF 10ng/mL Pam3CSK 100ng/mL supriya: 03 index: 3
        % TNF 10ng/mL CpG 100nM supriya: 03 index: 5
        % TNF 10ng/mL LPS 10ng/mL supriya: 04 index: 8
        % TNF 10ng/mL PolyIC 100ug/mL supriya: 04 index: 9
        % Pam3CSK 100ng/mL CpG 100nM supriya: does not apply
        % Pam3CSK 100ng/mL LPS 10ng/mL supriya: 23 index: 26
        % Pam3CSK 100ng/mL PolyIC 100ug/mL supriya: 13 index: 14
        % CpG 100nM LPS 10ng/mL supriya: 23 index: 25
        % CpG 100nM PolyIC 100ug/mL supriya: 21 index: 18
        % LPS 10ng/mL PolyIC 100ug/mL supriya: 21 index: 19
        
        % ade_index = [3;18;11;6;16];
        supriya_index = [4;27;22;6;15];
        supriya_dual_ligand_index = [3;5;8;9;26;14;25;18;19];
        
        % sampled data
        data_filename = save_metric_name_sim;%'Sim5_codon_all5dose_metric.mat';
        load(strcat(data_save_file_path_exp,data_filename))
        sample_dual_ligand_index = [6:15];
        sample_index = [1:5];
        
        % check whether 1:9:end or else
        metric_fields_name = fieldnames(metrics{1});
        for i_metric = 1:length(metrics)
            
            for i_metric_fields = 1:length(metric_fields_name)
                metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
            end
        end
        
        data_sample = data;
        % metrics_sample = metrics;
        
        metrics = [metrics_supriya(supriya_index),... Supriya's single liagand data
            metrics_sample(sample_index),... Sampling single liagand data
            metrics_supriya(supriya_dual_ligand_index),... Supriya's dual ligand data
            metrics_sample(sample_dual_ligand_index)]; % sampling dual ligand data
        
        
        data_info.info_ligand = [data_supriya.info_ligand(supriya_index),... Supriya's single liagand data
            data_sample.info_ligand(sample_index),... Sampling single liagand data
            data_supriya.info_ligand(supriya_dual_ligand_index),... Supriya's dual ligand data
            data_sample.info_ligand(sample_dual_ligand_index)]; % sampling dual ligand data
        
        
        for i_info_ligand = 1:length(data_info.info_ligand)
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'polyIC','PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_LPS','LPS_PolyIC');
        end
        
        
        data_info.info_dose_str = [data_supriya.info_dose_str(supriya_index),... Supriya's single liagand data
            data_sample.info_dose_str(sample_index),... Sampling single liagand data
            data_supriya.info_dose_str(supriya_dual_ligand_index),... Supriya's dual ligand data
            data_sample.info_dose_str(sample_dual_ligand_index)]; % sampling dual ligand data
        
        data_info.data_label(1:length(supriya_index)) = {'Exp_Supriya'};
        data_info.data_label(end+1:end+length(sample_index)) = {'Sampling'};
        data_info.data_label(end+1:end+length(supriya_dual_ligand_index)) = {'Exp_Supriya_dual_ligand'};
        data_info.data_label(end+1:end+length(sample_dual_ligand_index)) = {'Sampling_dual_ligand'};
        
        % [collect_feature_vects,metrics_new] = calculate_codon_from_metric2023(data_info,metrics); %,  parameter
        
        [collect_feature_vects,metrics_new] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics); %,  parameter
        metrics = metrics_new;
        %
        collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
        collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
        collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
        %
        % save(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'),'collect_feature_vects','metrics','data_info')
        save(strcat(data_save_file_path_sim,'All_codon_dual_ligand_202511_nonminmaxscaled.mat'),'collect_feature_vects','metrics','data_info')
        
        % goal: visualize supriya, ade, fitting, sampling (2 seprate plots) as
        % benchmarking; then visualize supriya dual ligand vs sampling.
    else
        % load(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'))
        % load(strcat(data_save_file_path_1,'All_codon_dual_ligand_202306.mat'))
        load(strcat(data_save_file_path_sim,'All_codon_dual_ligand_202511_nonminmaxscaled.mat'))
        
    end
end

if cal_codon_dual_ligand % successfully run: 11/18/2025
    %%  calculate/load the codon for dual ligand
    if cal_codon
        
        % Supriya's data
        load(strcat(data_save_file_path_exp,'Exp_data_metrics_2024.mat'))
        data_supriya = data;
        metrics_supriya = metrics;
        
        % TNF 10ng/mL       supriya: 03 index: 4;       ade index: 3    sample: 13
        % Pam3CSK 100ng/mL  supriya: 24 index: 27;      ade index: 18   sample: 12
        % CpG 100nM         supriya: 23 index: 22;      ade index: 11   sample: 11
        % LPS 10ng/mL       supriya: 04 index: 6;       ade index: 6    sample: 14
        % PolyIC 100ug/mL   supriya: 20 index: 15;      ade index: 16   sample: 15
        
        % goal: supriya's data, ade's data, fitting data, and sampled data for
        % benchamarking
        % compare the prediction of smapled data with supriya's data
        
        % TNF 10ng/mL Pam3CSK 100ng/mL supriya: 03 index: 3
        % TNF 10ng/mL CpG 100nM supriya: 03 index: 5
        % TNF 10ng/mL LPS 10ng/mL supriya: 04 index: 8
        % TNF 10ng/mL PolyIC 100ug/mL supriya: 04 index: 9
        % Pam3CSK 100ng/mL CpG 100nM supriya: does not apply
        % Pam3CSK 100ng/mL LPS 10ng/mL supriya: 23 index: 26
        % Pam3CSK 100ng/mL PolyIC 100ug/mL supriya: 13 index: 14
        % CpG 100nM LPS 10ng/mL supriya: 23 index: 25
        % CpG 100nM PolyIC 100ug/mL supriya: 21 index: 18
        % LPS 10ng/mL PolyIC 100ug/mL supriya: 21 index: 19

        supriya_index = [4;27;22;6;15];
        supriya_dual_ligand_index = [3;5;8;9;26;14;25;18;19];
        
        % sampled data
        
        data_filename = save_metric_name_sim; %'Sim5_codon_all5dose_metric.mat';
        load(strcat(data_save_file_path_exp,data_filename))
        sample_dual_ligand_index = [6:15];
        sample_index = [1:5];
        
        metric_fields_name = fieldnames(metrics{1});
        for i_metric = 1:length(metrics)
            
            for i_metric_fields = 1:length(metric_fields_name)
                metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
            end
        end
        
        data_sample = data;
        % metrics_sample = metrics;
        
        metrics = [metrics_supriya(supriya_index),... Supriya's single liagand data
            metrics_sample(sample_index),... Sampling single liagand data
            metrics_supriya(supriya_dual_ligand_index),... Supriya's dual ligand data
            metrics_sample(sample_dual_ligand_index)]; % sampling dual ligand data
        
 
        data_info.info_ligand = [data_supriya.info_ligand(supriya_index),... Supriya's single liagand data
            data_sample.info_ligand(sample_index),... Sampling single liagand data
            data_supriya.info_ligand(supriya_dual_ligand_index),... Supriya's dual ligand data
            data_sample.info_ligand(sample_dual_ligand_index)]; % sampling dual ligand data
        
        
        for i_info_ligand = 1:length(data_info.info_ligand)
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'polyIC','PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_LPS','LPS_PolyIC');
        end
        
        
        data_info.info_dose_str = [data_supriya.info_dose_str(supriya_index),... Supriya's single liagand data
            data_sample.info_dose_str(sample_index),... Sampling single liagand data
            data_supriya.info_dose_str(supriya_dual_ligand_index),... Supriya's dual ligand data
            data_sample.info_dose_str(sample_dual_ligand_index)]; % sampling dual ligand data
        
        data_info.data_label = {''};
        data_info.data_label(1:length(supriya_index)) = {'Exp_Supriya'};
        data_info.data_label(end+1:end+length(sample_index)) = {'Sampling'};
        data_info.data_label(end+1:end+length(supriya_dual_ligand_index)) = {'Exp_Supriya_dual_ligand'};
        data_info.data_label(end+1:end+length(sample_dual_ligand_index)) = {'Sampling_dual_ligand'};        

        [collect_feature_vects,metrics_new] = calculate_codon_from_metric202511_nonminmaxscaled(data_info,metrics); %,  parameter
        metrics = metrics_new;
        %
        collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
        collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
        collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
        %
        % save(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'),'collect_feature_vects','metrics','data_info')
        save(strcat(data_save_file_path_sim,'All_exp_sim_codon_dual_ligand_202511_nonminmaxscaled.mat'),'collect_feature_vects','metrics','data_info')
        
        % goal: visualize supriya, ade, fitting, sampling (2 seprate plots) as
        % benchmarking; then visualize supriya dual ligand vs sampling.
    else
        % load(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'))
        % load(strcat(data_save_file_path_1,'All_codon_dual_ligand_202306.mat'))
        load(strcat(data_save_file_path_sim,'All_exp_sim_codon_dual_ligand_202511_nonminmaxscaled.mat'))
        
    end
end

%% different data set index setting
index_supriya = 1:5;
% index_ade = 6:10;
% index_fitting = 11:15;
index_sampling = 6:10;
index_supriya_dual = 11:19;
index_sampling_dual = [20:23,25:29];


%%
if 1 % Figure 2
    %% violin plot of codon (duration and total activity) distribution under co stim vs single stim for comparison, only picked three
% 

% A + B : [supriya A, supriya B, supriya A + B],  [sampling A, sampling B, sampling A + B]
% TNF + Pam3CSK : [1,2,11],  [6,7,20]
% TNF + CpG : [1,3,12],  [6,8,21]
% TNF + LPS : [1,4,13],  [6,9,22]
% TNF + PolyIC : [1,5,14],  [6,10,23]
% Pam3CSK + CpG : x,  [7,8,24]
% Pam3CSK + LPS : [2,4,15],  [7,9,25]
% Pam3CSK + PolyIC : [2,5,16],  [7,10,26]
% CpG + LPS : [3,4,17],  [8,9,27]
% CpG + PolyIC : [3,5,18],  [8,10,28]
% LPS + PolyIC : [4,5,19],  [9,10,29]

    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    clear index_violin collect_feature_vects_draw
    codon_fields = fieldnames(collect_feature_vects);
    index_violin_exp = {[1,4,13]%pick TNF LPS and TNF-LPS
        [3,5,18]%pick CpG polyIC and CpG-pIC
        [4,5,19]};%pick LPS polyIC and LPS-pIC
    
    index_violin_sampling = {[6,7,20]%pick TNF LPS and TNF-LPS
        [8,9,27]%pick CpG polyIC and CpG-pIC
        [7,9,25]};%pick LPS polyIC and LPS-pIC
    
    index_violin_vec= {index_violin_exp,index_violin_sampling};
    color_vec = {[0,0,0],[0,1,0]};
    color_box_vec = {[0.3,0.3,0.3],[0,0.5,0]};
    data_type = {'exp','sampling'};
    for i_index = 1:length(index_violin_vec)
        index_violin = index_violin_vec{i_index};
        for i_index_violin = 1:length(index_violin)
            clear collect_feature_vects_draw
            
            for i_codon_fields = 1:length(codon_fields)
                collect_feature_vects_draw.(codon_fields{i_codon_fields}) = collect_feature_vects.(codon_fields{i_codon_fields})(index_violin{i_index_violin});
            end
            
            for i_codon = 1:length(codon_list)
                collect_feature_vects_draw.(codon_list{i_codon}) = combined_cell_zscore(collect_feature_vects_draw.(codon_list{i_codon}));
            end
            
            
            for i_data_types = 1:length(collect_feature_vects_draw.info_data_type)
                collect_feature_vects_draw.info_data_type{i_data_types} = data_type{i_index};
            end
            
            fig_opt.save_file = strcat(fig_save_path,'Co_Stim_picked3_',data_type{i_index},'_compare_codon_distrib_202308');
            %     fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
            %     fig_opt.paper_opt.papersize=[220 180]*3;
            fig_opt.paper_opt.paperpos=[0,0,100,50];
            fig_opt.paper_opt.papersize=[100 50];
            fig_opt.distri_color = color_vec(i_index);
            fig_opt.codon = {'Duration','TotalActivity'};
                     
            violin_plot_codon_compare_co_sti_202308(collect_feature_vects_draw,fig_opt)

            if 1
                codon_list = {'Duration','TotalActivity'};
                
                for i_codon=1:length(codon_list)
                    
                    vects = collect_feature_vects_draw.(codon_list{i_codon});
                    median_codon = cellfun(@median,vects(1:3));
                    codon_list{i_codon}
                    collect_feature_vects_draw.info_ligand{3}
                    collect_feature_vects_draw.info_data_type{3}
                    diff = median_codon(3) - max(median_codon(1:2))
                    if 1 %  One-sided Mann–Whitney U tests
                        [p1l, ~, ~] = ranksum(vects{1}, vects{3}, 'tail','left')
                        [p1r, ~, ~] = ranksum(vects{1}, vects{3}, 'tail','right')

                        [p2l, ~, ~] = ranksum(vects{2}, vects{3}, 'tail','left')
                        [p2r, ~, ~] = ranksum(vects{2}, vects{3}, 'tail','right')
                    end
                end
            end
            
        end
    end
    
    % [0,1,0],
    
end

%% funciton zscore

function zscores = combined_cell_zscore(dataCells)
% Concatenate all data from the cells to compute the overall mean and standard deviation
allData = vertcat(dataCells{:});
mu = nanmean(allData);
sigma = nanstd(allData);

% Initialize the output cell array with the same size as the input
zscores = cell(size(dataCells));

% Loop through each cell and compute the z-scores for the data inside using the overall mean and std
for i = 1:numel(dataCells)
    zscores{i} = (dataCells{i} - mu) / sigma;
end
end


