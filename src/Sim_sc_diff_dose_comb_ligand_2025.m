function [] = Sim_sc_diff_dose_comb_ligand_2025(data_save_file_path,input_paras,cal_codon,save_metric_name)

% This code are for simulation of each individual cells responses to
% different doses of the ligand combination.
% input:
   % data_save_file_path: the file path where simulation data will be saved
   % input_paras
   % input_paras.ligand = {{L1,L2},{L1,L2},(repeating),{L1,L2}} gives the
   % ligand information
   % input_paras.dose = {{L1D1,L2D1},{L1D1,L2D2},...,{L1Dn,L2Dn}} gives the
   % doses infomration
   % cal_codon = 1, will caculate the codon and dynamic metrics
   % save_metric_name: the datafile name that will be saved
% output:
   % simulaiton and signaling codon files (if applicable) will be saved
   % under [data_save_file_path], named as [save_metric_name]

if nargin < 4
    cal_codon = 0;
end

% todo: check whether ligand is {{L1,L2},{L1,L2},(repeating),{L1,L2}}
% repeating/same for all conditions, otherwise report error

proj_num_vec = input_paras.proj_num_vec;
proj_ligand_vec = input_paras.proj_ligand_vec;
proj_dose_str_vec = input_paras.proj_dose_str_vec;
proj_dose_val_vec = input_paras.proj_dose_val_vec;

if isfield(input_paras,'Num_sample')
    Num_sample = input_paras.Num_sample;
else
    Num_sample = 500;
end

if length(Num_sample) == length(proj_num_vec)
    Num_sample = Num_sample;
else
    Num_sample = Num_sample * ones(size(proj_num_vec));
end

if isfield(input_paras,'var_fold_mat')
    if length(input_paras.var_fold_mat) == length(proj_num_vec)
        var_fold_mat = input_paras.var_fold_mat;
    else
        var_fold_mat = input_paras.var_fold_vec * ones(size(proj_num_vec));
    end
else
    var_fold_mat = ones(size(proj_num_vec));
end


gene_info.gene_type = {'wt'};
gene_info.gene_parameter_value_vec_genotype = cell(0);

%% parameter from the data
% todo: can be simplified to two lines, reading single-cell parameters
if 1
    data_save_file_path_1 = '../ref_data2025/';%_fay_parameter/';

    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));

    data_fitting = data;
    clear data
    for i_data = 1:length(data_fitting.pred_mode)

        data_fitting.pred_mode_cv{i_data} = std(data_fitting.pred_mode_filter_nan{i_data},[],2)./mean(data_fitting.pred_mode_filter_nan{i_data},2);
        [~,data_fitting.pred_mode_cv_order{i_data}] = sort(data_fitting.pred_mode_cv{i_data},'descend');

    end

    % input
    thresh_TNF=0.33;
    thresh_field = {'pred_mode_cv_order'};%pred_mode_cv_order osc
    data_fitting.parameters_mode_filter_TNF = cell(2,1);
    for i_data = 1:3
        index = data_fitting.(thresh_field{1}){i_data}(1: ceil(length(data_fitting.(thresh_field{1}){i_data}) * thresh_TNF));
        data_fitting.parameters_mode_filter_TNF{i_data} = data_fitting.parameters_mode_nan{i_data}(index ,:);
    end

    for i_data = 4:length(data_fitting.parameters_mode_nan)
        index = data_fitting.(thresh_field{1}){i_data}(1: ceil(length(data_fitting.(thresh_field{1}){i_data}) * thresh_TNF));
        data_fitting.parameters_mode_filter_TNF{i_data} = data_fitting.parameters_mode_nan{i_data}(: ,:);
    end
    % change the parameters for TNF only for high CV
    data_fitting.parameters_mode_nan = data_fitting.parameters_mode_filter_TNF;

    data_field_names = fieldnames(data_fitting);

    data_index = [1:6,10:12,14:19];

    for i_data_field = 1:length(data_field_names)
        data_new.(data_field_names{i_data_field}) = data_fitting.(data_field_names{i_data_field}) (data_index);
    end

    data_fitting_all = data_fitting;
    data_fitting = data_new;

end


%% sampling parameters

old_str = {'p','p','p','p','p','p','p','p','p','p'};
new_str = {'params','params','params','params','params','params','params','params','params','params'};
gene_info_all = cell(2,1);

for i_proj_num_vec = 1:length(proj_num_vec)

    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end

    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    var_fold = var_fold_mat(i_proj_num_vec);

    if i_proj_num_vec == 1
        switch length(proj_num_vec{i_proj_num_vec})
            case 1
                i_data = find(strcmp(data_fitting_all.info_ligand,sim_info.ligand));
                para_val_ele = [];
                for i_i_data= 1:length(i_data)
                    para_val_ele = [para_val_ele;data_fitting_all.parameters_mode_nan{i_data(i_i_data)} ] ; % [cellnumber x parameter]
                end
                i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1);

                rpt_time = ceil(Num_sample(i_proj_num_vec)*5/size(para_val_ele,1));
                para_val_total = [];
                for i_rpt = 1:rpt_time
                    para_val_total = [para_val_total;para_val_ele];
                end
                para_val = para_val_total(randperm(size(para_val_total,1),Num_sample(i_proj_num_vec)),:);
                est.name = cellfun(@replace,data_fitting.para_name{i_data},old_str(1:length(data_fitting.para_name{i_data})),new_str(1:length(data_fitting.para_name{i_data})),'UniformOutput',false);
                NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
                shift_index = find(strcmp(est.name,'shift'));
                parameter_index = setdiff(setdiff(1:length(est.name),NFkB_index,'stable'),shift_index,'stable');

                % non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
                gene_info.parameter_name_vec = {est.name(parameter_index)};
                gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,parameter_index)';
                gene_info.species_name_vec= {{'NFkB'}};
                gene_info.species_value_vec_genotype{1}{1} = para_val(:,NFkB_index)';

            otherwise
                % sample bi-stimulation

                [para_sample_multi_ligand,estimates_2ligand] = sample_fitting_multi_receptor_2023_05(proj_num_vec{i_proj_num_vec},Num_sample(i_proj_num_vec),data_fitting,data_fitting_all,sim_info,'alldose');

                gene_info.parameter_name_vec = {estimates_2ligand.name(estimates_2ligand.non_NFkB_index)};
                gene_info.gene_parameter_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.non_NFkB_index)';
                gene_info.species_name_vec= {{'NFkB'}};
                gene_info.species_value_vec_genotype{1}{1} = para_sample_multi_ligand(:,estimates_2ligand.NFkB_index)';
        end
    else
    end

    % stimili info

    save_filename = strcat('202305_para_sampled_fitting_NFkBinit_',replace(num2str(var_fold),'.','p'),'varfoldchange');

    for i_ligand = 1:length(sim_info.ligand)
        save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand},...
            '_',replace(replace(sim_info.dose_str{i_ligand},'/',''),'.','p'));
    end

    % species that will be saved
    % must be r x 1, for each cell i must be ri x 1
    data_info.species_outputname = {'nucNFkB';'IKK';'TAK1';'TRAF6';'TRIF';'MyD88';...
        'C1';'TLR4LPSen';'TLR4LPS';'TLR9_CpG';'TLR3_polyIC';'TLR2_P3CSK'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'IKK'};{'TAK1'};{'TRAF6'};{'TRIF'};{'MyD88'};...
        {'C1'};{'TLR4LPSen'};{'TLR4LPS'};{'TLR9_CpG'};{'TLR3_polyIC'};{'TLR2_P3CSK'}};
    %
    % data_info.species_outputname = {'nucNFkB';'MyD88';'TLR4';'TLR2';'TLR3';'TLR9';'IKK';'TAK1';'IkBamRNA'};
    % data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'MyD88'};{'TLR4'};{'TLR2'};{'TLR3'};{'TLR9'};{'IKK'};{'TAK1'};{'IkBat'}};
    % data_info.save_file_name = save_filename; % only beginning, no .mat

    if isfield(input_paras,'IKK_act_hill')

        gene_info.parameter_name_vec{1}(end+1) = {'params67n2'};
        gene_info.gene_parameter_value_vec_genotype{1}{1}(end+1,:)   = ...
            input_paras.IKK_act_hill * ones(1,size(gene_info.gene_parameter_value_vec_genotype{1}{1},2));

    end

    gene_info_all{i_proj_num_vec} = gene_info;
    
    % History version: before 2025/09/15
    % sim_data_tbl = genotype_sim_save_2025_CpGpolyIC_compete(sim_info,data_info,gene_info);

    % version 2025/09/15: update CpG pIC SR induced endosomal
    sim_data_tbl = genotype_sim_save_2025_CpGpolyIC_SR(sim_info,data_info,gene_info);

    i_ligand = 1;
    ligand_str= proj_ligand_vec{i_proj_num_vec}{i_ligand};
    dose_str = proj_dose_str_vec{i_proj_num_vec}{i_ligand};

    for i_ligand = 2:length(proj_ligand_vec{i_proj_num_vec})
        ligand_str = strcat(ligand_str,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand});
        dose_str = strcat(dose_str,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand});
    end


    data.info_ligand{i_proj_num_vec} = ligand_str;
    data.info_dose_index{i_proj_num_vec} = 1;
    data.info_dose_str{i_proj_num_vec} = dose_str;
    data.info_num_cells{i_proj_num_vec} = size(sim_data_tbl.trajectory,1);
    [~, data.order{i_proj_num_vec}] = sort(max(sim_data_tbl.trajectory,[],2),'descend');


    if cal_codon
        data.model_sim{i_proj_num_vec} = sim_data_tbl.trajectory(:,1:5:end);
    else
        data.model_sim{i_proj_num_vec} = sim_data_tbl.trajectory(:,1:1:end);
    end
    % save(save_metric_name,'data');
end

    data.info_dose_val = proj_dose_val_vec;

if cal_codon
    data_cal_codon = data;

    for i_cond = 1:length(data.model_sim)
        data_cal_codon.model_sim{i_cond} = data.model_sim{i_cond}(1:length(data_info.species_outputname):size(data.model_sim{i_cond},1),:);

    end

    data_cal_codon.exp = data_cal_codon.model_sim;

    % cd(codon_info.codon_path)
    vis_data_field = {'model_sim'};%,'sample'};
    data_label = {'simulation'};%,'sample'};
    [collect_feature_vects,metrics] = calculate_codon_2023(data_cal_codon,vis_data_field,data_label);%,  parameter

    % for i_metrics = length(metrics)
    % pk_time(i_metrics) = mean( metrics{i_metrics}.pk1_time );
    % osc_ratio(i_metrics) = sum(metrics{i_metrics}.oscpower>2e-5)/length(metrics{i_metrics}.oscpower);
    % end

    save(strcat(data_save_file_path,save_metric_name),'data','metrics','collect_feature_vects','data_info','gene_info_all');

else
    save(strcat(data_save_file_path,save_metric_name),'data','data_info','gene_info_all');

end


