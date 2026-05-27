addpath('./lib/')
addpath('./src/')
addpath('./bin/')

%% [run sim single-cell matched]

if 1

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


    for i_pair = 1:10

        ligand1 = all_pairs{i_pair}{1};
        ligand2 = all_pairs{i_pair}{2};
        cell_num = 1000;

        input_paras = initialize_input_paras_pair_ligand(ligand1,ligand2,[],[],cell_num/2);
        index_dose2 = [1;3;13;15];

        input_paras.Num_sample = input_paras.Num_sample(index_dose2);
        input_paras.proj_dose_str_vec = input_paras.proj_dose_str_vec(index_dose2);
        input_paras.proj_dose_val_vec = input_paras.proj_dose_val_vec(index_dose2);
        input_paras.proj_ligand_vec = input_paras.proj_ligand_vec(index_dose2);
        input_paras.proj_num_vec = input_paras.proj_num_vec(index_dose2);

        cal_codon =1;

        data_save_file_path = '../test_data/';

        tic
        save_metric_name = strcat('Sim1008_',ligand1,'_',ligand2,'_sc_dose2_IKK_hill2.mat');
        Sim_sc_diff_dose_comb_ligand_2025_CpGpIC_SR(data_save_file_path,input_paras,cal_codon,save_metric_name)
        toc

    end
end