%% initializing
data_save_file_path = '../test_data/';
fig_save_path = '../subfigures/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path);
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

if 1
    % run all conditions
    Initialize_all_combinatorial_ligands

    save_metric_name = 'Sim2025_all_comb_SR.mat';
    cal_codon =1;
    
    % the differences between this version and before: para_sample_fitting_sim_sti_doses_var_2023_05
    % pIC-TLR3 degradation increases x 2
    % pIC-CpG competing for endosmal transport resources, using Machaelis
    % Menten kinetics
    % pIC MW updated: 1000kDa(+) from 200kDa(+)
    para_sample_fitting_sim_sti_doses_var_2025_CpGpIC_SR(data_save_file_path,input_paras,cal_codon,save_metric_name)

end
