function [input_paras] = initialize_input_paras_pair_ligand(ligand1,ligand2, index1, index2,cell_number)

index_default = [1,1+(6:2:15)];

if nargin<5
    cell_number = 500;
end
% % debug
% ligand1 = 'TNF';
% ligand2 = 'CpG';
%     index1 = index_default;
%     index2 = index_default;
% 
if (nargin <3) | isempty(index1)
    index1 = index_default;
end

if (nargin < 4) | isempty(index2)
    index2 = index_default;
end

if nargin >5
    error('too many inputs')
end

ligand_all = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};    % your cell array
index_ligand1 = find(strcmp(ligand_all,ligand1),1); 
index_ligand2 = find(strcmp(ligand_all,ligand2),1); 
sim_para_distribution_index = [2,3,4,5,6];

if isempty(index_ligand1)
    error(strcat(ligand1,'not found'));
elseif isempty(index_ligand2)
    error(strcat(ligand2,'not found'));
end

dose_length1 = length(index1);
dose_length2 = length(index2);
stim_number_total = dose_length2 * dose_length1;


%% initialize input_paras_all to include all ligands information, for later assign values for each pair.
dose_ratio = [0,10.^[-2.5:0.25:2.5]];
    dose_length = length(dose_ratio);

    input_paras_all.proj_num_vec = [mat2cell(2*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(3*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(4*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(5*ones(dose_length,1),ones(dose_length,1),[1]);
        mat2cell(6*ones(dose_length,1),ones(dose_length,1),[1])];
    input_paras_all.Num_sample = [1*ones(dose_length,1);
        1*ones(dose_length,1);
        1*ones(dose_length,1);
        1*ones(dose_length,1);
        1*ones(dose_length,1)]*600;%600;

    % TNF : index 6-15 (14)
    % Pam : index 6-15 (14)
    % CpG : index 6-15 (16)
    % LPS : index 6-15 (13)
    % pIC : index 6-15 (15)



    input_paras_all.proj_ligand_vec = {{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};
        {'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};{'TNF'};
        {'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};
        {'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};{'LPS'};
        {'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};
        {'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};{'CpG'};
        {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};
        {'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};{'PolyIC'}; {'PolyIC'}; {'PolyIC'};
        {'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};
        {'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'};};
    input_paras_all.proj_dose_str_vec = {{'d0'};{'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d0'};{'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d0'};{'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d0'};{'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'};
        {'d0'};{'d1'};{'d2'};{'d3'};{'d4'};{'d5'};{'d6'};{'d7'};{'d8'};{'d9'};{'d10'};{'d11'};{'d12'};
        {'d13'};{'d14'};{'d15'};{'d16'};{'d17'};{'d18'};{'d19'};{'d20'};{'d21'}};
    dose_val_a = [mat2cell((1*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((1*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((100*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((10000*dose_ratio)',ones(dose_length,1),[1]);
        mat2cell((100*dose_ratio)',ones(dose_length,1),[1])];
    % 
    % [mat2cell([0;(1*dose_ratio)'],ones(dose_length+1,1),[1]);
    %     mat2cell([0;(1*dose_ratio)'],ones(dose_length+1,1),[1]);
    %     mat2cell([0;(100*dose_ratio)'],ones(dose_length+1,1),[1]);
    %     mat2cell([0;(10000*dose_ratio)'],ones(dose_length+1,1),[1]);
    %     mat2cell([0;(100*dose_ratio)'],ones(dose_length+1,1),[1])];

    input_paras_all.proj_dose_val_vec = cell(size(dose_val_a));
    for i_proj = 1:length(input_paras_all.proj_dose_val_vec)
        input_paras_all.proj_dose_val_vec{i_proj} = dose_val_a(i_proj);
    end

    input_paras_all.var_fold_mat = 1;
    input_paras_all.var_fold_vec = 1;

    cal_codon =1;
    input_paras_filed_names = fieldnames(input_paras_all);


    %% assign values to input_paras

input_paras.proj_num_vec = mat2cell([index_ligand1*ones(stim_number_total,1),index_ligand2*ones(stim_number_total,1)],...
    ones(stim_number_total,1),[2]);
input_paras.Num_sample = cell_number *ones(stim_number_total,1);


ligand_pair = {ligand1,ligand2};
input_paras.proj_ligand_vec = repmat({ligand_pair}, stim_number_total, 1);
input_paras.var_fold_mat = 1;
input_paras.var_fold_vec = 1;

dose_str_ligand1 = input_paras_all.proj_dose_str_vec((index_ligand1-1) *  dose_length  + index1);
dose_str_ligand2 = input_paras_all.proj_dose_str_vec((index_ligand2-1) *  dose_length  + index2);
dose_val_ligand1 = dose_val_a((index_ligand1-1) *  dose_length  + index1); 
dose_val_ligand2 = dose_val_a((index_ligand2-1) *  dose_length  + index2);

input_paras.proj_dose_str_vec = cell(2,1);
input_paras.proj_dose_val_vec = cell(2,1);
i_sti =1;
for i_index_ligand1 = 1:dose_length1
    for i_index_ligand2 = 1: dose_length2

        input_paras.proj_dose_str_vec{i_sti} = [dose_str_ligand1{i_index_ligand1},dose_str_ligand2{i_index_ligand2}];
        input_paras.proj_dose_val_vec{i_sti} = {dose_val_ligand1{i_index_ligand1},dose_val_ligand2{i_index_ligand2}};
        i_sti = i_sti +1;
    end
end