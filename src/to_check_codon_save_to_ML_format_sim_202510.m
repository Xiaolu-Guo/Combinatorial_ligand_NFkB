function [] = codon_save_to_ML_format_sim_202510(data_save_file_path,data_file,codon_save_path,index_to_save)%,save_name_tag
% this is for experimental data, if simulation data, then change line 54:
% metric_cal{i_id}.(field_list{i_field}) = metrics{i_id}.(field_list{i_field})(1:9:end,:);

if nargin < 3
    codon_save_path = data_save_file_path;
end

if nargin < 4
    index_to_save = [4,6,22,15,11,8,5,9,3,25,19,26,18,14,20,21,28];%31
    index_to_save = {[4];% 1-TNF
            [6]; %2- lps
            [22]; %3 cpg
            [15]; %4
            [11]; %5
            [8];%6
            [5];%7
            [9];%8
            [3];%9
            [25];%10
            [19];%11
            [26];%12
            [18];%13
            [];%14
            [14];%15
            [];%16
            [];%17
            [];%18
            [20];%19
            [];%20
            [21];%21
            [];%22
            [];%23
            [];%24
            [];%25
            [];%26
            [];%27
            [];%28
            [];%29
            [];%30
            [28]};%31

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

        metric_cal{i_id}.(field_list{i_field}) = metrics{i_id}.(field_list{i_field})(:,:);
    end
end

collect_feature_vects = calculate_codon_from_metric2025(data_info,metric_cal);

nfkb_codon_all = [];
nfkb_id_all = [];
index_all = 1:length(collect_feature_vects.(codon_list{1}));
for i_index = 1:length(index_to_save)% in total length(index_data) conditions
    % i_ligand = index_to_save{i_index};
    if ~isempty(index_to_save{i_index})
        i_ligand = index_to_save{i_index};
    nfkb_codon = [];
    nfkb_id = [];
    for i_codon =1:length(codon_list)
        nfkb_codon = [nfkb_codon,collect_feature_vects.(codon_list{i_codon}){i_ligand}];
    end
    nfkb_id = (i_index-1)*ones(size(nfkb_codon,1),1);
    nfkb_codon_all = [nfkb_codon_all;nfkb_codon];
    nfkb_id_all = [nfkb_id_all;nfkb_id];
    %[~, order_data] = sort(max(metrics{index_data(i_ligand)}.time_series,[],2),'descend')
    % h=heatmap(metrics{index_data(i_ligand)}.time_series(order_data,:) ,'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
    end
end

writematrix(nfkb_codon_all,strcat(codon_save_path,'Comb_ligands_exp_conditions_X_codon_','.csv'));
writematrix(nfkb_id_all,strcat(codon_save_path,'Comb_ligands_exp_conditions_y_codon_','.csv'));

end

