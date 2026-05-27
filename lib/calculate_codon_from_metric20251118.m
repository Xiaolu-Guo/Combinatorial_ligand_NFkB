function [collect_feature_vects] = calculate_codon_from_metric20251118(data_info,metric_structs) %, parameters)

% load the excel saving the metric lists that's used for
% each signal codon calculation
FeatureListFile = 'FeatureList_v3_Dur3.xlsx';
CodewordListFileSheet = 'Codewords';
FeatureListTable = readtable(FeatureListFile, 'Sheet', CodewordListFileSheet);
FeatureList = table2cell(FeatureListTable(:,1));
FeatureSpecifiers = table2cell(FeatureListTable(:,2));
CodewordList = table2cell(FeatureListTable(:,3));

% initialize the signaling codons
collect_feature_vects = struct;
collect_feature_vects.info_ligand = cell(2,1);
collect_feature_vects.info_dose_str = cell(2,1);
collect_feature_vects.info_data_type = cell(2,1);

for i_ids=1:length(metric_structs)
    collect_feature_vects.info_ligand{i_ids} = data_info.info_ligand{i_ids};
    collect_feature_vects.info_dose_str{i_ids} = data_info.info_dose_str{i_ids};
    collect_feature_vects.info_data_type{i_ids} = data_info.data_label{i_ids};
end

%% generate violin plots of codewords
nanzscore = @(x)(x-nanmean(x, 1))./nanstd(x, 0, 1);

codewords = unique(CodewordList);

for i = 1:length(codewords)
    vects = cell(1, length(metric_structs));
    for j= 1:length(FeatureList)
        if strcmp(CodewordList{j}, codewords{i})
            for k=1:length(metric_structs)
                metric_struct = metric_structs{k};
                feature = metric_struct.(FeatureList{j});
                if abs(FeatureSpecifiers{j}) ~= 0
                    feature = feature(:, abs(FeatureSpecifiers{j}));
                end
                if FeatureSpecifiers{j} < 0
                    feature = feature*-1;
                end
                % if parameters.calcResponders
                %     feature = feature(logical(metric_struct.responder_index));
                % end
                vects{k} = [vects{k} feature];
            end
        end
    end
    if size(vects{1}, 2) > 1  %%% for codewords that are composed of more than one metric-->take z score and then take average
        vects_zscore=vects(1:end);
        all_data_zscore = cell2mat(vects_zscore(:));
        mean_zscore= nanmean(all_data_zscore,1);
        std_zscore=nanstd(all_data_zscore);
        all_data = cell2mat(vects(:));
        all_data = (all_data-mean_zscore)./std_zscore;
        
        count = 1;
        for k = 1:length(metric_structs)
            vects{k} = nanmean(all_data(count:(size(vects{k}, 1)+count-1), :), 2);
            count = count + size(vects{k}, 1);
        end
    end
    all_data = cell2mat(vects(:)); %%% normalize each codeword
    all_data = nanzscore(all_data);
    % all_data = (all_data-prctile(all_data,0.5))/(prctile(all_data,99.5)-prctile(all_data,0.5)); %normalized except those extrema
    count = 1;
    for k = 1:length(metric_structs)
        vects{k} = all_data(count:(size(vects{k}, 1)+count-1), :);
        count = count + size(vects{k}, 1);
    end
    collect_feature_vects.(codewords{i}) = vects';
end

end
