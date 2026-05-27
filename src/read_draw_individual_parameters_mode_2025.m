%plot distribution

%

function [parameters,para_name] = read_draw_individual_parameters_mode_2025(proj_num,data_parent_path)


[proj_path,~] = subfunc_get_proj_path_2(proj_num,data_parent_path,'XG2025');

result_filepath= strcat(proj_path,'results/');

% estimates_cell=readtable(strcat(result_filepath,'summary.txt'));%,opts);
para_simul_cell=readtable(strcat(result_filepath,'IndividualParameters/simulatedIndividualParameters.txt'));
para_pred_cell=readtable(strcat(result_filepath,'IndividualParameters/estimatedIndividualParameters.txt'));

para_num = length(para_simul_cell.Properties.VariableNames)-2; % the first two columns are not parameters 

draw_para_num = para_num;%length(estimates.params);

for ii=1:draw_para_num
    
    para_name{ii} = erase(para_pred_cell.Properties.VariableNames{1+para_num*2+ii}, '_mode');% para_pred_cell.Properties.VariableNames{1+para_num*2+ii};
    %  histogram(para_pred_cell.(para_pred_cell.Properties.VariableNames{1+para_num*2+ii}),100,'Normalization','pdf')
    parameters(:,ii) = para_pred_cell.(para_pred_cell.Properties.VariableNames{1+para_num*2+ii});%mode
    
end
