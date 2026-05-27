function [] = data_save_monolix_format(data_save_file_path,data_file)%,save_name_tag

if nargin < 3
    codon_save_path = data_save_file_path;
end
% if nargin < 4
%     index_to_save = [4,6,22,15,11,8,5,9,3,25,19,26,18,14,20,21,28];%31
%     index_to_save = {[4];% 1-TNF
%             [6]; %2- lps
%             [22]; %3 cpg
%             [15]; %4
%             [11]; %5
%             [8];%6
%             [5];%7
%             [9];%8
%             [3];%9
%             [25];%10
%             [19];%11
%             [26];%12
%             [18];%13
%             [];%14
%             [14];%15
%             [];%16
%             [];%17
%             [];%18
%             [20];%19
%             [];%20
%             [21];%21
%             [];%22
%             [];%23
%             [];%24
%             [];%25
%             [];%26
%             [];%27
%             [];%28
%             [];%29
%             [];%30
%             [28]};%31
%
% end
load(strcat(data_save_file_path,data_file))

% index_to_save_monolix = [6,24,...LPS
%     1,10,22,...CpG
%     25];%CpG+LPS

index_to_save_monolix = [13,16,18];%PolyIC + CpG

%% the path to save the data
data_save_path = '../revision/data_SAEM_format/';

%% rescaling the data and save

if 1
    % the macrophage cell volume is 6 pl, measured by Stefanie Luecke
    NFkB_max_range = [0.25, 1.75]/6;

    % We assume cells in response to LPS 100ng can reach
    % the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
    % (Pam3CSK might not be able to reach the highest Nuc NFkB conc)

    data_for_rescale = data.exp{6};
    data_for_rescale = [data_for_rescale;data.exp{24}]; % LPS 10 ng
    rescale_factor = NFkB_max_range(2) / prctile(max(data_for_rescale,[],2),85) ;%for 100ng LPS, 90 percentile, for 10ng LPS, 99 percentile


    % sti_info = 'LPS';
    % sti_AMT = 10;
    % sti_ADM = 3;
    % data_to_rescale = [data.exp{6};data.exp{24}];

    % sti_info = 'CpG';
    % sti_AMT = 100;
    % sti_ADM = 2;
    % data_to_rescale = [data.exp{1};data.exp{10};data.exp{22}];

    % sti_info = 'LPS-CpG';
    % sti_AMT = [10,100];
    % sti_ADM = [3,2];
    % data_to_rescale = [data.exp{25}];

    sti_info = 'pIC_CpG';
    sti_AMT = [100000,100];
    sti_ADM = [5,2];
    data_to_rescale = [data.exp{13};data.exp{16};data.exp{18}];

    data_to_monolix_format = data_to_rescale * rescale_factor;

end

%% SAEM data format
write_SAEM_data_filename = strcat(data_save_path,'SAEM_',sti_info,'.txt');
TimePtsWeight = 1;
time_each_frame = 5;%min

if 1
    if isfile(write_SAEM_data_filename)
        delete(write_SAEM_data_filename);
    end

    file_to_write=fopen(write_SAEM_data_filename,'w'); %create the txt data file

    fprintf(file_to_write,'%s \t %s \t %12s \t %12s \t %12s \t %12s \r\n','ID','TIME','AMT','ADM','Y','TimePtWeight'); %write the txt data file

    cell_id=1;

    time_stimuli=0.15;
    time_stimuli_2 = 0.151;

    cell_num = size(data_to_monolix_format,1);
    stim_info.ADM = sti_ADM;
    stim_info.dose_val = sti_AMT;

    exp_data= data_to_monolix_format;

    time_matrix=ones(cell_num,1) * (-5:time_each_frame: 480);

    exp_data_wr=exp_data';

    time_matrix_wr=time_matrix';

    for i_cell=1:cell_num

        fprintf(file_to_write,'%d \t %d \t %12s \t %12s \t %12.5g \t %12.5g \r\n',cell_id,time_matrix_wr(1,i_cell),'.','.',exp_data_wr(1,i_cell),1); %write the txt data file

        fprintf(file_to_write,'%d \t %d \t %12s \t %12s \t %12.5g \t %12.5g \r\n',cell_id,time_matrix_wr(2,i_cell),'.','.',exp_data_wr(2,i_cell),1); %write the txt data file

        if isscalar(stim_info.dose_val)
            fprintf(file_to_write,'%d \t %.2f \t %12.6g \t %12d \t %12s \t %12s \r\n',cell_id,time_stimuli,stim_info.dose_val,stim_info.ADM,'.','.'); %write the txt data file
        elseif length(stim_info.dose_val) == 2
            fprintf(file_to_write,'%d \t %.2f \t %12.6g \t %12d \t %12s \t %12s \r\n',cell_id,time_stimuli,stim_info.dose_val(1),stim_info.ADM(1),'.','.'); %write the txt data file

            fprintf(file_to_write,'%d \t %.3f \t %12.6g \t %12d \t %12s \t %12s \r\n',cell_id,time_stimuli_2,stim_info.dose_val(2),stim_info.ADM(2),'.','.'); %write the txt data file

        else
            error('wrong number of stim')
        end

        for k=3: size(data_to_monolix_format,2)

            fprintf(file_to_write,'%d \t %d \t %12s \t %12s \t %12.5g \t %12.5g \r\n',cell_id,time_matrix_wr(k,i_cell),'.','.',exp_data_wr(k,i_cell),TimePtsWeight);
        end

        cell_id=cell_id+1;
    end


    fclose(file_to_write);

end

