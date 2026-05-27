% to move to src, if decided to use for publication
% data_save_file_path,data_file_list,fig_save_path
hill_coeff = 2;

%% initialize
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


save_fig_path = fig_save_path;


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


total_pair = 10;

%% draw 
if 0 % draw the pair and single-ligand
    %%
for i_pair = 1:total_pair

    %i_pair = 7;
    ligand1 = all_pairs{i_pair}{1};
    ligand2 = all_pairs{i_pair}{2};

    load(strcat(data_save_file_path,data_file_list{i_pair}))
    for i_species = 1
        %i_species = 5

        % calculate y limits for each speicies
        i_cond = 1;
        data_traj_all = [];
        for dose_ligand1 = 1:2
            for dose_ligand2 = 1:2
                data_traj_all = [data_traj_all;data.model_sim{i_cond}(i_species:12:end,:)];
                i_cond = i_cond+1;
            end
        end
        ylim_species = [0,0];
        ylim_species(1) = prctile(data_traj_all(:),0.25);
        ylim_species(2) = prctile(data_traj_all(:),99.75);
        %i_species
        %ylim_species


        % draw the heatmap
        i_cond = 1;
        paperpos=[0,0,100,130]*6;
        papersize=[100 130]*6;
        draw_pos=[10,10,90,120]*6;

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        for dose_ligand1 = 1:2
            for dose_ligand2 = 1:2
                species_name = data_info.species_outputname{i_species};
                data_traj = data.model_sim{i_cond}(i_species:12:end,:);

                subplot(2,2,dose_ligand1 + (dose_ligand2-1)*2)
                plot_heatmap_subplot(data_traj,ylim_species);

                i_cond = i_cond+1;
            end
        end
        save_name_file = strcat('heatmap_hill',num2str(hill_coeff),'_',num2str(i_pair),'_',ligand1,'_',ligand2,'_',species_name);

        saveas(gcf,strcat(save_fig_path,save_name_file),'pdf');
        close
    end

end

end

if 0 % only draw pairs
    %%
for i_pair = 1:total_pair

    %i_pair = 7;
    ligand1 = all_pairs{i_pair}{1};
    ligand2 = all_pairs{i_pair}{2};

    load(strcat(data_save_file_path,data_file_list{i_pair}))
    for i_species = 1
        %i_species = 5

        % calculate y limits for each speicies
        i_cond = 1;
        data_traj_all = [];
        for dose_ligand1 = 1:2
            for dose_ligand2 = 1:2
                data_traj_all = [data_traj_all;data.model_sim{i_cond}(i_species:12:end,:)];
                i_cond = i_cond+1;
            end
        end
        ylim_species = [0,0];
        ylim_species(1) = prctile(data_traj_all(:),0.25);
        ylim_species(2) = prctile(data_traj_all(:),99.75);
        %i_species
        %ylim_species


        % draw the heatmap
        
        paperpos=[0,0,100,130]*3;
        papersize=[100 130]*3;
        draw_pos=[10,10,90,120]*3;

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        i_cond = 4;
        species_name = data_info.species_outputname{i_species};
        data_traj = data.model_sim{i_cond}(i_species:12:end,:);

        plot_heatmap_subplot(data_traj,ylim_species);

        save_name_file = strcat('heatmap_hill',num2str(hill_coeff),'_L1L2_',num2str(i_pair),'_',ligand1,'_',ligand2,'_',species_name);

        saveas(gcf,strcat(save_fig_path,save_name_file),'pdf');
        close
    end

end

end

if 1 % only draw single-ligand
    %%
    draw_index = [1,1,2,3,4];
    draw_i_cond = [3,2,2,2,2];
for i_draw = 1:length(draw_index)

    i_pair = draw_index(i_draw);
    ligand1 = all_pairs{i_pair}{1};
    ligand2 = all_pairs{i_pair}{2};

    load(strcat(data_save_file_path,data_file_list{i_pair}))
    for i_species = 1
        %i_species = 5

        % calculate y limits for each speicies
        i_cond = 1;
        data_traj_all = [];
        for dose_ligand1 = 1:2
            for dose_ligand2 = 1:2
                data_traj_all = [data_traj_all;data.model_sim{i_cond}(i_species:12:end,:)];
                i_cond = i_cond+1;
            end
        end
        ylim_species = [0,0];
        ylim_species(1) = prctile(data_traj_all(:),0.25);
        ylim_species(2) = prctile(data_traj_all(:),99.75);
        %i_species
        %ylim_species


        % draw the heatmap
        
        paperpos=[0,0,100,130]*3;
        papersize=[100 130]*3;
        draw_pos=[10,10,90,120]*3;

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

        i_cond = draw_i_cond(i_draw);
        species_name = data_info.species_outputname{i_species};
        data_traj = data.model_sim{i_cond}(i_species:12:end,:);

        plot_heatmap_subplot(data_traj,ylim_species);

        save_name_file = strcat('heatmap_hill',num2str(hill_coeff),'_single_',num2str(i_pair),'_',all_pairs{i_pair}{4-draw_i_cond(i_draw)},'_',species_name);

        saveas(gcf,strcat(save_fig_path,save_name_file),'pdf');
        close
    end

end

end

%% subfuncitons for plot heatmaps
function fig_gcf = plot_heatmap_subplot(data_traj,ylim_vec)
% be default plot by integral
[~,data_order_index] = sort(sum(data_traj,2),"descend");


% figure()

if nargin ==1
    color_limits = [-0.001,0.3];
elseif nargin ==2
    color_limits = ylim_vec;
else
    error('too many/too few inputs')
end

cell_num=size(data_traj,1);

% subplot(1,length(vis_data_field),i_data_field)
h=heatmap(data_traj(data_order_index,:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limits);%[-0.001,0.2] for TNF

XLabels = 0:5:((size(data_traj,2)-1)*5);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/60);
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:cell_num;
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.YDisplayLabels = YCustomXLabels;

% xlabel('Time (hours)');
% ylabel(vis_data_field{i_data_field});
% clb=colorbar;
% clb.Label.String = 'NFkB(A.U.)';
colorbar('off')

set(gca,'fontsize',14,'fontname','Arial');
fig_gcf = gcf;

end