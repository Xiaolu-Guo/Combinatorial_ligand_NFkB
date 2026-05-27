
% relpace/ delete / or move to debug: the original codes single_cell_plot_script.m
if 1
    %% synergy condition: under 5 by 5 condtion
    switch metric_name
        case'peak'
            i_row_plot_cond = 2;
            i_col_plot_cond = 1;
        case 'total'
            i_row_plot_cond = 2;
            i_col_plot_cond = 2;
    end
    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
    responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

    % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
    synergy_score_nfkb_cond = Synergy_score(:,i_row_plot_cond,i_col_plot_cond);

    i_cond = i_row_plot_cond * 6 + i_col_plot_cond + 1;

    params_cond = gene_info_all{i_cond}.gene_parameter_value_vec_genotype{1}{1}';
    params_names = gene_info_all{i_cond}.parameter_name_vec{1};

    % remove extra useless parameters: first p52n2, overwrite by second p52n2
    params_cond = params_cond(:,[1:2,4:end]);
    params_names = params_names([1:2,4:end]);

    params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);
    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));

    % generate random variable for control
    rng(19);
    para_rand = rand(1000,1);

    % PAWN analysis ,interval_bounds, ft
    [KS] = PAWN_given_data(synergy_score_nfkb_cond(responder_cells_nfkb_cond,:) ,[para_rand(responder_cells_nfkb_cond,:),params_cond(responder_cells_nfkb_cond,:)] );

    fold_zoom_1 = 0.6;
    paperpos=[0,0,180,100];
    papersize=[180 100];
    draw_pos=[10,10,160,80];
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)

    b = bar(max(KS,[],2),'FaceColor','flat');

    % Default color for all bars
    numBars = size(b.CData,1);   % number of bars
    b.CData = repmat([0 0.4470 0.7410], numBars, 1);
    
    % First bar = gray
    b.CData(1,:) = [0.5 0.5 0.5];

    hold on;
    XL = xlim();
    plot(XL,[max(KS(1,:),[],2),max(KS(1,:),[],2)],'r--') % plot the base line for neg control 
    ylim([0,0.6])
    saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_PAWN_Sens_para_',metric_name,'_withlabel'),'pdf')
    set(gca,'YTick',0:0.2:0.6,'XTickLabel',{},'YTickLabel',{})

    saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_PAWN_Sens_para_',metric_name),'pdf')
    close()

    % using descrete situation:
    % synergy_or_not = -2 * (synergy_score_nfkb_cond<=1 ) + 1* (synergy_score_nfkb_cond>= 1.25 );
    % figure(3)
    % [KS_descrete] = PAWN_given_data(synergy_or_not(responder_cells_nfkb_cond,:) ,params_cond(responder_cells_nfkb_cond,:) );
    % bar(max(KS_descrete,[],2))

end



% update name, move to /src
% relpace/ delete / or move to debug: the original codes single_cell_plot_script.m
if 1
    %% antag condition: under 5 by 5 condtion

    switch metric_name
        case'peak'
            i_row_plot_cond = 4;
            i_col_plot_cond = 5;
        case 'total'
            i_row_plot_cond = 5;
            i_col_plot_cond = 5;
    end
    
    responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

    antag_score_nfkb_cond = Antag_score(:,i_row_plot_cond,i_col_plot_cond);

    i_cond = i_row_plot_cond * 6 + i_col_plot_cond + 1;

    params_cond = gene_info_all{i_cond}.gene_parameter_value_vec_genotype{1}{1}';
    params_names = gene_info_all{i_cond}.parameter_name_vec{1};

    % remove extra useless parameters: first p52n2, overwrite by second p52n2
    params_cond = params_cond(:,[1:2,4:end]);
    params_names = params_names([1:2,4:end]);

    params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);
    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));

    % generate random variable for control
    rng(16);
    para_rand = rand(1000,1);

    % PAWN analysis ,interval_bounds, ft
    [KS] = PAWN_given_data(antag_score_nfkb_cond(responder_cells_nfkb_cond,:) ,[para_rand(responder_cells_nfkb_cond,:),params_cond(responder_cells_nfkb_cond,:)] );

    fold_zoom_1 = 0.6;
    paperpos=[0,0,180,100];
    papersize=[180 100];
    draw_pos=[10,10,160,80];
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos*fold_zoom_1,'PaperSize', papersize*fold_zoom_1,'Position',draw_pos*fold_zoom_1)

    b = bar(max(KS,[],2),'FaceColor','flat');

    % Default color for all bars
    numBars = size(b.CData,1);   % number of bars
    b.CData = repmat([0 0.4470 0.7410], numBars, 1);
    
    % First bar = gray
    b.CData(1,:) = [0.5 0.5 0.5];

    hold on;
    XL = xlim();
    plot(XL,[max(KS(1,:),[],2),max(KS(1,:),[],2)],'r--') % plot the base line for neg control 

    ylim([0,0.6])
    saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_PAWN_Sens_para_',metric_name,'_withlabel'),'pdf')
    set(gca,'YTick',0:0.2:0.6,'XTickLabel',{},'YTickLabel',{})

    saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_PAWN_Sens_para_',metric_name),'pdf')
    close()

    % using descrete situation:
    % antag_or_not = -1 * (antag_score_nfkb_cond>=1.1 ) + 1* (antag_score_nfkb_cond<= 0.9 );
    % figure(3)
    % [KS_descrete] = PAWN_given_data(antag_or_not(responder_cells_nfkb_cond,:) ,params_cond(responder_cells_nfkb_cond,:) );
    % bar(max(KS_descrete,[],2))

end


