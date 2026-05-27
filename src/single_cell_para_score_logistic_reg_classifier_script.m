
fold_zoom_2_class = 1.5;%1.5
fold_zoom_3 = 1.5;%2

paperpos=[0,0,100,100];
papersize=[100 100];
draw_pos=[10,10,90,90];

paperpos_wide=[0,0,120,100];
papersize_wide=[120 100];
draw_pos_wide=[10,10,110,90];

if 1
    %% synergy condition: under 5 by 5 condtion

    figure(1)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos*fold_zoom_2_class,'PaperSize', papersize*fold_zoom_2_class,'Position',draw_pos*fold_zoom_2_class)


    switch metric_name
        case'peak'
            i_row_plot_cond = 2;
            i_col_plot_cond = 1;

        case 'total'
            i_row_plot_cond = 2;
            i_col_plot_cond = 2;
    end
    i_cond = i_row_plot_cond * 6 + i_col_plot_cond + 1;
    params_cond = gene_info_all{i_cond}.gene_parameter_value_vec_genotype{1}{1}';
    params_names = gene_info_all{i_cond}.parameter_name_vec{1};
    params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);
    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
    responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

    % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
    synergy_score_nfkb_cond = Synergy_score(:,i_row_plot_cond,i_col_plot_cond);
    plot_cell_index_cond = find(synergy_score_nfkb_cond>1.25 & responder_cells_nfkb_cond);
    plot_cell_index_cond = find(synergy_score_nfkb_cond>1.25 & responder_cells_nfkb_cond);
    plot_cell_index_cond = find(synergy_score_nfkb_cond>1.25 & responder_cells_nfkb_cond);


    if 1
        %% classifier on non-synergistic (<1) vs synergistic (>1.25);
        data_for_regression_all = [params_cond(responder_cells_nfkb_cond,:),synergy_score_nfkb_cond(responder_cells_nfkb_cond)];
        index_for_non_speci_vs_speci = data_for_regression_all(:,end)>1.25 | data_for_regression_all(:,end)<1;
        data_for_classifier = [data_for_regression_all(index_for_non_speci_vs_speci,1:end-1),...
            data_for_regression_all(index_for_non_speci_vs_speci,end)>1.25];

        % Load data (10 columns: 9 features + 1 label)
        data = data_for_classifier;

        X = data(:, 1:end-1);   % Features
        X = X(:,[1:2,4:end]);

        rescale_or_not = 1;
        X_rescale_new = zeros(size(X));

        if rescale_or_not
            X_rescale = log(X);
            for i_para = 1:8
                % figure(1)
                % subplot(2,4,i_para)
                % histogram(X_rescale(:,i_para))

                X_rescale_new(:,i_para) = (X_rescale(:,i_para) - min(X_rescale(:,i_para)))...
                    /(max(X_rescale(:,i_para)) - min(X_rescale(:,i_para)));
                % figure(2)
                % subplot(2,4,i_para)
                % histogram(X_rescale_new(:,i_para))

            end
            %             figure(1)
            % close()
            % figure(2)
            % close()
        end

        X = X_rescale_new;

        y = data(:, end);    % Labels (0 or 1)

        rng(42);  % for reproducibility
        % rng(21)
        cv = cvpartition(y,'HoldOut',0.2);

        trainIdx = training(cv);   % logical index for 80% train
        testIdx  = test(cv);       % logical index for 20% test

        X_train = X(trainIdx,:);
        y_train = y(trainIdx);

        X_test  = X(testIdx,:);
        y_test  = y(testIdx);

        % Fit logistic regression model
        mdl = fitglm(X_train, y_train, 'Distribution', 'binomial');

        % Predict on same dataset (or split into test set if needed)
        y_prob = predict(mdl, X_test);
        y_pred = round(y_prob);  % Convert probabilities to class labels

        if 1
            % Evaluate
            % confusionchart(y, y_pred);
            % Evaluate using confusion matrix
            figure(1)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_2_class,'PaperSize', papersize*fold_zoom_2_class,'Position',draw_pos*fold_zoom_2_class)

            confusionchart(y_test, y_pred, 'Normalization','row-normalized');
            %set(gca,'XtickLabel',{'non-syn','syn'},'YtickLabel',{'non-syn','syn'});
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_norm_byrow_',metric_name),'pdf')
            close()

            figure(2)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_2_class,'PaperSize', papersize*fold_zoom_2_class,'Position',draw_pos*fold_zoom_2_class)

            confusionchart(y_test, y_pred, 'Normalization','column-normalized');
            % set(gca,'XtickLabel',{'non-syn','syn'},'YtickLabel',{'non-syn','syn'});
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_norm_bycolumn_',metric_name),'pdf')
            close()


            figure(3)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            scale_factor = max(data_for_regression_all(:,[1:2,4:end-1])',[],2)-min(data_for_regression_all(:,[1:2,4:end-1])',[],2);
            % bar(log10(abs(mdl.Coefficients.Estimate(2:9,:))./scale_factor));
            bar(abs(mdl.Coefficients.Estimate(2:9,:)));

            set(gca,'XtickLabel',params_names_label([1:2,4:end]));
            ylabel('Coeffieciens')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_feature_coefficients_',metric_name),'pdf')
            close()

            % figure(3)
            %
            % set(gcf, 'PaperUnits','points')
            % set(gcf, 'PaperPosition', paperpos*fold_zoom_3,'PaperSize', papersize*fold_zoom_3,'Position',draw_pos*fold_zoom_3)
            % bar(-log10(mdl.Coefficients.pValue(2:9,:)));
            % set(gca,'XtickLabel',params_names_label([1:2,4:end]));
            % ylabel('-log10(pvals)')
            % saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_feature_coefficient_pvals_',metric_name),'pdf')
            % close()
        end

        accuracy_vec = zeros(1,size(X,2)+1);
        accuracy_vec(1) = mean(y_test == y_pred);
        accuracy_vec_shuffle = accuracy_vec(1);

        for i_para = 1:size(X,2)
            X_train_rmv_para = X_train(:,setdiff(1:size(X_train,2),i_para));
            X_test_rmv_para = X_test(:,setdiff(1:size(X_test,2),i_para));
            mdl_train_rmv_para = fitglm(X_train_rmv_para, y_train, 'Distribution', 'binomial');
            % Predict on same dataset (or split into test set if needed)
            y_prob_test_rmv_para = predict(mdl_train_rmv_para, X_test_rmv_para);
            y_pred_test_rmv_para = round(y_prob_test_rmv_para);  % Convert probabilities to class labels
            accuracy_vec(i_para+1) = mean(y_test == y_pred_test_rmv_para);

            X_train_shuffle_para = X_train;
            X_train_shuffle_para(:,i_para) = X_train(randperm(size(X_train,1)),i_para);
            X_test_shuffle_para = X_test;
            X_test_shuffle_para(:,i_para) = X_test(randperm(size(X_test,1)),i_para);
            mdl_train_shuffle_para = fitglm(X_train_shuffle_para, y_train, 'Distribution', 'binomial');
            % Predict on same dataset (or split into test set if needed)
            y_prob_train_shuffle_para = predict(mdl_train_shuffle_para, X_test_shuffle_para);
            y_pred_train_shuffle_para = round(y_prob_train_shuffle_para);  % Convert probabilities to class labels
            accuracy_vec_shuffle(i_para+1) = mean(y_test == y_pred_train_shuffle_para);

        end

        if 1
            figure(3)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            bar(accuracy_vec)
            set(gca,'XtickLabel',['all',params_names_label([1:2,4:end])]);
            ylabel('Accuracy')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_feature_scores_',metric_name),'pdf')
            close()

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            bar(accuracy_vec(1) - accuracy_vec(2:end))
            set(gca,'XtickLabel',params_names_label([1:2,4:end]));
            ylabel('Accuracy Difference')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_feature_scores_Acc_diff_',metric_name),'pdf')
            close()
        end

        if 1
            figure(3)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            bar(accuracy_vec_shuffle)
            set(gca,'XtickLabel',['all',params_names_label([1:2,4:end])]);
            ylabel('Accuracy')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_shuffle_feature_scores_',metric_name),'pdf')
            close()

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            bar(accuracy_vec_shuffle(1) - accuracy_vec_shuffle(2:end))
            set(gca,'XtickLabel',params_names_label([1:2,4:end]));
            ylabel('Accuracy Difference')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_2classifier_shuffle_feature_scores_Acc_diff_',metric_name),'pdf')
            close()
        end

        if 0
            %% leaving out the group of parameters
            index_vec = {[1,2,6],[3,4,5],[7,8]};

            for i_para = (size(X,2)+1):(size(X,2)+3)

                X_train_rmv_para = X_train(:,setdiff(1:size(X_train,2),index_vec{i_para-size(X,2)}));
                X_test_rmv_para = X_test(:,setdiff(1:size(X_test,2),index_vec{i_para-size(X,2)}));
                mdl_train_rmv_para = fitglm(X_train_rmv_para, y_train, 'Distribution', 'binomial');
                % Predict on same dataset (or split into test set if needed)
                y_prob_test_rmv_para = predict(mdl_train_rmv_para, X_test_rmv_para);
                y_pred_test_rmv_para = round(y_prob_test_rmv_para);  % Convert probabilities to class labels
                accuracy_vec(i_para+1) = mean(y_test == y_pred_test_rmv_para);

                X_train_shuffle_para = X_train;
                X_test_shuffle_para = X_test;
                for j_para = 1:length(index_vec{i_para-size(X,2)})
                    k_para = index_vec{i_para-size(X,2)}(j_para);
                    X_train_shuffle_para(:,k_para) = X_train(randperm(size(X_train,1)),k_para);

                    X_test_shuffle_para(:,k_para) = X_test(randperm(size(X_test,1)),k_para);
                end
                mdl_train_shuffle_para = fitglm(X_train_shuffle_para, y_train, 'Distribution', 'binomial');
                % Predict on same dataset (or split into test set if needed)
                y_prob_train_shuffle_para = predict(mdl_train_shuffle_para, X_test_shuffle_para);
                y_pred_train_shuffle_para = round(y_prob_train_shuffle_para);  % Convert probabilities to class labels
                accuracy_vec_shuffle(i_para+1) = mean(y_test == y_pred_train_shuffle_para);

            end

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_3,'PaperSize', papersize*fold_zoom_3,'Position',draw_pos*fold_zoom_3)

            bar(accuracy_vec(1) - accuracy_vec(2:end))
            set(gca,'XtickLabel',[params_names_label([1:2,4:end]),'Core','LPSRCP','PamRCP']);
            ylabel('Accuracy Diff')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_grouppara_2classifier_feature_scores_Acc_diff_',metric_name),'pdf')
            close()

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_3,'PaperSize', papersize*fold_zoom_3,'Position',draw_pos*fold_zoom_3)

            bar(accuracy_vec_shuffle(1) - accuracy_vec_shuffle(2:end))
            set(gca,'XtickLabel',[params_names_label([1:2,4:end]),'Core','LPSRCP','PamRCP']);
            ylabel('Accuracy Difference')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_synergy_single_cell_para_grouppara_2classifier_shuffle_feature_scores_Acc_diff_',metric_name),'pdf')
            close()

        end


    end

end

if 1
    %% antag condition: under 5 by 5 condtion


    figure(1)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos*fold_zoom_2_class,'PaperSize', papersize*fold_zoom_2_class,'Position',draw_pos*fold_zoom_2_class)

    switch metric_name
        case'peak'
            i_row_plot_cond = 4;
            i_col_plot_cond = 5;
        case 'total'
            i_row_plot_cond = 5;
            i_col_plot_cond = 5;
    end

    i_cond = i_row_plot_cond * 6 + i_col_plot_cond + 1;
    params_cond = gene_info_all{i_cond}.gene_parameter_value_vec_genotype{1}{1}';
    params_names = gene_info_all{i_cond}.parameter_name_vec{1};
    params_names_label = cellfun(@(s) replace(s, 'params', 'k'), params_names, 'UniformOutput', false);

    %synergy_nfkb_ratio = squeeze(synergy_ratio(1,:,:));
    responder_cells_nfkb_cond = Valid_SNR(:,i_row_plot_cond,i_col_plot_cond);

    % score_to_plot = Synergy_score(Valid_SNR(:,i_row,i_col),i_row,i_col);
    antag_score_nfkb_cond = Antag_score(:,i_row_plot_cond,i_col_plot_cond);
    plot_cell_index_cond = find(antag_score_nfkb_cond<0.9 & responder_cells_nfkb_cond);


    if 1
        %% classifier on non-synergistic (<1) vs synergistic (>1.25);

        data_for_regression_all = [params_cond(responder_cells_nfkb_cond,:),antag_score_nfkb_cond(responder_cells_nfkb_cond)];
        index_for_non_speci_vs_speci = data_for_regression_all(:,end)<0.9 | data_for_regression_all(:,end)>1;
        data_for_classifier = [data_for_regression_all(index_for_non_speci_vs_speci,1:end-1),...
            data_for_regression_all(index_for_non_speci_vs_speci,end)<0.9];

        % Load data (10 columns: 9 features + 1 label)
        data = data_for_classifier;

        X = data(:, 1:end-1);   % Features
        X = X(:,[1:2,4:end]);

        rescale_or_not = 1;
        X_rescale_new = zeros(size(X));
        if rescale_or_not
            X_rescale = log(X);
            for i_para = 1:8
                % figure(1)
                % subplot(2,4,i_para)
                % histogram(X_rescale(:,i_para))

                X_rescale_new(:,i_para) = (X_rescale(:,i_para) - min(X_rescale(:,i_para)))...
                    /(max(X_rescale(:,i_para)) - min(X_rescale(:,i_para)));
                % figure(2)
                % subplot(2,4,i_para)
                % histogram(X_rescale_new(:,i_para))

            end
            %             figure(1)
            % close()
            % figure(2)
            % close()
        end

        X = X_rescale_new;

        y = data(:, end);    % Labels (0 or 1)

        %rng(42);  % for reproducibility
        %rng(21); 14,17, 47, 66, 59, 82,96
        rng(42)
        cv = cvpartition(y,'HoldOut',0.2);

        trainIdx = training(cv);   % logical index for 80% train
        testIdx  = test(cv);       % logical index for 20% test

        X_train = X(trainIdx,:);
        y_train = y(trainIdx);

        X_test  = X(testIdx,:);
        y_test  = y(testIdx);


        % Fit logistic regression model
        mdl = fitglm(X_train, y_train, 'Distribution', 'binomial');

        % Predict on same dataset (or split into test set if needed)
        y_prob = predict(mdl, X_test);
        y_pred = round(y_prob);  % Convert probabilities to class labels

        % Evaluate
        % confusionchart(y, y_pred);
        % Evaluate using confusion matrix
        figure(1)
        confusionchart(y_test, y_pred, 'Normalization','row-normalized');
        %set(gca,'XtickLabel',{'non-syn','syn'},'YtickLabel',{'non-syn','syn'});
        saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_norm_byrow_',metric_name),'pdf')
        close()

        figure(2)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos*fold_zoom_2_class,'PaperSize', papersize*fold_zoom_2_class,'Position',draw_pos*fold_zoom_2_class)

        confusionchart(y_test, y_pred, 'Normalization','column-normalized');
        %set(gca,'XtickLabel',{'non-syn','syn'},'YtickLabel',{'non-syn','syn'});
        saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_norm_bycolumn_',metric_name),'pdf')
        close()


        figure(3)

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

        scale_factor = max(data_for_regression_all(:,[1:2,4:end-1])',[],2)-min(data_for_regression_all(:,[1:2,4:end-1])',[],2);
        %bar(log10(abs(mdl.Coefficients.Estimate(2:9,:))./scale_factor));
        bar(abs(mdl.Coefficients.Estimate(2:9,:)));
        set(gca,'XtickLabel',params_names_label([1:2,4:end]));
        ylabel('Coeffieciens')
        saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_feature_coefficients_',metric_name),'pdf')
        close()


        figure(3)

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)
        bar(-log10(mdl.Coefficients.pValue(2:9,:)));
        set(gca,'XtickLabel',params_names_label([1:2,4:end]));
        ylabel('-log10(pvals)')
        saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_feature_coefficient_pvals_',metric_name),'pdf')
        close()


        accuracy_vec = zeros(1,size(X,2)+1);
        accuracy_vec(1) = mean(y_test == y_pred);
        accuracy_vec_shuffle = accuracy_vec(1);

        for i_para = 1:size(X,2)

            X_train_rmv_para = X_train(:,setdiff(1:size(X_train,2),i_para));
            X_test_rmv_para = X_test(:,setdiff(1:size(X_test,2),i_para));
            mdl_train_rmv_para = fitglm(X_train_rmv_para, y_train, 'Distribution', 'binomial');
            % Predict on same dataset (or split into test set if needed)
            y_prob_test_rmv_para = predict(mdl_train_rmv_para, X_test_rmv_para);
            y_pred_test_rmv_para = round(y_prob_test_rmv_para);  % Convert probabilities to class labels
            accuracy_vec(i_para+1) = mean(y_test == y_pred_test_rmv_para);

            X_train_shuffle_para = X_train;
            X_train_shuffle_para(:,i_para) = X_train(randperm(size(X_train,1)),i_para);
            X_test_shuffle_para = X_test;
            X_test_shuffle_para(:,i_para) = X_test(randperm(size(X_test,1)),i_para);
            mdl_train_shuffle_para = fitglm(X_train_shuffle_para, y_train, 'Distribution', 'binomial');
            % Predict on same dataset (or split into test set if needed)
            y_prob_train_shuffle_para = predict(mdl_train_shuffle_para, X_test_shuffle_para);
            y_pred_train_shuffle_para = round(y_prob_train_shuffle_para);  % Convert probabilities to class labels
            accuracy_vec_shuffle(i_para+1) = mean(y_test == y_pred_train_shuffle_para);

        end
        figure(3)

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

        bar(accuracy_vec)
        set(gca,'XtickLabel',['all',params_names_label([1:2,4:end])]);
        ylabel('Accuracy')
        saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_feature_scores_',metric_name),'pdf')
        close()

        figure(4)

        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

        bar(accuracy_vec(1) - accuracy_vec(2:end))
        set(gca,'XtickLabel',params_names_label([1:2,4:end]));
        ylabel('Accuracy Diff')
        saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_feature_scores_Acc_diff_',metric_name),'pdf')
        close()


        if 1
            figure(3)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            bar(accuracy_vec_shuffle)
            set(gca,'XtickLabel',['all',params_names_label([1:2,4:end])]);
            ylabel('Accuracy')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_shuffle_feature_scores_',metric_name),'pdf')
            close()

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos_wide*fold_zoom_3,'PaperSize', papersize_wide*fold_zoom_3,'Position',draw_pos_wide*fold_zoom_3)

            bar(accuracy_vec_shuffle(1) - accuracy_vec_shuffle(2:end))
            set(gca,'XtickLabel',params_names_label([1:2,4:end]));
            ylabel('Accuracy Difference')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_2classifier_shuffle_feature_scores_Acc_diff_',metric_name),'pdf')
            close()
        end


        if 0
            %% leaving out the group of parameters
            index_vec = {[1,2,6],[3,4,5],[7,8]};

            for i_para = (size(X,2)+1):(size(X,2)+3)

                X_train_rmv_para = X_train(:,setdiff(1:size(X_train,2),index_vec{i_para-size(X,2)}));
                X_test_rmv_para = X_test(:,setdiff(1:size(X_test,2),index_vec{i_para-size(X,2)}));
                mdl_train_rmv_para = fitglm(X_train_rmv_para, y_train, 'Distribution', 'binomial');
                % Predict on same dataset (or split into test set if needed)
                y_prob_test_rmv_para = predict(mdl_train_rmv_para, X_test_rmv_para);
                y_pred_test_rmv_para = round(y_prob_test_rmv_para);  % Convert probabilities to class labels
                accuracy_vec(i_para+1) = mean(y_test == y_pred_test_rmv_para);

                X_train_shuffle_para = X_train;
                X_test_shuffle_para = X_test;
                for j_para = 1:length(index_vec{i_para-size(X,2)})
                    k_para = index_vec{i_para-size(X,2)}(j_para);
                    X_train_shuffle_para(:,k_para) = X_train(randperm(size(X_train,1)),k_para);

                    X_test_shuffle_para(:,k_para) = X_test(randperm(size(X_test,1)),k_para);
                end
                mdl_train_shuffle_para = fitglm(X_train_shuffle_para, y_train, 'Distribution', 'binomial');
                % Predict on same dataset (or split into test set if needed)
                y_prob_train_shuffle_para = predict(mdl_train_shuffle_para, X_test_shuffle_para);
                y_pred_train_shuffle_para = round(y_prob_train_shuffle_para);  % Convert probabilities to class labels
                accuracy_vec_shuffle(i_para+1) = mean(y_test == y_pred_train_shuffle_para);

            end

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_3,'PaperSize', papersize*fold_zoom_3,'Position',draw_pos*fold_zoom_3)

            bar(accuracy_vec(1) - accuracy_vec(2:end))
            set(gca,'XtickLabel',[params_names_label([1:2,4:end]),'Core','LPSRCP','PamRCP']);
            ylabel('Accuracy Diff')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_grouppara_2classifier_feature_scores_Acc_diff_',metric_name),'pdf')
            close()

            figure(4)

            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos*fold_zoom_3,'PaperSize', papersize*fold_zoom_3,'Position',draw_pos*fold_zoom_3)

            bar(accuracy_vec_shuffle(1) - accuracy_vec_shuffle(2:end))
            set(gca,'XtickLabel',[params_names_label([1:2,4:end]),'Core','LPSRCP','PamRCP']);
            ylabel('Accuracy Difference')
            saveas(gcf,strcat(save_fig_path,ligand1,ligand2,'_antag_single_cell_para_grouppara_2classifier_shuffle_feature_scores_Acc_diff_',metric_name),'pdf')
            close()
        end

    end



end