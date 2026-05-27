% if decided to publish, to move to src.

%% plot score heatmap
% no need to plot scores for other species, not informative

figure(i_metric_name)

paperpos = [0,0,500,300];
papersize = [500,300];

draw_pos=[10,10,480,280];

set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)

subplot(total_species,total_pair, (i_species-1) * total_pair + i_pair)

if synergy_1_antag_0
    ratio_plot = squeeze(synergy_ratio(i_species,:,:));
    cell_num_plot = squeeze(synergy_responder_cell_num(i_species,:,:));

    color_limits = [-0.15 0.15];
else
    ratio_plot = squeeze(antag_ratio(i_species,:,:));
    cell_num_plot = squeeze(antag_responder_cell_num(i_species,:,:));

    color_limits = [-0.3 0.3];
end


if devided_responder_or_all_cells

    n = 13;
    xq = linspace(-0.3, 0.3, n);
    % define key colors at –0.3, 0, +0.3
    keyX  = [-0.3 , 0 , 0.3];
    keyRGB = [1  0  0    % blue
        1  1  1    % white
        0  0  1 ]; % red
    % interpolate each channel
    cmap = interp1(keyX, keyRGB, xq, 'linear');

    % create heatmap
    % h = heatmap(ratio_plot, ...
    %     'Colormap', cmap, ...
    %     'ColorLimits', [-0.3 0.3]);

    % figure(2)
    % subplot(total_species,total_pair, (i_species-1) * total_pair + i_pair)
    X_loc = ones(5,1) * [1:5]
    Y_loc = (ones(5,1) * [5:-1:1] )';
    bc = bubblechart(X_loc(:),Y_loc(:),cell_num_plot(:),ratio_plot(:));
    colormap(gca,cmap)
    clim(color_limits)
    bc.MarkerEdgeColor = [0.3 0.3 0.3];
    bc.LineWidth = 0.25;
    bubblesize([1 58]*0.1)
    bubblelim([10,500])

    colorbar off
    set(gca,'XTickLabel',{},'YTickLabel',{})
    xlim([0.45,5.55])
    ylim([0.45,5.55])

    if 0 % generate figure legends

        subplot(total_species,total_pair, (i_species-1) * total_pair + i_pair)

        X_loc = ones(5,1) * [1:5];
        Y_loc = (ones(5,1) * [5:-1:1] )';
        cell_num_plot = ones(5,1) * [5,5,100,250,500];
        ratio_plot = ones(5,1) * [0,0,0,0,0];
        bc = bubblechart(X_loc(:),Y_loc(:),cell_num_plot(:),ratio_plot(:));
        colormap(gca,cmap)
        clim(color_limits)
        bc.MarkerEdgeColor = [0.3 0.3 0.3];
        bc.LineWidth = 0.25;
        bubblesize([1 58]*0.1)
        bubblelim([10,500])

        colorbar off
        set(gca,'XTickLabel',{},'YTickLabel',{})
        xlim([0.45,5.55])
        ylim([0.45,5.55])
        saveas(gcf,strcat(save_fig_path,'synergy_antag_ratio_bubble_size_legends',metric_name),'pdf')


        % pause here to tune which size to generate
        a = 1;
    end
end


% h.XDisplayLabels = repmat("", 1, 5);
% h.YDisplayLabels = repmat("", 1, 5);


