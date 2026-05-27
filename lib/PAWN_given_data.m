function [KS] = PAWN_given_data(y,X)
% ,interval_bounds, ft
%PAWN Run PAWN for Global Sensitivity Analysis of a supplied model
%
% given [X,y] data set, calculate KS for each x_i; with xvals as interval
% bounds

% ## Refernces
% [1]: Pianosi, F., Wagener, T., 2015. A simple and efficient method for
% global sensitivity analysis based on cumulative distribution functions.
% Environ. Model. Softw. 67, 1-11. doi:10.1016/j.envsoft.2015.01.004
% [2]: Distribution-based sensitivity analysis from a generic input-output sample
% Author links open overlay panelFrancesca Pianosi a b
% , Thorsten Wagener a b

y_u = y;
npts = 100;
n_interval = 10; % or 15, or 20

% Find bounds of the model outputs
% m1 = min([y_u']);
% m2 = max([y_u']);
% Evaluate the CDF with kernel density for unconditioned samples
% [f,~] = ksdensity(y_u, linspace(m1,m2,npts), 'Function', 'cdf');
% [f, x] = ecdf(y_u);

% Evaluate the CDF with kernel density for conditioned samples and use
% that to find the KS statistic (Eqn 4 in the paper).
KS = -ones(size(X,2),n_interval);
ft = -ones(size(X,2),npts,n_interval);
interval_bounds = NaN(size(X,2),n_interval+1);

for i_X = 1:size(X,2)

    interval_bounds(i_X,:) = prctile(X(:,i_X),0:10:100);
    % linspace(min(X(:,i_X))*0.999,max(X(:,i_X)),n_interval+1);
    for i_interval=1:n_interval
        lower_bound = interval_bounds(i_X,i_interval);
        up_bound = interval_bounds(i_X,i_interval+1);

        % Temporarily store the current conditioned samples
        y_t = y((X(:,i_X)<= up_bound) & (X(:,i_X)> lower_bound));
        if length(y_t)<10

            KS(i_X,i_interval) = NaN;
        else
            % [ft(i_X,:,i_interval),~] = ksdensity(y_t, linspace(m1, m2,...
            %     npts), 'Function', 'cdf');
            % [ft_tmp,xt] = ecdf(y_t);
            % ft(i_X,:,i_interval) = ft_tmp;
            % KS(i_X,i_interval) = max(abs(ft(i_X,:,i_interval)-f)); % Eqn 4
            [~, ~, KS(i_X,i_interval)] = kstest2(y_u, y_t);

        end
    end
end
end