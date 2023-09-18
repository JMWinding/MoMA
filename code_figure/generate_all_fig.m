% matversion = "author";
matversion = "temp";

%%
cenoteFinal = 221;
cenoteOffset = cenoteFinal - 1;

%% 
% figure 6a,6b
figure_results_mainbar_nTx_thpt
% figure 8
figure_results_Lp_thpt
% figure 9
figure_results_missing
% figure 10 (include reference index)
figure_results_codecomplement
% figure 11
figure_results_ceLoss_nTx_ce % DIFFERENT
% figure 12a,12b
figure_results_moreMo_Tx_ce % DIFFERENT
% figure 13
figure_results_codetuple % DIFFERENT
% figure 14
figure_results_moreMo_rate_pd
% figure 15
figure_results_moreMo_pdorder

%% debug
if false
    % 
    figure_results_ceLoss_nTx_pdcluster
    % 
    figure_results_moreMo_Tx_pdcluster
end