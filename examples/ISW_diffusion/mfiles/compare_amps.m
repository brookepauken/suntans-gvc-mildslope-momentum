clear all; close all; 

modifiers = {'sigtest_10T_matchdz_superbee400', 'sigtest_10T_matchdz_QUICKsharp400', 'hybrid_10T_dzf3_qs_newgrid_200'};
modifiers = {'sigtest_10T_matchdz_superbee400', 'sigtest_10T_matchdz_QUICKsharp400', 'isotest_10T_dzf3_QUICKsharp400', 'isotest_10T_dzf6_QUICKsharp400'};
modifiers = {'sigtest_10T_matchdz_superbee400', 'sigtest_10T_matchdz_QUICKsharp400', 'hybrid_10T_dzf3_qs_newgrid_densfix200', 'hybrid_10T_dzf6_qs_newgrid_densfix200', 'hybrid_10T_dzf3_superbee_newgrid_densfix200'};

modifiers = {'copy_sigma_10T_tvd3200', 'copy_sigma_10T_tvd6200', 'hyrbid_10T_tvd3_dzf3200', 'hyrbid_10T_tvd6_dzf6200', 'copy_hybrid_10T_200_qs_uniformfix200', 'hybrid_10T_800_uniform_qs_smallestamp800'};
modifiers = {'copy_sigma_10T_tvd3200', 'copy_sigma_10T_tvd6200', 'hyrbid_10T_tvd3_dzf3200', 'hyrbid_10T_tvd6_dzf6200'};
%modifiers = {'hyrbid_10T_tvd3_dzf3200', 'hyrbid_10T_tvd6_dzf6200', 'copy_hybrid_10T_200_qs_uniformfix200', 'hybrid_10T_800_uniform_qs_smallestamp800'};
%modifiers = {'hyrbid_10T_tvd6_dzf6200', 'copy_hybrid_10T_200_qs_uniformfix200', 'copy_hybrid_10T_200_qs_uniformfix_quickinterp200', 'copy_hybrid_10T_200_nonuniformfix_quick200'};

modifiers = {'copy_sigma_10T_tvd6200', 'hyrbid_10T_tvd6_dzf6200'};


figure(1)



for filenum = 1:length(modifiers)
    loadtitle = [modifiers{filenum}, '__amps.mat'];
    load(loadtitle);

    timesave = timesave(1:1048);
    amps_max = amps_max(1:1048);

    if(filenum==1)
        ampsmax_1 = amps_max;
    elseif(filenum==2)
        ampsmax_2 = amps_max;
    elseif(filenum==3)
        ampsmax_3 = amps_max;
    else
        ampsmax_4 = amps_max;
    end

    figure(1)
    plot(timesave(amps_max>0), amps_max(amps_max>0))
    hold on

    figure(2)
    if(filenum<2)
    %if(filenum<3)
        plot(timesave(amps_max>0), (amps_max(amps_max>0)-amps_max(1))/amps_max(1)*100, '-o', 'Color', '#D95319', 'MarkerIndices',1:10:1048, 'MarkerSize', 3, 'LineWidth', 1.25)
    else
        plot(timesave(amps_max>0), (amps_max(amps_max>0)-amps_max(1))/amps_max(1)*100, '-', 'Color', '#7E2F8E', 'LineWidth', 1.5)
    end
    hold on



end
% 
% figure(1)
% grid on
% % legend('Sigma Nx400 Superbee', 'Sigma Nx400 Quick and Sharp', 'Hybrid Nx200 Superbee', 'FontSize', 14, 'Location', 'Best')
% % legend('Sigma Nx400 Superbee', 'Sigma Nx400 Quick and Sharp', 'Iso Nx400 dzf3', 'Iso Nx400 dzf6', 'FontSize', 14, 'Location', 'Best')
% legend('Sigma Nx400 Superbee', 'Sigma Nx400 Quick and Sharp', 'Hybrid Nx200 dzf3 QS', 'Hybrid Nx200 dzf6 QS', 'Hybrid Nx200 dzf3 Superbee',  'FontSize', 14, 'Location', 'Best')
% xlabel('Time (periods)', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('Amplitude (m)', 'Interpreter', 'latex', 'FontSize', 14)


figure(2)
grid on
% legend('Sigma Nx400 Superbee', 'Sigma Nx400 Quick and Sharp', 'Hybrid Nx200 Superbee', 'FontSize', 14, 'Location', 'Best')
%legend('Hybrid w/ Superbee', 'Hybrid w/ Quick and Sharp', 'Hybrid updated BCS', 'Hybrid 800 smaller amp', 'FontSize', 14, 'Location', 'Best')
legend('$\sigma$, Superbee', '$\sigma$, QUICK/SHARP', 'Hybrid, Superbee', 'Hybrid, QUICK/SHARP', 'FontSize', 18, 'Location', 'Best')
legend('$\sigma$, QUICK/SHARP', 'Hybrid, QUICK/SHARP', 'FontSize', 18, 'Location', 'Best')
%legend('Hybrid w/ Quick and Sharp', 'Hybrid w/ Quick and Sharp and updated BCS', 'Quick interp', 'Nonuniform all quick', 'FontSize', 14, 'Location', 'Best')
xlabel('Time (periods)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Percent change in amplitude', 'Interpreter', 'latex', 'FontSize', 14)
%xlim([0, 4])
ylim([-2, 4])

set(gca,'fontsize', 20) 


% figure(3)
% clf;
% grid on
% grid on
% plot(timesave, abs(ampsmax_2(1:length(ampsmax_4)) - ampsmax_4))
% xlabel('Time (periods)', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('abs(Updated BC amp - Quick interp amp) (m)', 'Interpreter', 'latex', 'FontSize', 14)
