function hf_dunn_index_dags_subplots(ks,dunn_ylim,dunn_index_non_eclipse_,dunn_index_eclipse_,dunn_index_all,dunn_k_to_plot_final);
%% Function to plot the Dunn Index for clusters of DAGs
% Dunn Index: The ratio of the minimum of inter-cluster distances and the
% maximum of intra-cluster distance.

% By Leila Hernandez Rodriguez, 2021

linea = 1.5;  
font_size=12;

%figure('units','normalized','outerposition',[0 0 1 1])
figure('Position', [100, 100, 700, 300])
subplot(1,3,1)
plot(ks,dunn_index_non_eclipse_,'-*b','LineWidth',linea)
xlabel('$k$','Interpreter','latex','FontSize', font_size)
ylabel('$\textit{Clear-sky day, August 18,2017}$','Interpreter','latex','FontSize', font_size)
xlim([1 length(ks)])
ylim([0 dunn_ylim])
grid on

subplot(1,3,2)
plot(ks,dunn_index_eclipse_,'-*r','LineWidth',linea)
xlabel('$k$','Interpreter','latex','FontSize', font_size)
ylabel('$\textit{Eclipse day, August 21,2017}$','Interpreter','latex','FontSize', font_size)
xlim([1 length(ks)])
ylim([0 dunn_ylim])
grid on

subplot(1,3,3)
plot(ks,dunn_index_all,'-*m','LineWidth',linea)
xlabel('$k$','Interpreter','latex','FontSize', font_size)
ylabel('$\textit{Clear-sky & eclipse days}$','Interpreter','latex','FontSize', font_size)
xlim([1 length(ks)])
ylim([0 dunn_ylim])
hold on
h=plot(dunn_k_to_plot_final,dunn_index_all(dunn_k_to_plot_final),'-mo','MarkerSize',8,'MarkerFaceColor','k');
hold off
leg = legend(h,'$k$ for k-means','Location', 'Best','FontSize', 20);
set(leg,'Interpreter','latex')

grid on

[ax1,h1]=suplabel('Dunn index','t');
set(h1,'FontSize',12,'Interpreter','latex')
saveas(gcf, 'hf_dunnindex_final.jpg');

