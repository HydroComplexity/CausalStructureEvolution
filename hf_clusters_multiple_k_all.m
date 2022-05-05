function hf_clusters_multiple_k_all(tc,D,x_datapoints,y_preds,keval,color_dags,min_n_clusters,k_to_plot,t_eclipse,sizelines)

%% Function to plot clustered DAGs for multiple k
%% By Leila Hernandez Rodriguez, 2021

% Plot clusters in time-series plots for August 21 and August 18 clusters 
% together, for various values of k
% Also, plot k=3, for clear sky day and eclipse day in two plots

label_fontsize=20;
axis_fontsize=14;

figure('units','normalized','outerposition',[0 0 1 1])
y_clusters = min_n_clusters;

for j = 1:length(keval)
    subplot(2,2,j)
    x = datenum(tc);                      
    y = y_preds(:,:,keval(j))';
    y1 = y(1:length(D));
    y2 = y(length(D)+1:end);
    s(1) = scatter(x,y1,'filled'); % non_eclipse 
    s(1).MarkerFaceColor = 'b';
    ylim([1 keval(j)])
    hold on
    s(2) = scatter(x,y2,'filled'); % eclipse
    s(2).MarkerFaceColor = 'r';
    datetick('x','HH:MM')
    hold on
    plot([t_eclipse t_eclipse],[1 keval(j)],'--r','LineWidth',sizelines);
    ylabel('\textit{Cluster ID}','Interpreter','latex','FontSize', 50)
    ax = gca;
    ax.FontSize = axis_fontsize; 
    caption = sprintf('$k = %d$', keval(j));
    title(caption, 'FontSize', 11,'Interpreter','latex','FontSize', label_fontsize); 
    set(gca, 'YTick', 0:y_clusters)
    y_clusters=y_clusters+1;
    grid on   
end
%legend(s,'non-eclipse','eclipse','Location','southoutside','Orientation','horizontal')

%[ax1,h1]=suplabel('Clusters of all DAGs for August 18 and 21','t');
%set(h1,'FontSize',12)
saveas(gcf, 'hf_clusters_multiple_k_all.jpg');

% Plot k=3, for clear sky day and eclipse day in two plots
y_clusters = min_n_clusters;
figure('Position', [100, 100, 800, 250])
x = datenum(tc);                      
y = y_preds(:,:,keval(1))';
y1 = y(1:length(D));
y2 = y(length(D)+1:end);
subplot(1,2,1)
s(1) = scatter(x,y1,'filled'); % non_eclipse 
s(1).MarkerFaceColor = 'b';
ylim([1 keval(1)])
datetick('x','HH:MM')
xlim([min(x) max(x)])
hold on
plot([t_eclipse t_eclipse],[1 keval(1)],'--r','LineWidth',sizelines);
grid on
ylabel('\textit{Cluster ID}','Interpreter','latex','FontSize', 11)
set(gca, 'YTick', 0:k_to_plot)
hold on
subplot(1,2,2)
s(2) = scatter(x,y2,'filled'); % eclipse
s(2).MarkerFaceColor = 'r';
ylim([1 keval(1)])
datetick('x','HH:MM')
xlim([min(x) max(x)])
hold on
plot([t_eclipse t_eclipse],[1 keval(1)],'--r','LineWidth',sizelines);
grid on
ylabel('\textit{Cluster ID}','Interpreter','latex','FontSize', 11)
set(gca, 'YTick', 0:k_to_plot)
saveas(gcf, 'hf_clusters_k3.jpg');
%ax = gca;
%ax.FontSize = axis_fontsize; 
%caption = sprintf('$k = %d$', keval(a));
%title(caption, 'FontSize', 11,'Interpreter','latex','FontSize', label_fontsize); 
   
%legend(s,'non-eclipse','eclipse','Location','southoutside','Orientation','horizontal')

%[ax1,h1]=suplabel('Clusters of all DAGs for August 18 and 21','t');
%set(h1,'FontSize',12)
%saveas(gcf, 'hf_clusters_k3.jpg');

% Plot another figure only to get the legend from there!
% figure('Position', [100, 100, 600, 300])
% s(1) = scatter(x,y1,'filled'); % non_eclipse 
% s(1).MarkerFaceColor = 'b';
% hold on
% s(2) = scatter(x,y2,'filled'); % eclipse
% s(2).MarkerFaceColor = 'r';
% leg = legend(s,'$\textit{non-eclipse}$','$\textit{eclipse}$','Location','southoutside','Orientation','horizontal');
% set(leg,'Interpreter','latex')
% saveas(gcf, 'hf_clusters_multiple_k_all_legend.jpg');