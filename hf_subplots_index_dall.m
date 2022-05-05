function [c_eclipse_all,c_non_eclipse_all]= hf_subplots_index_dall(eclipse,non_eclipse,...
    Dall,var_dag,ts_dag,pdist_eclipse,pdist_non_eclipse,tl_eclipse,...
    t,t_eclipse,t_total_e_start,t_total_e_end,t_part_e_start,t_part_e_end,...
    preds_all,keval,kevaln,p,q,tt1,tt2,k_to_plot)

%% Function to plot 4 subplots of time series from clustering analysis result
%% By Leila Hernandez Rodriguez, 2021

% Inputs: 
% tt1, tt2 = Initial and final datetime to plot in the x-axis

% Outputs:
% c_eclipse_all = Clusters colors for the eclipse, August 21, 2017
% c_non_eclipse_all = Clusters colors for August 18, 2017

%% Useful parameters
D = {eclipse{1:end-1}}; 
DD = {non_eclipse{1:end-1}};
n = var_dag;    % Number of variables
m = ts_dag;   % Number of time steps
sizelines = 1; % Size of the lines to plot
sizelines_eclipse = 1.5; % Size of the lines to plot
sizelines_noneclipse = 1.5;%0.5
transpline_eclipse = 0.8;
transpline_noneclipse=0.8;%0.25
patch_transp =0.7;             % Transparent area for partial eclipse period
n_dags = length(D);  % Number of dags for the eclipse day
full_dag = n_dags+1; % Fully connected DAG

marker_eclipse='-or';
marker_non_eclipse='-ob';

% Labels and ticks font size
label_fontsize=18;
axis_fontsize=16;

%% 4 plots to show clustered DAGs as time series

%% a) Distance between consecutive DAGs
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
set(gca,'FontSize',11)
Y = pdist_eclipse; % Aug 21
Y_18 = pdist_non_eclipse; % Aug 18
Ymin = min([min(Y),min(Y_18)]);
Ymax = max([max(Y),max(Y_18)]);
x=datenum(t);

% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse;  
hold on

% Plot a gray area for the partial eclipse period
xarea = [t_part_e_start t_part_e_end t_part_e_end t_part_e_start];
yarea = [Ymin Ymin Ymax Ymax];
patch(xarea,yarea,[0.9 0.9 0.9],'FaceAlpha',patch_transp,'LineStyle','none')

% % Plot a gray area for the total eclipse period (2 minutes and 40 seconds)
% xarea = [t_total_e_start t_total_e_end t_total_e_end t_total_e_start];
% yarea = [Ymin Ymin Ymax Ymax];
% patch(xarea,yarea,[0.5 0.5 0.5],'FaceAlpha',patch_transp,'LineStyle','none')

% plot lines and markers
p1 = plot(x, Y,marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse;  
hold on

% Plot scatter plots to show the clusters 
yp = preds_all(:,1:(length(D)-1),keval(kevaln))';  % non_eclipse
%gg = scatter(x,Y_18,yp);
color_k_to_plot=colormap(lines(k_to_plot)); % This line defines the color to all plot! IMPORTANT!PARULA, jet, hsv
gg= gscatter(x,Y_18,yp,color_k_to_plot(1:k_to_plot,:)); % Shows automatically the legend
hold on
ypp = preds_all(:,(length(D)+1):(length(Dall)-1),keval(kevaln))'; % eclipse
ggg = gscatter(x,Y,ypp,color_k_to_plot(1:k_to_plot,:));
ax = gca;
ax.FontSize = axis_fontsize; 
hold on

c_non_eclipse_all = get(gg,'Color'); % Colors assigned to each cluster - eclipse
c_eclipse_all = get(ggg,'Color'); % Colors are assigned to each cluster 

% Plot a vertical line to plot the peak of the eclipse
plot([t_eclipse t_eclipse],[Ymin Ymax],'--r','LineWidth',sizelines)
datetick('x', 'HH:MM')
hold off
grid on
xlim([tt1 tt2])
ylim([Ymin Ymax])
ylabel('\textit{Distance}','Interpreter','latex','FontSize', label_fontsize)
xlabel('')
title('\textit{Distance between consecutive DAGs}','Interpreter','latex','FontSize', label_fontsize)
dim = [.04 .04 .9 .9];
str = {'p =',num2str(p),'q =',num2str(q)};
legend('boxoff')  

%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%saveas(gcf, 'distances_consecutive.png');

%% b) Distance between a sparse DAG and the fully connected DAG 
% Estimate distances between each DAG and a theoretical fully connected DAG

for j1 = 1:n_dags
    pdist_mat(j1) = dag_dist(eclipse{j1},eclipse{full_dag},p,q,n,m);
    pdist_mat_18(j1) = dag_dist(non_eclipse{j1},non_eclipse{full_dag},p,q,n,m);
end

% Plot distances
%f = figure;  
%f.Position = [10 10 1000 300];
subplot(2,2,2)
set(gca,'FontSize',11)
Y = pdist_mat(1:length(D)-1);
Y_18 = pdist_mat_18(1:length(D)-1);
Ymin = min([min(Y),min(Y_18)]);
Ymax = max([max(Y),max(Y_18)]);

% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse; 
hold on

% Plot a gray area for the partial eclipse period
xarea = [t_part_e_start t_part_e_end t_part_e_end t_part_e_start];
yarea = [Ymin Ymin Ymax Ymax];
patch(xarea,yarea,[0.9 0.9 0.9],'FaceAlpha',patch_transp,'LineStyle','none')

% % Plot a gray area for the total eclipse period (2 minutes and 40 seconds)
% xarea = [t_total_e_start t_total_e_end t_total_e_end t_total_e_start];
% yarea = [Ymin Ymin Ymax Ymax];
% patch(xarea,yarea,[0.5 0.5 0.5],'FaceAlpha',patch_transp,'LineStyle','none')

% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse; 
hold on

% Plot scatter plots to show the clusters 
yp = preds_all(:,1:(length(D)-1),keval(kevaln))';
gscatter(x,Y_18,yp,color_k_to_plot(1:k_to_plot,:))
hold on
ypp = preds_all(:,(length(D)+1):(length(Dall)-1),keval(kevaln))'; 
gscatter(x,Y,ypp,color_k_to_plot(1:k_to_plot,:))
ax = gca;
ax.FontSize = axis_fontsize; 
hold on

% Plot a vertical line to plot the peak of the eclipse
plot([t_eclipse t_eclipse],[Ymin Ymax],'--r','LineWidth',sizelines)
datetick('x', 'HH:MM')
hold off
grid on
xlim([tt1 tt2])
ylim([Ymin Ymax])
xlabel('')
ylabel('\textit{Distance}','Interpreter','latex','FontSize', label_fontsize)
title('\textit{Distance of a DAG from fully connected DAG}','Interpreter','latex','FontSize', label_fontsize)

%dim = [.92 .92 .005 .005];
%str = {'p =',num2str(p),'q =',num2str(q)};
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%saveas(gcf, 'distances_DAG.png');

%% c) Degree of relative change, Dr

% Estimate Index (The fraction of change in number of links in B=DAG(t) as compared to A=DAG(t-1))
% n= total number of links in to A=DAG(t-1)
% p1 = pdist_eclipse./tl_eclipse;
% p2 = pdist_non_eclipse./tl_eclipse;

%p1 = 1-(pdist_eclipse./tl_eclipse);      % Prof. Kumar says is more
%meaninful deleting the 1 minus
%p2 = 1-(pdist_non_eclipse./tl_eclipse);

p1 = (pdist_eclipse./tl_eclipse);
p2 = (pdist_non_eclipse./tl_eclipse);

% Plot distances
%f = figure;  
%f.Position = [10 10 1000 500];
subplot(2,2,3)
set(gca,'FontSize',11)
Y = p1;
Y_18 = p2;
Ymin = min([min(Y),min(Y_18)]);
Ymax = max([max(Y),max(Y_18)]);
% Ymin = -1.4; % Leila: this values match the y-lims for the MIWTR case in
% FIgure 9 in the paper
% Ymax = 0.6;
% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse; 
hold on

% Plot a gray area for the partial eclipse period
xarea = [t_part_e_start t_part_e_end t_part_e_end t_part_e_start];
yarea = [Ymin Ymin Ymax Ymax];
patch(xarea,yarea,[0.9 0.9 0.9],'FaceAlpha',patch_transp,'LineStyle','none')

% % Plot a gray area for the total eclipse period (2 minutes and 40 seconds)
% xarea = [t_total_e_start t_total_e_end t_total_e_end t_total_e_start];
% yarea = [Ymin Ymin Ymax Ymax];
% patch(xarea,yarea,[0.5 0.5 0.5],'FaceAlpha',patch_transp,'LineStyle','none')

% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18,marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse; 
hold on

% Plot scatter plots to show the clusters 
yp = preds_all(:,1:(length(D)-1),keval(kevaln))';
gscatter(x,Y_18,yp,color_k_to_plot(1:k_to_plot,:))
hold on
ypp = preds_all(:,(length(D)+1):(length(Dall)-1),keval(kevaln))'; 
gscatter(x,Y,ypp,color_k_to_plot(1:k_to_plot,:))
ax = gca;
ax.FontSize = axis_fontsize; 
hold on

% Plot a vertical line to plot the peak of the eclipse
plot([t_eclipse t_eclipse],[Ymin Ymax],'--r','LineWidth',sizelines)
datetick('x', 'HH:MM')
hold off
grid on
xlim([tt1 tt2])
ylim([Ymin Ymax])
ylabel('$D_{r}[-]$','interpreter','latex','FontSize', label_fontsize)
%title('\textit{Fraction change in number of links}','Interpreter','latex','FontSize', label_fontsize)
title('\textit{Degree of relative change}','Interpreter','latex','FontSize', label_fontsize)
xlabel('')
%str = {'p =',num2str(p),'q =',num2str(q)};
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%saveas(gcf, 'distances_consecutive.png');

%% Degree of connectivity, Dc
e = n*(m*(m-1)/2)+n*(n-1)*(m*(m-1)/2); % By Leila 2021, number of edges for a fully connected DAG with no contenporaneous edges
dc = 1-(pdist_mat./e);   
dc_18 = 1-(pdist_mat_18./e);  

% Plot degree of connectivity
subplot(2,2,4)
set(gca,'FontSize',11)
Y = dc(1:length(D)-1);  % 65 pairs of DAGs
Y_18 = dc_18(1:length(D)-1);  % 65 pairs of DAGs
Ymin = min([min(Y),min(Y_18)]);
Ymax = max([max(Y),max(Y_18)]);

% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse; 
hold on

% Plot a gray area for the partial eclipse period
xarea = [t_part_e_start t_part_e_end t_part_e_end t_part_e_start];
yarea = [Ymin Ymin Ymax Ymax];
patch(xarea,yarea,[0.9 0.9 0.9],'FaceAlpha',patch_transp,'LineStyle','none')

% % Plot a gray area for the total eclipse period (2 minutes and 40 seconds)
% xarea = [t_total_e_start t_total_e_end t_total_e_end t_total_e_start];
% yarea = [Ymin Ymin Ymax Ymax];
% patch(xarea,yarea,[0.5 0.5 0.5],'FaceAlpha',patch_transp,'LineStyle','none')

% plot lines and markers
p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
p1.Color(4) = transpline_eclipse;  
hold on
p2 = plot(x, Y_18,marker_non_eclipse,'LineWidth',sizelines_noneclipse);
p2.Color(4) = transpline_noneclipse; 
hold on

% Plot scatter plots to show the clusters 
yp = preds_all(:,1:(length(D)-1),keval(kevaln))';
gscatter(x,Y_18,yp,color_k_to_plot(1:k_to_plot,:))
hold on
ypp = preds_all(:,(length(D)+1):(length(Dall)-1),keval(kevaln))'; 
gscatter(x,Y,ypp,color_k_to_plot(1:k_to_plot,:))
ax = gca;
ax.FontSize = axis_fontsize; 
hold on

% Plot a vertical line to plot the peak of the eclipse
plot([t_eclipse t_eclipse],[Ymin Ymax],'--r','LineWidth',sizelines)
datetick('x', 'HH:MM')
xlim([tt1 tt2])
ylim([Ymin Ymax])
hold off
grid on
ylabel('$D_{c}[-]$','interpreter','latex','FontSize', label_fontsize)
xlabel('')
title('\textit{Degree of connectivity}','Interpreter','latex','FontSize', label_fontsize)
%[ax1,h1]=suplabel('Clustering eclipse and non-eclipse data','t');
%set(h1,'FontSize',12)

% erase all legends
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')

saveas(gcf, 'hf_subplots_index_final.jpg');


%close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This figure is ONLY useful get the legend from it 
% 
% figure
% set(gca,'FontSize',11)
% Y = dc(1:length(D)-1);  % 65 pairs of DAGs
% Y_18 = dc_18(1:length(D)-1);  % 65 pairs of DAGs
% Ymin = min([min(Y),min(Y_18)]);
% Ymax = max([max(Y),max(Y_18)]);
% 
% % plot lines and markers
% p1 = plot(x, Y, marker_eclipse,'LineWidth', sizelines_eclipse);
% p1.Color(4) = transpline_eclipse;  
% hold on
% p2 = plot(x, Y_18, marker_non_eclipse,'LineWidth',sizelines_noneclipse);
% p2.Color(4) = transpline_noneclipse; 
% hold on
% 
% % Plot scatter plots to show the clusters (THIS CHANGES! INCLUDE THE ALL DATA RESULTS HERE)
% yp = preds_all(:,1:(length(D)-1),keval(kevaln))';
% gscatter(x,Y_18,yp,color_k_to_plot(1:k_to_plot,:))
% hold on
% ypp = preds_all(:,(length(D)+1):(length(Dall)-1),keval(kevaln))'; 
% gscatter(x,Y,ypp,color_k_to_plot(1:k_to_plot,:))
% hold on
% %print(gcf, 'figure.jpg', '-jpg', '-r600');
% saveas(gcf, 'hf_subplots_index_final_legend.jpg');

%close all

