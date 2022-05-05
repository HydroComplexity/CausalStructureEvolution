%% Function to find the closest DAG of each centroid to the DAG made with data for all data
% incluiding the eclipse day and the non-eclipse day

% By Leila Hernandez Rodriguez, 2021

% For the filled areas in the plot:(NOT USED!) 
% Javier Montalt Tordera (2021). Filled area plot (https://www.mathworks.com/matlabcentral...
%/fileexchange/69652-filled-area-plot), MATLAB Central File Exchange. Retrieved September 10, 2021.

function [dist_k_vs_all_dags,labels,k_number,real_Dall_min_distance,real_Dall_min_index,real_Dall_min,real_centroid_non_eclipse_index,real_centroid_eclipse_index] = hf_kcentroids_closest_to_all_data(k_to_plot,...
centers_all,Dall,Dall_n,p,q,n,m,c_eclipse_all,preds_all,n_dags_per_cluster_all,n_dags)

k = k_to_plot;                    % number of clusters
kk = Dall_n;                      % number of DAGs
k_number = preds_all(:,:,k_to_plot)';     % location/index of clusters
dist_k_vs_all_dags=[];

% Labels and ticks font size
label_fontsize=18;
axis_fontsize=16;

%% FOR ALL DAGs
% Estimate distances between each cluster and all DAGs
for i=1:k;
    for j=1:kk;
       kcenter_one= centers_all(:,:,i);   % Clusters k (1 to 6)
       kcenter_two= Dall{1,j};            % All DAGs from obs (132 DAGs)
       dist_k_vs_all_dags(j,i) = dag_dist(kcenter_one,kcenter_two,p,q,n,m);
       labels{j,i} = sprintf('k_{%d}_{}DAG_{%d}',i,j);      
    end
end

% % Test to find more similar DAGs to clusters in all DAGs.
% % It does not worked! I found that one DAG was the more similar to 3 cluster's
% % centroid...
% min_dist_vs_all_dags=[];
% cluster_min_dist_vs_all_dags=[];  % to verify the cluster number
% index_min_dist_vs_all_dags=[];    % the corresponding index
% for i=1:k;
%     [min_dist_vs_all_dags(1,i),index_min_dist_vs_all_dags(1,i)]=min(dist_k_vs_all_dags(:,i));    
% end

%% For each cluster, find the more similar DAG of all clusters 
% it can be from the eclipse or the non-eclipse day

i2=[];
for i=1:k;
    i2=find(k_number==i);  %for i=1, find the index in DAll of all DAGs in cluster 1
    [a,b]=min(dist_k_vs_all_dags(i2,i)); % find the value and index of the min distance between any DAG and the centroid in cluster 1   
    real_Dall_min_distance(1,i) = a;    % minimum distance between the closest DAG and the centroid for all clusters 
    real_Dall_min_index(1,i) = i2(b);   % index of distance between the closest DAG and the centroid for all clusters 
    real_Dall_min(1,i) = Dall(i2(b));   % vector of real closest DAGs to cluster centroids, eclipe or non-eclipse. 
    
    real_Dall_mean_distance(1,i) = mean(dist_k_vs_all_dags(i2,i)); % find the average distance between any DAG and the centroid in cluster 1      
    
    [c,d]=max(dist_k_vs_all_dags(i2,i)); % find the value and index of the min distance between any DAG and the centroid in cluster 1   
    real_Dall_max_distance(1,i) = c;    % distance between the farthest DAG and the centroid for all clusters 
    real_Dall_max_index(1,i) = i2(d);   % index of distance between the farthest DAG and the centroid for all clusters 
    real_Dall_max(1,i) = Dall(i2(d));   % vector of real farthest DAGs to cluster centroids, eclipe or non-eclipse. 
    
    real_std_Dall(1,i)=std(dist_k_vs_all_dags(i2,i)); % find the value of the standard deviation
      
end
% To find the more similar DAG for the eclipse day and non-eclipse day
% First 66 dags are NON-ECLIPSE
non_eclipse_dist = dist_k_vs_all_dags(1:n_dags,:);
eclipse_dist = dist_k_vs_all_dags(n_dags+1:Dall_n,:);
i3=[];
i4=[];
k_number_non_eclipse=k_number(1:n_dags);
k_number_eclipse=k_number(n_dags+1:Dall_n);

for i=1:k;
    % non_eclipse
    i3=find(k_number_non_eclipse==i);  %for i=1, find the index in DAll of all DAGs in cluster 1
    [aa,bb]=min(non_eclipse_dist(i3,i)); % find the value and index of the min distance between any DAG and the centroid in cluster 1   
    if aa>0
        real_non_eclipse_min_distance(1,i) = aa;    % distance to the closest DAG non-eclipse 
        real_centroid_non_eclipse_index(1,i) = i3(bb);   % index of closest DAG non-eclipse
        real_non_eclipse_dags(1,i) = non_eclipse_dist(i3(bb));   % vector of real closest DAGs to cluster centroids for the non-eclipse. 
    else
        real_non_eclipse_min_distance(1,i) = nan;
        real_non_eclipse_min_index(1,i) = nan;   % index of closest DAG non-eclipse    
    end
    % eclipse
     i4=find(k_number_eclipse==i);  %for i=1, find the index in DAll of all DAGs in cluster 1
     [cc,dd]=min(eclipse_dist(i4,i)); % find the value and index of the min distance between any DAG and the centroid in cluster 1   
%     real_eclipse_min_distance(1,i) = cc;    % distance to the closest DAG eclipse 
%     real_centroid_eclipse_index(1,i) = i4(dd);   % index of closest DAG eclipse
%     real_eclipse_dags(1,i) = Dall(i4(dd));   % vector of real closest DAGs to cluster centroids for the eclipse. 

    if cc>0
        real_eclipse_min_distance(1,i) = cc;    % distance to the closest DAG non-eclipse 
        real_centroid_eclipse_index(1,i) = i4(dd);   % index of closest DAG non-eclipse
        real_eclipse_dags(1,i) = eclipse_dist(i4(dd));   % vector of real closest DAGs to cluster centroids for the non-eclipse. 
    else
        real_eclipse_min_distance(1,i) = nan;
        real_eclipse_min_index(1,i) = nan;   % index of closest DAG non-eclipse    
    end
end

% Plot distances for each cluster
figure('Position', [100, 100, 900, 500])

% Plot standard deviation as error bars
x=1:k_to_plot;
xx= k_to_plot+1; 
y=real_Dall_mean_distance;
ymax= y  + real_std_Dall;
ymin= y - real_std_Dall;

h1 = axes;
errorbar(x,(ymin+ymax)/2,(ymax-ymin)/2,'Color', 'k')
set(h1, 'YAxisLocation', 'Left')
set(h1, 'TickLength',[0 0])
set(h1, 'TickLength',[0 0])
set(h1, 'XTick', 0:k_to_plot+1)
%set(h1, 'Xtick', [])
hold on

% plot the shaded area between the max and min values
curve1 = real_Dall_max_distance;
curve2 = real_Dall_min_distance;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
hstd=fill(x2, inBetween,'k');
set(hstd,'facealpha',.1)
hold on;
plot(x, y, '-.k', 'LineWidth', 2);
hold on;

% This "for loop" is for the filled dots and to use the colors of the clusters using "c_eclipse_all"
for i=1:k_to_plot   
    lp=c_eclipse_all{i,1};
    %scatter(x(i), real_Dall_min_distance(i), 150, lp, 'filled'); #colors
    %from clusters
    scatter(x(i), real_eclipse_min_distance(i), 60, 'r', 'filled');
    scatter(x(i), real_non_eclipse_min_distance(i), 60,'b','filled');
    hold on
end

ax = gca;
ax.FontSize = axis_fontsize; 

hold on

ylabel('$\textit{Distance between centroid and DAGs}$','Interpreter','latex','FontSize', label_fontsize)
grid on
ax = gca;
ax.GridColor = [0.1, 0.1, 0.1]; 
hold on

h2 = axes;
tf = isstruct(n_dags_per_cluster_all);% check if n_dags_per_cluster_all is a vector or a structure
if tf==1;
    yy=n_dags_per_cluster_all(k_to_plot);
    plot(x,yy.ki,'--b','LineWidth',1.5)
else    
    yy=n_dags_per_cluster_all;
    plot(x,yy,'--b','LineWidth',1.5)
end
set(h2, 'YAxisLocation', 'Right')
set(h2, 'XLim', get(h1, 'XLim'))
set(h2, 'Color', 'None')
set(h2, 'Xtick', [])
set(h2, 'YColor', 'b' )
%set(h2,'yticklabels', yy)
box off

ax = gca;
ax.FontSize = axis_fontsize; 

xlabel('$Cluster$','Interpreter','latex','FontSize', label_fontsize,'Position',[-0.5 0.5])
ylabel('$\textit{Number of DAGs}$','Interpreter','latex','FontSize', label_fontsize)
%title('$\textit{Distance between centroid and clustered DAGs}$','Interpreter','latex')
print('-r600', 'hf_kcentroids_closest.jpg', '-djpeg');

% % Create an excel file of the real closest DAGs to plot them in the Jupyter notebook 
% %filename = sprintf('real_cluster_centers_all_k=%d_%s.xls',k,datestr(now,'mm-dd-yyyy HH-MM'));  
% filename = sprintf('real_cluster_centers_all_k_%d.xls',k);  
%   for i=1:k
%       sheetname=sprintf('k%i',k,i);
%       xlswrite(filename,real_Dall_dags{1,i},sheetname);
%   end
% 
% % Create an excel file of the real closest NON-ECLIPSE DAGs to plot them in the Jupyter notebook 
% %filename = sprintf('real_cluster_centers_non_eclipse_k_%d_%s.xls',k,datestr(now,'mm-dd-yyyy HH-MM'));  
% filename = sprintf('real_cluster_centers_non_eclipse_k_%d.xls',k);  
%   for i=1:k
%       sheetname=sprintf('k%i',k,i);
%       xlswrite(filename,real_non_eclipe_dags{1,i},sheetname);
%   end
% 
% % Create an excel file of the real closest ECLIPSE DAGs to plot them in the Jupyter notebook 
% filename = sprintf('real_cluster_centers_eclipse_k_%d_%s.xls',k,datestr(now,'mm-dd-yyyy_HH-MM'));  
% %filename = sprintf('real_cluster_centers_eclipse_k_%d.xls',k);  
%   for i=1:k
%       sheetname=sprintf('k%i',k,i);
%       xlswrite(filename,real_Dall_dags{1,i},sheetname);
%   end
