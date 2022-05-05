%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Causal Structure Evolution:                                                               %
% Evolution of the causal structure of a high-frequency multivariate ecohydrological system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created by Leila C Hernandez Rodriguez, 2021.
% University of Illinois at Urbana-Champaign

% Credits to "Beyond rankings: comparing directed acyclic graphs"
% Malmi, Tatti, and Gionis, 2015. eric.malmi@aalto.fi
% https://github.com/ekQ/dag-comparisons
% Copyright (c) 2015 Eric Malmi

%%
clear all
close all
clc

% NOTE: "USER FRIENDLY" denotes that the variables or parameters can be
% changed by the user. The suggested values are included. 

%% ---------------------------------------------------------------------
%% Preliminaries: 
%% 1. Store the excel files with the adjacency matrices in the location of this file
%% 2. Uncomment this seccion the first time you use the code only to upload the DAGs

% %Load data from August 21, 2017. Aka solar-eclipse-day
% [~,sheet_name]=xlsfinfo('adj_matrix_DAG_2017_08_21_eclipse.xlsx');
% for k=1:numel(sheet_name)
% 
%  eclipse{k}=xlsread('adj_matrix_DAG_2017_08_21_eclipse.xlsx',sheet_name{k});  % Matrix of dags
% end 
% save eclipse
% 
% %Load data from August 18, 2017. Aka Clear-sky-day
% [~,sheet_name]=xlsfinfo('adj_matrix_DAG_2017_08_18.xlsx');
% for k=1:numel(sheet_name)
%  non_eclipse{k}=xlsread('adj_matrix_DAG_2017_08_18.xlsx',sheet_name{k});  % Matrix of dags
% end 
% save non_eclipse

%% ---------------------------------------------------------------------
%% Load the adjacency matrices storaged in excel files in the same folder as this file

%% Load adjacency matrices after MIWTR of the DAGs created with 3-mins data each
load eclipse                           % Load data from August 21 - eclipse day   % USER FRIENDLY
load non_eclipse                       % Load data from August 18 - clear-sky day % USER FRIENDLY

n_restarts = 5;                        % Number of restarts for DAGs (Default: 5) % USER FRIENDLY

% Select data to be analyzed 
D = {eclipse{1:end-1}}; 
DD = {non_eclipse{1:end-1}};
% The last value is the fully connected DAG (useful for some analysis)
D_full = eclipse;
DD_full = non_eclipse;

%% ---------------------------------------------------------------------
%% Define parameters and methods

% Distance measure parameters (Check Malmi et al, 2015, before deciding on the values of p and q)
p = 1.0;  % Parameter in Malmi et al,2015. % USER FRIENDLY
q = 0.0;  % Parameter in Malmi et al,2015. % USER FRIENDLY. More than zero is not an option for a DAG.

% k number of clusters
min_n_clusters = 3;   % Minimum number of clusters % USER FRIENDLY
max_n_clusters = 6;   % Maximum number of clusters % USER FRIENDLY

ks = [1:max_n_clusters]; % For analysis
keval = [min_n_clusters:max_n_clusters]; % For selecting k (do not delete!)

% DAG structure
var_dag = 5;         % Number of variables in a DAG  % USER FRIENDLY 
ts_dag = 11;         % Number of time-steps in a DAG % USER FRIENDLY 
var_name = {'WS','Uz','T','CO2','H2O'}; % Your variable names % USER FRIENDLY 

% Define parameters and methods for k-means cluster 
n_dags = length(D);  % Number of dags for the eclipse day
full_dag = n_dags+1; % Fully connected DAG
L = var_dag*ts_dag;  % Number of nodes in a DAG
e = var_dag*(ts_dag*(ts_dag-1)/2)+var_dag*(var_dag-1)*(ts_dag*(ts_dag-1)/2); % Number of edges in a fully connected DAG with no contenporaneous edges
n = var_dag;    % Number of variables
m = ts_dag;     % Number of time steps
sizelines = 1;  % Size of the lines to plot

%% Plotting settings 
t1 = datetime(2017,08,21,10,46,53); % USER FRIENDLY  
t2 = datetime(2017,08,21,13,58,53); % USER FRIENDLY  
t3 = datetime(2017,08,21,14,01,53); % USER FRIENDLY  
t = t1:minutes(3):t2;
tc = t1:minutes(3):t3;
tt1=datenum(t1);
tt2=datenum(t2);

% To plot lines for the partial eclipse % USER FRIENDLY:
t_eclipse = datenum(datetime(2017,08,21,12,19,53));        % peak of solar eclipse 1219_1222 12:19:53 CDT (94% of obscurity)
% The longest duration of the eclipse was near Carbondale, Illinois, 
% where the sun was completely covered for two minutes and 40 seconds.
% We include a shaded section of 80 seconds before and after the time of
% the peak at the flux towe
t_total_e_start = datenum(datetime(2017,08,21,12,18,33));   % almost totality eclipse begins at 12:18:33 CST
t_total_e_end = datenum(datetime(2017,08,21,12,21,13));     % partial eclipse ends at 12:21:13 CST
% Also, the partial eclipse starts and ends
t_part_e_start = datenum(datetime(2017,08,21,10,52,45));   % partial eclipse starts at 10:52:45 CST
t_part_e_end = datenum(datetime(2017,08,21,13,44,21));     % partial eclipse ends at 13:44:21 CST
sizelines_eclipse = 4;               % Size of the lines to plot
sizelines_noneclipse = 0.5;
transpline_eclipse = 0.8;
transpline_noneclipse=0.3;
%patch_transp =0.9;                  % Transparent area for partial eclipse period

%% DAG Clustering 
% Method for DAG agregation and clustering
methods = {'greedy'};   % k-means type algorithm which updates the cluster centroids (Malmi et al, 2015)
max_iter = 10;          % USER FRIENDLY % Maximum number of iterations (Default: 10)

%% ---------------------------------------------------------------------
%% DEAR USER: DO NOT EDIT FROM HERE BELOW %%
%--------------------------------------------------------------------------------------------------------
% Clustering for the eclipse day - August 21, 2017
[preds_eclipse,res_eclipse,avg_c_n_eclipse,pred_eclipse,err_eclipse,n_iter_eclipse,is_removed_eclipse,centers_eclipse,n_dags_per_cluster_eclipse,...
intra_cluster_eclipse,intra_clustery_eclipse,clusters_centers_eclipse]= dag_hf_func_clusters_v2(n_dags,ks,methods,n_restarts,D,p,q,var_dag,ts_dag,...
max_iter,[],[],false);

% preds = List of clustered DAG, identifiers for all k
% avg_c_n = Average number of links in the DAG clusters centroids

% % Save in an excel file the clusters centers for the ECLIPSE, August 21,
% 2017 = clusters_centers_eclipse
% % Export clusters centers adjacency matrices for each k to an independent excel file
% % (each k in one file, each cluster centroid in one excel sheet) 
% % By Leila Hernandez Rodriguez, 2021
% for kk=1:length(ks)
%     %filename = sprintf('cluster_centers_eclipse_k=%d_%s.xls',kk,datestr(now,'mm-dd-yyyy HH-MM'));  
%     filename = sprintf('cluster_centers_eclipse_k%d.xls',kk);  
%   for ii=1:kk
%       sheetname=sprintf('k%d_%d',kk,ii);
%       xlswrite(filename,clusters_centers_eclipse(kk).k(:,:,ii),sheetname);
%   end
% end

% Clear all outputs from using the last function
%clear preds res avg_c_n pred err n_iter is_removed centers n_dags_per_cluster intra_cluster intra_clustery clusters_centers

% %--------------------------------------------------------------------------------------------------------
% Clustering for the clear-sky day - August 18, 2017
[preds_non_eclipse,res_non_eclipse,avg_c_non_eclipse,pred_non_eclipse,err_non_eclipse,n_iter_non_eclipse,is_removed_non_eclipse,centers_non_eclipse,...
    n_dags_per_cluster_non_eclipse,intra_cluster_non_eclipse,intra_clustery_non_eclipse,clusters_centers_non_eclipse]= dag_hf_func_clusters_v2(n_dags,ks,...
    methods,n_restarts,DD,p,q,var_dag,ts_dag,max_iter,[],[],false);

% % Save in an excel file the clusters centers for August 18, 2017 = clusters_centers 
% % Export clusters centers adjacency matrices for each k to an independent excel file
% % (each k in one file, each cluster centroid in one excel sheet) 
% % By Leila Hernandez Rodriguez, 2021
% for kk=1:length(ks)
%     filename = sprintf('cluster_centers_k%d.xls',kk);  
%   for ii=1:kk
%       sheetname=sprintf('k%d_%d',kk,ii);
%       xlswrite(filename,clusters_centers(kk).k(:,:,ii),sheetname);
%   end
% end

%% ---------------------------------------------------------------------
% Deciding the best number of clusters, k 
% Dunn Index: Is the ratio of the minimum of inter-cluster distances 
% and maximum of intra-cluster distances.
% We want to maximize the Dunn index. 
% The larger the value of the Dunn index is, the better the clusters are.
% Inter-cluster distance should be the highest and the intra-cluster
% distance the smallest to have optimum clusters.
% ---------------------------------------------------------------------

% For clear-sky day, August 18, 2017
[dunn_index_non_eclipse_]=hf_dunn_index_dags(ks,n_dags,preds_non_eclipse,clusters_centers_non_eclipse,eclipse,p,q,n,m);
[dunn_value_non_eclipse, dunn_index_non_eclipse] = max(dunn_index_non_eclipse_(:,ks(keval)));

% For eclipse day - August 21, 2017
[dunn_index_eclipse_]=hf_dunn_index_dags(ks,n_dags,preds_eclipse,clusters_centers_eclipse,eclipse,p,q,n,m);
[dunn_value_eclipse, dunn_index_eclipse] = max(dunn_index_eclipse_(:,ks(keval)));

% Find optimum k using Dunn Index for the eclipse day = maximun Dunn Index 
dunn_k_to_plot = dunn_index_eclipse;

%% ---------------------------------------------------------------------
% Estimate distances between consecutive DAGs
% ---------------------------------------------------------------------
for j1 = 1:n_dags
    for j2 = j1+1:n_dags        
        pdist_eclipse(j1) = dag_dist(eclipse{j1},eclipse{j2},p,q,n,m);
        tl_eclipse(j1) = dag_dist_tl(eclipse{j1},eclipse{j2},p,q,n,m);
        pdist_non_eclipse(j1) = dag_dist(non_eclipse{j1},non_eclipse{j2},p,q,n,m);
        tl_non_eclipse(j1) = dag_dist_tl(non_eclipse{j1},non_eclipse{j2},p,q,n,m);        
    end
end

%% ---------------------------------------------------------------------
% Estimate clusters combining eclipse (e) and clear-sky (n) data
% ---------------------------------------------------------------------
% 1. All data together and make clusters. 
Dall = (vertcat(DD',D'))';           % D = eclipse data, August 21, 2017; DD = August 18, 2017 
Dall_n = length(Dall);               % Number of DAGs when appending the data from eclipse and clear-sky 

[preds_all,res_all,avg_c_n_all,pred_all,err_all,n_iter_all,is_removed_all,centers_all,n_dags_per_cluster_all,intra_cluster_all,...
    intra_clustery_all,clusters_centers_all]= dag_hf_func_clusters_v2(Dall_n,ks,methods,n_restarts,Dall,p,q,var_dag,ts_dag,max_iter,[],[],false);

% ---------------------------------------------------------------------
% 2. Find the optimal number of clusters, k
% Nominal-k method
nominal_k_all = int16((Dall_n/2)^0.5);   % For n_dags=66, nominal_k = 5.7 ~ 6

% ---------------------------------------------------------------------
% 3. Dunn Index
[dunn_index_all]=hf_dunn_index_dags(ks,Dall_n,preds_all,clusters_centers_all,Dall,p,q,n,m);

% Find optimum k using Dunn Index for all data = maximun Dunn Index 
[dunn_value_final, dunn_index_final] = max(dunn_index_all(:,ks(keval)));

% Optimum k should be larger than 2, always 
dunn_k_to_plot_final = dunn_index_final+2;
fprintf('Final Dunn Index, k = %d', dunn_k_to_plot_final);

% Select which Dunn Index to use as final k 
kevaln = find(keval==dunn_k_to_plot_final); % Using k from all data 
x_all=1:length(preds_all)';
k_to_plot = dunn_k_to_plot_final;  % This is the selected k

%k_to_plot=3;                                                              
dunn_ylim = max([dunn_value_non_eclipse dunn_value_eclipse dunn_value_final]);

%% ---------------------------------------------------------------------
% Plotting the results 
% ---------------------------------------------------------------------
 
% 1. Plot Dunn Index for clear-sky, eclipse, and data all together
%hf_dunn_index_dags_subplots(ks,dunn_ylim,dunn_index_non_eclipse_,dunn_index_eclipse_,dunn_index_all,dunn_k_to_plot_final);
hf_dunn_index_dags_subplots(ks,dunn_ylim,dunn_index_non_eclipse_,dunn_index_eclipse_,dunn_index_all,k_to_plot);

% 2. Plot clusters in time-series plots for August 21 and August 18 clusters 
x_datapoints=x_all;
y_preds=preds_all; 
color_dags = 'm';
hf_clusters_multiple_k_all(tc,D,x_datapoints,y_preds,keval,color_dags,min_n_clusters,k_to_plot,t_eclipse,sizelines)

% 3. For the optimal k, make subplots for the all clustered DAGs - August 21 and August 18
[c_eclipse_all,c_non_eclipse_all]= hf_subplots_index_dall(eclipse,non_eclipse,Dall,var_dag,ts_dag,pdist_eclipse,...
    pdist_non_eclipse,tl_eclipse,t,t_eclipse,t_total_e_start,t_total_e_end,t_part_e_start,t_part_e_end,...
    preds_all,keval,kevaln,p,q,tt1,tt2,k_to_plot);

% 4. For each cluster's centroid, find the more similar DAG created with
% data
[dist_k_vs_all_dags,labels,k_number,real_Dall_min_distance,real_Dall_min_index,real_Dall_min,real_centroid_non_eclipse_index,real_centroid_eclipse_index] = hf_kcentroids_closest_to_all_data(k_to_plot,...
centers_all,Dall,Dall_n,p,q,n,m,c_eclipse_all,preds_all,n_dags_per_cluster_all,n_dags);

%% Save the workspace with a time stemp
time = datestr(now,'mm_dd_yyyy_HH_MM');
filename = sprintf('k=%d_%s.mat',k_to_plot,time);

%% Save in an excel file the clusters centers 
% Export clusters centers adjacency matrices for each k to an independent excel file
% (each k in one file, each cluster centroid in one excel sheet) 
% By Leila Hernandez Rodriguez, 2021

d=filename; % Name of the file 
k=k_to_plot;     % Number of clusters 

filenames = sprintf('cluster_centers_all_k=%d_%s.xls',k,datestr(now,'mm-dd-yyyy_HH-MM'));  
%filenames = sprintf('without_cluster_centers_all_k_%d.xls',k);  
for i=1:k
   sheetname=sprintf('k%i',k,i);
   xlswrite(filenames,centers_all(:,:,i),sheetname);
end

%% NOTE: The centroid can be find using the index. 
%% Eclipse day: real_centroid_eclipse_index
%% Clear-sky day: real_centroid_non_eclipse_index
%% Use the excel file ``files_organization.xlsx" to find the centroid DAGs.

close all  

   

