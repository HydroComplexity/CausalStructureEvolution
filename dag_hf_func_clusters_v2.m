function [preds,res,avg_c_n,pred,err,n_iter,is_removed,centers,n_dags_per_clusters,intra_cluster,...
    intra_clustery,clusters_centers] = dag_hf_func_clusters(n_dags,ks,methods,n_restarts,As, p, q, var_dag, ts_dag, max_iter, true_clusters,seedAs, debug)

% "dag_hf_func_clusters" function nested with the k-means fuction
% "graph_k_means"
% Estimates DAGs clusters using k-means type of algorithm using Malmi et al(2015) 

% Input:
% ks = vector of k number of clusters
% methods = vector of methods for DAG agregation and clustering
% n_restarts = Number of restarts in k-means  

% Output:
% preds = list of clustered DAG, identifiers for all k
% res = 
% avg_c_n = average number of links in the DAG clusters centroids

res = zeros(length(methods),length(ks));
preds = zeros(length(methods),n_dags,length(ks));
res_all = zeros(length(methods),n_restarts);
times = zeros(length(methods), n_restarts);
iters = zeros(length(methods), n_restarts);

% Set a consistent set of initial conditions, Leila 2021
nominal_k = int16((n_dags/2)^0.5);   % For n_dags=66, nominal_k = 5.7 ~ 6
seedy=nominal_k;
rng(seedy,'twister') % To set the seed by the user. 0 is used in Matlab per default
s = rng;
save('s')
%rng(s)   % Added by Leila, 2021
% Other options to set a random initialization
%rng('default') % To use the default 
%rng('shuffle')  % To set the seed with the date-time

tic
for ki = 1:length(ks)
    k = ks(ki);
    for mm = 1:length(methods)
        method = methods{mm};
        % Cluster
        min_err = Inf;
        min_pred = -1;
        min_n_iter = -1;
        min_active_clusters = -1;
        min_centers = -1;
        %rng(s)   % Added by Leila, 2021
        %initialization = randi(k,n_dags,1); % Random initial cluster assignments for k-means % Leila moved this here
        for j = 1:n_restarts    % n_restarts = Number of restarts in k-means
            initialization = randi(k,n_dags,1); % Random initial cluster assignments for k-means
            t_beg = tic;
            [pred,err,n_iter,is_removed,centers,n_dags_per_cluster,intra_cluster,intra_clustery] = graph_k_means_v2(...
                As, k, p, q, var_dag, ts_dag, max_iter, initialization, [], method, [], false);
            times(mm,j) = toc(t_beg);
            res_all(mm,j) = err;
            iters(mm,j) = n_iter;
            n_dags_per_clusters(ki).ki = n_dags_per_cluster; % Leila, number of dags for each cluster            
            if err < min_err
                min_err = err;
                min_pred = pred;
                min_n_iter = n_iter;
                min_active_clusters = sum(~is_removed);
                min_centers = centers;   
            end
        end
        preds(mm,:,ki) = min_pred'; % IMPORTANT! Clustered DAGs!
        avg_center_nodes = 0;
        for j = 1:k
            avg_center_nodes = avg_center_nodes + ...
                sum(sum(min_centers(:,:,j)));
        end
        avg_center_nodes = avg_center_nodes / k;
        fprintf(['%s\terr: %.1f, steps: %d, clusters: %d, ' ...
                 'avg_c_n: %.1f\n'], methods{mm}, min_err, ...
                 min_n_iter, min_active_clusters, avg_center_nodes);
        res(mm,ki) = min_err;
        avg_c_n(mm,ki) = avg_center_nodes;
        % By Leila
        clusters_centers(k).k = centers;           
   end
end
toc