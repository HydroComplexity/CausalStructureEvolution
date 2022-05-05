function [dunn_index]=hf_dunn_index_dags(ks,n_dags,y_intra,clusters_centers,data,p,q,n,m)

%% Function to calculate the Dunn Index for clusters of DAGs
% Dunn Index: The ratio of the minimum of inter-cluster distances and the
% maximum of intra-cluster distance.

% By Leila Hernandez Rodriguez, 2021

clear kcenter_dist intra_cluster

kcenter_dist = [];
intra_cluster = [];

% Intra-cluster distance 
for i=1:length(ks)
    preds_e_k = y_intra(:,1:n_dags,i)';  
    for j = 1:length(preds_e_k)
        kcenter_e = clusters_centers(i).k(:,:,preds_e_k(j)); 
        kcenter_dist_e(j) = dag_dist(data{j},kcenter_e,p,q,n,m);         
    end  
    intra_cluster_max_e(i)=max(kcenter_dist_e);   % Intra-cluster distance        
end

% Inter-cluster distance 
for k=2:length(ks)
    for i=1:k
        for j=1:k
            if i~=j
               kcenter_uno_e=clusters_centers(k).k(:,:,i); 
               kcenter_dos_e=clusters_centers(k).k(:,:,j);
               kinter_dist_e(i) = dag_dist(kcenter_uno_e,kcenter_dos_e,p,q,n,m);         
            end
        end
    end
    inter_cluster_min_e(k)=min(kinter_dist_e);
end

dunn_index=inter_cluster_min_e./intra_cluster_max_e;

