function  dist = dag_dist(graph1,graph2,p,q,var_dag,ts_dag,edges1,edges2)
%DAG_DIST Distance between two directed acyclic graphs (DAGs).
%   DIST = DAG_DIST(graph1,graph2,p) returns the distance between two dags
%   using the given parameter value p and q=0. Dags are given as adjacency
%   matrices.
%
%   DIST = DAG_DIST(graph1,graph2,p,q) returns the distance between two
%   dags using the given parameter values p and q. Dags are given as
%   adjacency matrices.
%
%   DIST = DAG_DIST(graph1,graph2,p,q,edges1,edges2) returns the distance
%   between two dags given as edge lists (M-by-2 matrices). For sparse
%   graphs this option is more efficient.

% Copyright (c) 2015 Eric Malmi

% Updated by Leila Hernandez, 2021
% Estimate the number of edges in a fully connected DAG based on n and m
% for DAGs with no contemporary edges
% n = Number of variables in the DAG
% m = Number of time steps in the DAG
% e = Number of edges in the DAG

if nargin == 3
    q = 0;
end
if nargin <= 6   % Originally was 4. Leila added var_dag,ts_dag as inputs.
    temp = graph1.*graph2';
    % Discordant pairs
    D = sum(temp(:));
    deltaE = sum(abs(graph1(:) - graph2(:)));
    temp2 = graph1.*graph2;
    % Concordant pairs
    C = sum(temp2(:));
    %nn = size(graph1,1);
    %tl = sum(graph1(:) == 1);  % By Leila 2021, counts the number of edges in the first DAG
    e = var_dag*(ts_dag*(ts_dag-1)/2)+var_dag*(var_dag-1)*(ts_dag*(ts_dag-1)/2); % By Leila 2021, number of edges for a fully connected DAG with no contenporaneous edges
    %e = 1375;
    dist = D + p*(deltaE-2*D) + q*(e -(deltaE-D+C)); % By Leila 2021, based on Malmi at al 2015
    %fprintf('%d %d\n', e, deltaE-D+C);
    %dist = p*(deltaE-2*D) + D + q*(n*(n-1)/2-(deltaE-D+C)); Malmi 2015
    %fprintf('%d %d %d %d %d\n', C, D, (deltaE-2*D), e, (e-(deltaE-D+C)));    
else % Edge lists are available
    n1 = size(edges1,1);
    n2 = size(edges2,1);
    if n1 < n2
        edges = edges1;
        graph = graph2;
    else
        edges = edges2;
        graph = graph1;
    end
    case1 = 0; % Concordant pairs
    case2 = 0; % Discordant pairs
    for i = 1:size(edges,1)
        if graph(edges(i,1),edges(i,2))
            case1 = case1 + 1;
        elseif graph(edges(i,2),edges(i,1))
            case2 = case2 + 1;
        end
    end
    case3 = n1-(case1+case2) + n2-(case1+case2);
    nv = size(graph1,1);
    case4 = (nv*nv-nv)/2 - case3 - case2 - case1; % Semi-ambiguous pairs
    dist = case2 + p*case3 + q*case4; % Ambiguous pairs
    %fprintf('%d %d %d %d\n', case1,case2,case3,case4);
end
