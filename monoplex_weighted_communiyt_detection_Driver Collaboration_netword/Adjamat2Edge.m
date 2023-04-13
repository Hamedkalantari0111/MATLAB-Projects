function [edges,edgeCount]=Adjamat2Edge(adj)
n=length(adj); % number of nodes in L1
[ELa ELb]=find(adj>0)
edges=[ELa ELb]
% edges(:,3)=1
edges(:,3)=ELa+ELb;
edgeCount=length(edges)