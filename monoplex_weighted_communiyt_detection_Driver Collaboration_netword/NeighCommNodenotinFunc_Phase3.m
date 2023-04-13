function[Neighbourcomun_node_notin,Neighbourcomun_node_notin_size]= NeighCommNodenotinFunc_Phase3(Final_Community,Final_Community_size,nodes_not_in_com,Adja_Mat)
k=0;
w=1;
Neighbourcomun_node_notin=zeros(size(nodes_not_in_com,1),size(nodes_not_in_com,1));
if nodes_not_in_com==0
    Neighbourcomun_node_notin=0;
else
    for i=1:size(nodes_not_in_com,1)
        for j=1:size(Final_Community,1) %j=coms
            for k=1:Final_Community_size(j,1)
                if Adja_Mat(Final_Community(j,k),nodes_not_in_com(i,1))>0 && ismember(j,Neighbourcomun_node_notin(i,:))==0
                    Neighbourcomun_node_notin(i,w)=j; % Neoghbour community of "nodes not in communtiy" ("neighbour Communities" of the "nodes_not_in")
                    w=w+1;
                end
            end
        end
        w=1;
    end
end
if Neighbourcomun_node_notin==0
    Neighbourcomun_node_notin_size=0;
else
    for i=1:size(Neighbourcomun_node_notin,1)
        Neighbourcomun_node_notin_size(i,1)=size(Neighbourcomun_node_notin,2)-sum(Neighbourcomun_node_notin(i,:)==0);
    end
end