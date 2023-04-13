function[WOS,NCA,NCA_p1,NCA_p2]=WOS_NCA_Function(LC,AdjaMat,Input,InClm,LayerInputNo,F,Nodes,Alfa,Gama)
for i=1:size(LC,1)
    Commun_size(i,1)=size(LC,2)-sum(LC(i,:)==0); % size of the communities
end
[Loc_edge_index,Comm_loc_edge_index]= LocalEdgeFunc(LC,Commun_size,AdjaMat,Input,InClm,LayerInputNo,F); %Calculatinf the "edge index" in each local community

[Subscr_node,Subscript_node]= SubscriptNodeFunc(LC,Commun_size); %calculate the Subscription of the nodes % nodes eshteraki beyne har local community

[Subscr_edge,Subscript_edge]= SubscriptEdgesFunc(LC,Comm_loc_edge_index,Loc_edge_index); %calculate the Sucscription of the edges

[NeighNodesLocalCom]= NeighNodeLocalCommFunc(LC,Nodes,AdjaMat);%calculating the neighbours of the nodes in each community

[NodesInCommunity,NodesOutCommunity]= NodeInOutCommFunc(LC,Nodes,NeighNodesLocalCom); %Calculate Nodes In and Out of the Community

[Nodes_NotinCom_i_incom_j]= NodeNotInComiInComjFunction(LC,NodesOutCommunity,NodesInCommunity); % i= column and j=row; Calculatinh the "fi" parameter
 
  for i=1:size(NodesInCommunity,1)
     for j=1:size(NodesInCommunity,1)
         if i~=j 
             [c d]=ismember(NodesInCommunity(i,1:size(NodesInCommunity,2)-sum(NodesInCommunity(i,:)==0)),NodesInCommunity(j,1:size(NodesInCommunity,2)-sum(NodesInCommunity(j,:)==0))); % Nodes in comm i which are in the community j
             c_size=size(c,2);
             CCC(i,1:c_size,j)=c; %if "node in com i exist in com j"=1 , "node in com i was not exist in com j"=0
             DDD(i,1:c_size,j)=d;
         end
     end
  end  
 %============="Weighted Overlapping Score (WOS)" and finding the final "Communties"=============
 % WOS bayad Iterative bashad va nabayad yek bar mohasebe shavad
%[WOS,WOS_p1,WOS_p2]=WOSFunc(local_Community,Commun_size,Comm_loc_edge_index,Subscript_node,Alfa,Subscript_edge); % Using "WOSFunc.m" for calculation WOS
%%%========= WOS=the aboved WOS is the 1.8.17 papers WOS %%%method========
[WOS,NCA,NCA_p1,NCA_p2]=WOSFunc_New(LC,Commun_size,Comm_loc_edge_index,Subscript_node,Alfa,Subscript_edge,Gama,Nodes_NotinCom_i_incom_j,AdjaMat,NodesInCommunity,CCC); % Using "WOSFunc_part2.m" for second part of calculation WOS_NCA