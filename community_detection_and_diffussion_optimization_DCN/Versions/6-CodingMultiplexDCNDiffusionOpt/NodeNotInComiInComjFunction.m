function[Nodes_NotinCom_i_incom_j]= NodeNotInComiInComjFunction(LC,NodesOutCommunity,NodesInCommunity)
 for k=1:size(LC,1)
     for m=1:size(LC,1)
         if k~=m 
             [a b]=ismember(NodesOutCommunity(k,1:size(NodesOutCommunity,2)-sum(NodesOutCommunity(k,:)==0)),NodesInCommunity(m,1:size(NodesInCommunity,2)-sum(NodesInCommunity(m,:)==0)));
             x=1;
             for n=1:size(b,2)
                 if b(1,n)~=0
                     M(k,m,x)=NodesInCommunity(m,b(1,n));
                     x=x+1;
                 elseif b(1,size(b,2))==0
                     M(k,m,x)=0;
                 end
             end
             Nodes_NotinCom_i_incom_j(m,k)=x-1; % "Nodes_NotinCom_i_incom_j"=the number of nodes in community "k" (and not in community m!)which is related to the nodes in community "m" (and note in comm k!)
             if Nodes_NotinCom_i_incom_j(m,k)~=0
                 fi(m,k)=Nodes_NotinCom_i_incom_j(m,k)/Nodes_NotinCom_i_incom_j(m,k); % if at least one nodes in community "k" is related to at least one node in community "m"
             else
                 fi(m,k)=0; % if at least one node in community "k" is related to at least one node in community "m"
             end
         end
     end
 end