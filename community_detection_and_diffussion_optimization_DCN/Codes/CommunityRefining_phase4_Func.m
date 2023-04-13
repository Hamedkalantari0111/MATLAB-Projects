function [FinalCommunity_phase3,FinalCommunity_phase3_size,NodesExistinAllComs,Outlayer]=CommunityRefining_phase4_Func(Final_Community,Final_Community_size,Adja,Nodes)

%========================"Nodes in communities" neighbours=======================
[N_in,N_in_size]= NodeInCommFunc_Phase3(Final_Community,Final_Community_size,Adja,Nodes); % N=the neighbour of node i (node mojod dar har com)i=com, j=node in com, k= numbers of the neighbours

delta=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
landa=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
for i=1:size(Final_Community,1) %i=community
    for j=1:Final_Community_size(i,1) % j=node in final community of phase 1
            for k=1:N_in_size(i,j)
                if sum(ismember(N_in(i,k,j),Final_Community(i,:)))==1 && N_in(i,k,j)~=0
                    delta(Final_Community(i,j),N_in(i,k,j),i)=1;
                    delta(N_in(i,k,j),Final_Community(i,j),i)=1;
                    landa(Final_Community(i,j),N_in(i,k,j),i)=0; %i=community , j=is the node index in community
                else if sum(ismember(N_in(i,k,j),Final_Community(i,:)))==0 && N_in(i,k,j)~=0
                        delta(Final_Community(i,j),N_in(i,k,j),i)=0;
                        landa(Final_Community(i,j),N_in(i,k,j),i)=1;
                        landa(N_in(i,k,j),Final_Community(i,j),i)=1;
                    end
                end
            end
    end
end
%===========================calculate M measure for final community in phase 2 =====
[M_Measure]= Mmeasure_phase3Func(Final_Community,Nodes,Adja,delta,landa);
%===========mohasebe node Namojood dar community ha============
[nodes_not_in_com]= NodesnotinCommsFunc_phase3(Final_Community,Final_Community_size,Nodes); % Nodes not in the communities
%===========NC=community hamsaye moshtarak ba node u==============
[Neighbourcomun_node_notin,Neighbourcomun_node_notin_size]= NeighCommNodenotinFunc_Phase3(Final_Community,Final_Community_size,nodes_not_in_com,Adja);
%================= Final Comm of the Phase3=======================
[FinalCommunity_phase3]= FinalCommPhase3Func(Final_Community,nodes_not_in_com,Nodes,Neighbourcomun_node_notin,Neighbourcomun_node_notin_size,Adja,M_Measure); % final Community of the Phase#3
%=================================================================
k=0;
for i=1: size(FinalCommunity_phase3,1)
    FinalCommunity_phase3_size(i,1)=size(FinalCommunity_phase3,2)-sum(FinalCommunity_phase3(i,:)==0); % final Community size of the Phase#3==============
    NodesExtincom(k+1:k+FinalCommunity_phase3_size(i,1),1)=FinalCommunity_phase3(i,1:FinalCommunity_phase3_size(i,1));
    k=k+FinalCommunity_phase3_size(i,1);
end
NodesExistinAllComs=unique(NodesExtincom(:,1));
w=1;
for l=1:Nodes
    if ismember(l,NodesExistinAllComs(:,1))==0
        Outlayer(w,1)= l;
        w=w+1;
    else
        Outlayer(w,1)=0;
    end
end