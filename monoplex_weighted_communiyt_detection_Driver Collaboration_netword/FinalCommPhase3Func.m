function[FinalCommunity_phase3]= FinalCommPhase3Func(Final_Community,nodes_not_in_com,Nodes,Neighbourcomun_node_notin,Neighbourcomun_node_notin_size,Adja_Mat,M_Measure)
PreFinal_Community=Final_Community;
Final_Community2=Final_Community;
% Final_Community_size2=Final_Community_size;
y=1;
for p=1:size(nodes_not_in_com,1)
    for z=1:size(Final_Community2,1) % counts of the communities from phase 2
        Final_Community2_size(z,1)=size(Final_Community2,2)-sum(Final_Community2(z,:)==0); % calculate the community size (final comm from phase 2)
    end
    if Neighbourcomun_node_notin_size(p,1)~=0 % sometimes "Neighbourcomun_node_notin_size" could be equale to zero!!!!!
        for s=1:Neighbourcomun_node_notin_size(p,1) % s= tedade community moshtarak baraye har node namojod dar communityha
            FF2=nodes_not_in_com(p,1);
            Final_Community2(Neighbourcomun_node_notin(p,s),Final_Community2_size(Neighbourcomun_node_notin(p,s),1)+1)=FF2;
            
            Final_Community2_size(Neighbourcomun_node_notin(p,s),1)=size(Final_Community2,2)-sum(Final_Community2(Neighbourcomun_node_notin(p,s),:)==0);
   %============"assigning Node in communtiy" neighbours=========
            k=0;
            N_in2=0;
%             for i=1:size(Final_Community2,1) %tabeiee az p va s,dar j hast
                for j=1:Final_Community2_size(Neighbourcomun_node_notin(p,s),1)
                    for l=1:Nodes
                        if Adja_Mat(Final_Community2(Neighbourcomun_node_notin(p,s),j),l)>0
                            k=k+1;
                            N_in2(j,k)=l; % N=hamsaye node j in community x of Final_Community2(node in coms)
                        end
                    end
                    k=0;
                end
                k=0;
%             end 
    %==========size of the "neighbour nodes in communtiy matrix"=======
            for i=1:size(N_in2,1)
                N_in2_size(i,1)=size(N_in2,2)-sum(N_in2(i,:)==0); % the size of the N_in2
            end 
            delta2=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
            landa2=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
            for j=1:Final_Community2_size(Neighbourcomun_node_notin(p,s),1) % j=node in community
                for k=1:N_in2_size(j,1)
                    if sum(ismember(N_in2(j,k),Final_Community2(Neighbourcomun_node_notin(p,s),:)))==1 && N_in2(j,k)~=0
                        delta2(Final_Community2(Neighbourcomun_node_notin(p,s),j),N_in2(j,k))=1;
                        delta2(N_in2(j,k),Final_Community2(Neighbourcomun_node_notin(p,s),j))=1;
                        landa2(Final_Community2(Neighbourcomun_node_notin(p,s),j),N_in2(j,k))=0; %i=community , j=is the node index in community
                    else if sum(ismember(N_in2(j,k),Final_Community2(Neighbourcomun_node_notin(p,s),:)))==0 && N_in2(j,k)~=0
                            delta2(Final_Community2(Neighbourcomun_node_notin(p,s),j),N_in2(j,k))=0;
                            landa2(Final_Community2(Neighbourcomun_node_notin(p,s),j),N_in2(j,k))=1;
                            landa2(N_in2(j,k),Final_Community2(Neighbourcomun_node_notin(p,s),j))=1;
                        end
                    end
                end
            end
            %=======calculate M measure for final community===============
            for i=1:Nodes
                for j=1:Nodes
                    MMMM_in(i,j)=Adja_Mat(i,j)*delta2(i,j);
                    MMMM_out(i,j)=Adja_Mat(i,j)*landa2(i,j);
                end
                MMM1_in(i,1)=sum(MMMM_in(i,:));
                MMM1_out(i,1)=sum(MMMM_out(i,:));
            end
            MMM_in=1*sum(MMM1_in(:,1));
            MMM_out=1*sum(MMM1_out(:,1));
            M_Measure_new(Neighbourcomun_node_notin(p,s),1)=(0.5*MMM_in)/(0.5*MMM_out); % M measur of the fina communities
            if M_Measure_new(Neighbourcomun_node_notin(p,s),1) > M_Measure(Neighbourcomun_node_notin(p,s),1)
                M_Measure(Neighbourcomun_node_notin(p,s),1)=M_Measure_new(Neighbourcomun_node_notin(p,s),1);
                Final_Community(Neighbourcomun_node_notin(p,s),1:size(Final_Community2,2)-sum(Final_Community2(Neighbourcomun_node_notin(p,s),:)==0))= Final_Community2(Neighbourcomun_node_notin(p,s),1:size(Final_Community2,2)-sum(Final_Community2(Neighbourcomun_node_notin(p,s),:)==0));
            else
                Final_Community2(Neighbourcomun_node_notin(p,s),1:size(Final_Community,2)-sum(Final_Community(Neighbourcomun_node_notin(p,s),:)==0)) = Final_Community (Neighbourcomun_node_notin(p,s),1:size(Final_Community,2)-sum(Final_Community(Neighbourcomun_node_notin(p,s),:)==0));
            end
        end
    end
end
FinalCommunity_phase3=Final_Community2;  %%% final Community of the Phase#3 =============================================================================