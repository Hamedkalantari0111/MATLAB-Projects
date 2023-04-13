clear all
clc
Input=dlmread('KZEdges.txt');
Nodes=max(max(Input));
[InRow, InClm]=size(Input);
NN=1;
for i=1:InRow
   Input(i,3)=1; % weight of the edges
   Input(i,4)=NN; % Index of the edges
   NN=NN+1;
end
%==============calculation the adjacency matrix==================
Adja_Mat = zeros(Nodes,Nodes); % Adjacency Matrix
for i=1:InRow
    Adja_Mat(Input(i,1),Input(i,2))= Input(i,3);
    %Adja_Mat(Input(i,2),Input(i,1))= Input(i,3); % if Graf is not Directed
end
%===========Calculating the "Centrality measure" for each Node===============
s = Input(:,1);
t = Input(:,2); 
G = graph(s,t);
CC = centrality(G,'betweenness'); % could use other related measures such as 'eigenvector' or 'pagerank' 
[a.C,a.D]= sort(CC,'desc'); % Nodes sorted as Descending, D is Node name!
 %===========================Algorithm_part1=========================
            %===============Phase 1=================
%--------------% sorting edges for creating local communities--------------
x=1;
tic
% for i=1:size(a.D,1)
%     for j=1: size(a.D,1)
%         if a.D(i,1)~=a.D(j,1)
%             for k=1:size(Input,1)
%                 if Input(k,1)==a.D(i,1) && a.D(j,1)==Input(k,2)
%                     edge_Sorted(x,1)=Input(k,4); % edges sorted with respect to the aboved "Centrality measure"!!!!
%                     x=x+1;
%                 end
%             end
%         end
%     end
% end
for i=1:size(a.D,1)
    for j=1: size(a.D,1)
        if Adja_Mat(a.D(i,1),a.D(j,1))==1 %&& x<=LayerInputNo(1,1)
            edge_Sorted(x,:)=[a.D(i,1),a.D(j,1)]; % edges sorted with respect to the aboved "Centrality measure"!!!! of the whole layers
            x=x+1;
        end
    end
end
for i=1:InRow
    Adja_Mat(Input(i,2),Input(i,1))= Input(i,3); % if Graf is not Directed
end
toc
Inputtt=Input;
WW=1;
XX=1;
ZZ=1;
LC_Matrix=0;
MM=0;
for w=1:size(edge_Sorted,1)
    NC=0;
    LC=zeros(1,2);
    LC(1,:)=edge_Sorted(w,1:2); % Local Comunity
    LC_index=edge_Sorted(w,1);
    k=0;
    for i=1:size(LC,2)
        for j=1:Nodes
            if Adja_Mat(LC(1,i),j)==1
                k=k+1;
                N(i,k,w)=j; % N=Neighbour nodes of the above LC(i,j)
            end
        end
        k=0;
    end
    %==================Find NC (Neighbour Community)=====================
    z=0;
    for i=1:size(N,2)
        for j=1:size(N,2)
            if N(1,i,w)==N(2,j,w) && N(1,i,w)~=0
                z=z+1;
                NC(1,z)=N(1,i,w); %NC Common Neighbor node of two nodes (hamsaye "moshtarake" beine 2 node)
                break
            end
        end
    end
    %==============find M measure===============
    Delta=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
    landa=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
    CC=zeros(size(LC,2),size(N,2));
    DD=zeros(size(LC,2),size(N,2));
    for i=1:size(LC,2)
        [C,D]=ismember(N(i,:,w),LC(1,:));
        CC(i,:)=C; % index of Neighbour node in "N"
        DD(i,:)=D; %Neighbour Node
    end
    for i=1:size(CC,1)
        for j=1:size(CC,2)
            if CC(i,j)==1 
                Delta(LC(1,i),LC(1,DD(i,j)))=1; %if both two nodes i,j (of the edge) were in the Community =1 , O.W=0;
                landa(LC(1,i),LC(1,DD(i,j)))=0; %if only one of the nodes (of the edge) was in the Community=1 , O.W=0;
            elseif N(i,j,w)~=0
                Delta(LC(1,i),N(i,j,w))=0;
                Delta(N(i,j,w),LC(1,i))=0;
                landa(LC(1,i),N(i,j,w))=1;
                landa(N(i,j,w),LC(1,i))=1;
            end
        end
    end
    %============== the M measure============
    for i=1:Nodes
        for j=1:Nodes
            MM_in(i,j)=Adja_Mat(i,j)*Delta(i,j);
            MM_out(i,j)=Adja_Mat(i,j)*landa(i,j);
        end
        M1_in(i)=0.5*sum(MM_in(i,:));
        M1_out(i)=0.5*sum(MM_out(i,:));
    end
    M_in=0.5*sum(M1_in(:));
    M_out=1*sum(M1_out(:));
    M(w,1)=M_in/M_out;
    if size(NC,2)>2
        LC_temp=zeros(1,1);
        for m=1:size(NC,2)
            LC_temp=LC;
            LC_temp(1,size(LC_temp,2)+1)=NC(1,m);
            k=0;
            for i=1:size(LC_temp,2)
                for j=1:Nodes
                    if Adja_Mat(LC_temp(1,i),j)==1
                        k=k+1;
                        N_temp(i,k,w)=j; % N_temp=hamsaye C(i,j) hamsayehaye jadid hastand ba ezaye ezafe kardane har NC be LC_temp (for ezafi ast)
                    end
                end
                k=0;
            end
            %===============find M measure================
            Delta_temp=zeros(Nodes,Nodes); %if both two nodes i,j were in the Community =1 , O.W=0;
            landa_temp=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
            CC_temp=zeros(size(LC_temp,2),size(N_temp,2));
            DD_temp=zeros(size(LC_temp,2),size(N_temp,2));
            for l=1:size(LC_temp,2)
                [C,D]=ismember(N_temp(l,:,w),LC_temp(1,:));
                CC_temp(l,:)=C;
                DD_temp(l,:)=D;
            end
            for n=1:size(CC_temp,1)
                for o=1:size(CC_temp,2)
                    if CC_temp(n,o)==1 
                        Delta_temp(LC_temp(1,n),LC_temp(1,DD_temp(n,o)))=1;
                    elseif N_temp(n,o,w)~=0
                        landa_temp(LC_temp(1,n),N_temp(n,o,w))=1;
                        landa_temp(N_temp(n,o,w),LC_temp(1,n))=1;
                    end
                end
            end
            %==========the M measure===========
            for n=1:Nodes
                for o=1:Nodes  
                    MM_in_temp(n,o)=Adja_Mat(n,o)*Delta_temp(n,o);
                    MM_out_temp(n,o)=Adja_Mat(n,o)*landa_temp(n,o);
                end
                M1_in_temp(n)=sum(MM_in_temp(n,:));
                M1_out_temp(n)=sum(MM_out_temp(n,:));
            end
            M_in_temp=0.5*sum(M1_in_temp(:));
            M_out_temp=0.5*sum(M1_out_temp(:));
            M_temp(m,w)=M_in_temp/M_out_temp;
            %==========assigning node to LC===========
            if  M_temp(m,w)>M(w,1)
                for y=1:size(LC_temp,2)
                    local_C(WW,y)=LC_temp(1,y);
                end
                LC=LC_temp;
                M(w,1)=M_temp(m,w);
                local_M(WW,1)=M_temp(m,w);
                WW=WW+1;
            end
        end
        for yy=1:size(local_C,2)
            local_Community(XX,yy)=local_C(WW-1,yy); % local_Community: Local Communities from phase 1 extracted
        end
        XX=XX+1;
        ZZ=XX;
    end
end
save LocalCommunities_Phase1.mat % Layer 1 is the Layer of the Drivers Network of which the edge is the "frighterID"
clear all
clc
load ('LocalCommunities_Phase1.mat','local_Community','Adja_Mat','Input','Nodes')
%========================================================================
            %===============Phase 2=================
     %========calculating the ""local community edges""===========
for i=1:size(local_Community,1)
    Commun_size(i,1)=size(local_Community,2)-sum(local_Community(i,:)==0); % size of the communities
end
No=1;
for i=1:size(local_Community,1) % for each community
    for j=1:Commun_size(i,1)% for each node in local community
        for k=1:Commun_size(i,1) % be ezaye har node digar an satr
            if Adja_Mat(local_Community(i,j),local_Community(i,k))==1
                for l=1:size(Input,1)
                    if local_Community(i,j)==Input(l,1) && local_Community(i,k)==Input(l,2)
                        Loc_edge_index(i,No)=Input(l,4); % edge index in each local community
                        No=No+1;
                        break
                    end                    
                end
            end
        end
    end
    No=1;
end
for i=1:size(Loc_edge_index,1)
    Comm_loc_edge_index(i,1)=size(Loc_edge_index,2)-sum(Loc_edge_index(i,:)==0); % the number of the edges in each local community
end
    %=============== calculate the Subscription of the nodes % nodes eshteraki beyne har local community===========
Sub=0;
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        for k=1:Commun_size(i,1)
            for l=1:Commun_size(j,1)
                if i~=j && local_Community(i,k)==local_Community(j,l)
                    Sub=Sub+1;
                end
            end
            Subscr_node(i,j,k)=Sub;
            Sub=0;
        end
        Subscript_node(i,j)=sum(Subscr_node(i,j,:)); % nodes eshteraki beine har 2 local community
    end
end
    %=============== calculate the Sucscription of the edges===========
Sub=0;
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        for k=1:Comm_loc_edge_index(i,1)
            for l=1:Comm_loc_edge_index(j,1) 
                if i~=j && Loc_edge_index(i,k)==Loc_edge_index(j,l)
                    Sub=Sub+1;
                end
            end
            Subscr_edge(i,j,k)=Sub;
            Sub=0;
        end
        Subscript_edge(i,j)=sum(Subscr_edge(i,j,:)); % tedade yalhaye eshteraki beine har 2 local community
    end
end
%======================================================================
% calculating the neighbours of the nodes in each community
k=0;
 for i=1:size(local_Community,1)
     for l=1:size(local_Community,2)
         for j=1:Nodes 
             if local_Community(i,l)~=0
                 if Adja_Mat(local_Community(i,l),j)==1
                     k=k+1;
                     AA(i,k)=j; 
                 end
             end
         end
     end
     k=0;
     BB=unique(AA(i,:));
     CC=size(BB,2);
     NeighNodesLocalCom(i,1:CC)=BB; % NeighNodesLocalCom= all of the neighbour node of the community nodes
 end
 Delta_temp=zeros(Nodes,Nodes); %if both two nodes i,j were in the Community =1 , O.W=0;
 landa_temp=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
 CC_temp=zeros(size(local_Community,2),size(NeighNodesLocalCom,2));
 DD_temp=zeros(size(local_Community,2),size(NeighNodesLocalCom,2));
 for i=1:size(local_Community,1)
     [C,D]=ismember(NeighNodesLocalCom(i,:),local_Community(i,1:size(local_Community,2)-sum(local_Community(i,:)==0)));
     CC_temp(i,:)=C;
     DD_temp(i,:)=D; % the index of the local_Community where the node of the "NeighNodesLocalCom" is there  
 end
 w=0;
 y=0;
 for k=1:size(DD_temp,1)
     for l=1:size(DD_temp,2)
         if DD_temp(k,l)==0 && NeighNodesLocalCom(k,l)~=0
             w=w+1;
             NodesOutCommunity(k,w)=NeighNodesLocalCom(k,l);
         else if DD_temp(k,l)~=0 && NeighNodesLocalCom(k,l)~=0
                 y=y+1;
                 NodesInCommunity(k,y)=NeighNodesLocalCom(k,l);
             end
         end
     end
     w=0;
     y=0;
 end
 fi_base=zeros(size(NodesOutCommunity,1),size(NodesOutCommunity,2));
 for i=1:size(NodesOutCommunity,1)
     for j=1:size(NodesOutCommunity,2)
         if NodesOutCommunity(i,j)~=0
             fi_base(i,NodesOutCommunity(i,j))=1; % Outer Nodes of the community i which have a realation with the nides in community i
         end
     end
 end
 for k=1:size(local_Community,1)
     for m=1:size(local_Community,1)
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
  AAAA=0;
  for i=1:size(NodesInCommunity,1)
      for j=1:size(NodesInCommunity,1)
          if i~=j
              for k=1:size(NodesInCommunity,2)-sum(NodesInCommunity(i,:)==0)
                  for o=1:size(NodesInCommunity,2)-sum(NodesInCommunity(j,:)==0)
                      if  CCC(i,k,j)==0 || CCC(j,o,i)==0 
                          AAAA=AAAA+Adja_Mat(NodesInCommunity(i,k),NodesInCommunity(j,o));
                      end
                  end
              end
          end
          SNA_p2_up(i,j)= AAAA;
          AAAA=0;
      end
  end
   
 %============="Weighted Overlapping Score (WOS)" and finding the final "Communties"=============
 % WOS bayad Iterative bashad va nabayad yek bar mohasebe shavad
Alfa=0.5;
Gama=0.8;
%[WOS,WOS_p1,WOS_p2]=WOSFunc(local_Community,Commun_size,Comm_loc_edge_index,Subscript_node,Alfa,Subscript_edge); % Using "WOSFunc.m" for calculation WOS
%%%========= WOS=the aboved WOS is the 1.8.17 papers WOS %%%method========
[WOS,NCA,NCA_p1,NCA_p2]=WOSFunc_part2(local_Community,Commun_size,Comm_loc_edge_index,Subscript_node,Alfa,Subscript_edge,Gama,fi,Nodes_NotinCom_i_incom_j,Adja_Mat,NodesInCommunity,CCC); % Using "WOSFunc_part2.m" for second part of calculation WOS_NCA
%==================WOS=the aboved WOS is out ne WOS_NCA method============
Beta=0.5;
k_in=1;
l_in=1;
T=1;
tabu_Index=zeros(1,1);   %tabu=communityhaye ghabli ke bayad hazf shavad dar interation badi
lenght=1;
q=1;
for i=1:size(WOS,1)
    if ismember(i,tabu_Index)==0
        Community(q,l_in)=i;          %=========final Communities in phase 2 =======
        l_in=l_in+1;
    for j=1:size(WOS,2)
        if ismember(j,tabu_Index)==0
            if  j>=i && WOS(i,j)>=Beta  %============= 
                Community(q,l_in)=j;
                l_in=l_in+1;
            else
                Community(q,l_in)=0;
            end
        else
            Community(q,l_in)=0;
        end
    end
    final_Community_set_size(q,1)=size(Community,2)-sum(Community(q,:)==0);
    for m=1:final_Community_set_size(q,1)
        if Community(q,m)~=0
        tabu_Index(1,lenght)= Community(q,m);       %tabu=communityhaye ghabli ke bayad hazf shavad dar interation badi
        lenght=lenght+1;
        end
    end
    l_in=1;
    q=q+1;
    end
end
w=1;
     %=====final Community=extracting the nodes of the final community sets=========
for i=1:size(Community,1)
    for j=1:final_Community_set_size(i,1)
        for k=1:size(local_Community,2)
            if local_Community(Community(i,j),k)~=0 
                Final_Comm_set(i,w)=local_Community(Community(i,j),k);      % nodes of the final community
                w=w+1;
            end
        end
    end
    w=1;
end
% final community size (nodes in community)
for i=1:size(Final_Comm_set,1)
    Final_Comm_size(i,1)=size(Final_Comm_set,2)-sum(Final_Comm_set(i,:)==0);
end
for i=1:size(Final_Comm_set,1)
    a=unique(Final_Comm_set(i,1:Final_Comm_size(i,1)));
    for j=1:size(Final_Comm_set,2)
        if j<=size(unique(Final_Comm_set(i,1:Final_Comm_size(i,1))),2)
            aa=a(j);
        else
            aa=0;
        end
            Final_Community(i,j)= aa; % Final Comminity nodes of the second phase 
    end
end
for i=1:size(Final_Community,1)
    Final_Community_size(i,1)=size(Final_Community,2)-sum(Final_Community(i,:)==0); % the size of the Final Community
end
save LocalCommunities_Phase2.mat 
clear all
clc
load ('LocalCommunities_Phase2.mat','Final_Community','Final_Community_size','Input','Nodes','Adja_Mat')
    %======================part3: refining=========================
    %========================"Nodes in communities" neighbours=======================
k=0;
N_in=0;
for i=1:size(Final_Community,1) % i= Comunity counter
    for j=1:Final_Community_size(i,1) % j= node counter in each communtiy
        for l=1:Nodes
            if Adja_Mat(Final_Community(i,j),l)==1
                k=k+1;
                N_in(i,k,j)=l; % N=the neighbour of node i (node mojod dar har com)i=com, j=node in com, k= numbers of the neighbours
            end
        end
        k=0;
    end
    k=0;
end
  %==========size of the "neighbour nodes in community matrix"=======
for i=1:size(N_in,1)
    for j=1:size(N_in,3)
        N_in_size(i,j)=size(N_in,2)-sum(N_in(i,:,j)==0); 
    end
end
delta=zeros(1,1); %if both two nodes i,j where in the Community =1 , O.W=0;
landa=zeros(1,1); %if only one of the nodes was in the Community=1 , O.W=0;
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
%===========================calculate M measure for final community=====================
for k=1:size(Final_Community,1)
    for i=1:Nodes
        for j=1:Nodes
            MMM_in(i,j,k)=Adja_Mat(i,j)*delta(i,j,k);
            MMM_out(i,j,k)=Adja_Mat(i,j)*landa(i,j,k);
        end
        MM1_in(i,k)=sum(MMM_in(i,:,k));
        MM1_out(i,k)=sum(MMM_out(i,:,k));
    end
    MM_in(k,1)=0.5*sum(MM1_in(:,k));
    MM_out(k,1)=0.5*sum(MM1_out(:,k));
    M_Measure(k,1)=MM_in(k,1)/MM_out(k,1); % M measur of the fina communities
end

%===========mohasebe node Namojood dar community ha============
w=1;
for i=1:size(Final_Community,1)
    for j=1:Final_Community_size(i,1)
        FCN(w,1)=Final_Community(i,j); %Final_community_Nodes  not unique
        w=w+1;
    end
end
Final_community_Nodes=unique(FCN,'sorted'); %Final_community_Nodes unique

[a bb]=ismember(1:Nodes,Final_community_Nodes);
[c dd]=sort(a,'descend');
ee=sum(c(1,:)==0);
ff=size(c,2)-sum(c(1,:)==0);
for k=1:ee
    nodes_not_in_com(k,1)=dd(1,k+ff); % Nodes not in the communities
end
       %===========NC=community hamsaye moshtarak ba node u==============
k=0;
w=1;
Neighbourcomun_node_notin=zeros(size(nodes_not_in_com,1),size(nodes_not_in_com,1));
for i=1:size(nodes_not_in_com,1)
    for j=1:size(Final_Community,1) %j=coms
        for k=1:Final_Community_size(j,1)
            if Adja_Mat(Final_Community(j,k),nodes_not_in_com(i,1))==1 && ismember(j,Neighbourcomun_node_notin(i,:))==0
                Neighbourcomun_node_notin(i,w)=j; % Neoghbour community of "nodes not in communtiy" ("neighbour Communities" of the "nodes_not_in")
                w=w+1;
            end
        end
    end
    w=1;
end
for i=1:size(Neighbourcomun_node_notin,1)
    Neighbourcomun_node_notin_size(i,1)=size(Neighbourcomun_node_notin,2)-sum(Neighbourcomun_node_notin(i,:)==0);
end
%====================================================================
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
                        if Adja_Mat(Final_Community2(Neighbourcomun_node_notin(p,s),j),l)==1
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
                Final_Community(Neighbourcomun_node_notin(p,s),:)=Final_Community2(Neighbourcomun_node_notin(p,s),:);
            else
                Final_Community2(Neighbourcomun_node_notin(p,s),:)=Final_Community(Neighbourcomun_node_notin(p,s),:);
            end
        end
    end
end
FinalCommunity_phase3=Final_Community2;  %%% final Community of the Phase#3 =============================================================================
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
    end
end
save LocalCommunities_Phase3.mat 
clear all
load ('LocalCommunities_Phase3.mat','FinalCommunity_phase3','Final_Community2','M_Measure','Input','Nodes','Adja_Mat','FinalCommunity_phase3_size','NodesExistinAllComs','Outlayer')
%=========================== plloting the graph======================
for i=1:size(FinalCommunity_phase3,1)
    k=1;
    for l=1:size(FinalCommunity_phase3,1)
        for j=1:FinalCommunity_phase3_size(i,1)
            if ismember(FinalCommunity_phase3(i,j),FinalCommunity_phase3(l,:))==1 && i~=l
                OL(i,k)=FinalCommunity_phase3(i,j);
                k=k+1;
            end
        end
    end
    
end
OverlappingNodes=zeros(size(FinalCommunity_phase3,1),size(OL,2));
for i=1:size(FinalCommunity_phase3,1)
      a=unique(OL(i,:));
      b=sort(a,'desc');
      OverlappingNodes(i,1:size(b,2))=b; % Overlapping Nodes in each Community
      OverlappingNodes_Size(i,1)=size(OverlappingNodes,2)-sum(OverlappingNodes(i,:)==0); % Number of the overlapping nodes in each Community
end

s = Input(:,1); %[1 1 1 1 1 1 2 3 4 5 6 7 7 7 7 8 9 10 11 8 6];
t = Input(:,2); %[2 3 4 5 6 7 3 4 5 6 2 8 9 10 11 10 10 11 8 1 11];
G = graph(s,t);
h = plot(G)
A=FinalCommunity_phase3(1,1:FinalCommunity_phase3_size(1,1)); % Frist Community Nodes
B=FinalCommunity_phase3(2,1:FinalCommunity_phase3_size(2,1)); % Second Community Nodes
C= Outlayer(:,1);                                             % is the Outlayer Nodes
highlight(h,A,'NodeColor','g')
highlight(h,B,'NodeColor','r')
highlight(h,C,'NodeColor','y')
for j=1:size(OverlappingNodes_Size,1)
    X=OverlappingNodes(j,1:OverlappingNodes_Size(j,1));
    highlight(h,X)
end

