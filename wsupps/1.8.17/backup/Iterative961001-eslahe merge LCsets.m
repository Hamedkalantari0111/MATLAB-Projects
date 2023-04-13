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
    Adja_Mat(Input(i,2),Input(i,1))= Input(i,3); % if Graf is not Directed
end
           %===========Algorithm_part1============
%========================================================================
            %===============Phase 1=================
%========================================================================
Inputtt=Input;
WW=1;
XX=1;
ZZ=1;
LC_Matrix=0;
MM=0;
for w=1:1000
    Init_node=floor(1+rand*(size(Inputtt,1)-1));
    NC=0;
    LC=zeros(1,2);
    LC(1,1)= Inputtt(Init_node,1);
    LC(1,2)= Inputtt(Init_node,2);
    LC_index=Init_node;
    k=0;
    for i=1:size(LC,2)
        for j=1:Nodes
            if Adja_Mat(LC(1,i),j)==1
                k=k+1;
                N(i,k,w)=j; % N=hamsaye C(i,j)
            end
        end
        k=0;
    end
    %==================Find NC=====================
    z=0;
    for i=1:size(N,2)
        for j=1:size(N,2)
            if N(1,i,w)==N(2,j,w) && N(1,i,w)~=0
                z=z+1;
                NC(1,z)=N(1,i,w); 
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
        CC(i,:)=C;
        DD(i,:)=D;
    end
    for i=1:size(CC,1)
        for j=1:size(CC,2)
            if CC(i,j)==1 
                Delta(LC(1,i),LC(1,DD(i,j)))=1;
                landa(LC(1,i),LC(1,DD(i,j)))=0;
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
                        N_temp(i,k,w)=j; % N_temp=hamsaye C(i,j) hamsayehaye jadid hastand be ezaye ezafe kardane har NC be LC_temp
                    end
                end
                k=0;
            end
            %===============find M measure================
            Delta_temp=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
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
                if XX==ZZ
                    for u=1:size(LC,2)
                        LC_Matrix(XX,u+2)=LC(1,u);
                    end
                    LC_Matrix(XX,2)=LC_index; % edge Number of the input
                    LC_Matrix(XX,1)=XX; % rows of the LC-matrix
                    f=1;
                    for e=1:size(Inputtt,1)
                        if ismember(Inputtt(e,4),LC_Matrix(:,2))==0
                            Inputt(f,1:3)=Inputtt(e,1:3);
                            Inputt(f,4)=f;
                            f=f+1;
                        end
                    end
                end
                ZZ=ZZ+1;
                LC=LC_temp;
                M(w,1)=M_temp(m,w);
                local_M(WW,1)=M_temp(m,w);
                WW=WW+1;
            end
        end
        Inputtt=zeros(1,1);
        Inputtt=Inputt;
        Inputt=zeros(1,1);
        for yy=1:size(local_C,2)
            local_Community(XX,yy)=local_C(WW-1,yy);
        end
        XX=XX+1;
        ZZ=XX;
    end
end
%========================================================================
            %===============Phase 2=================
%========================================================================
     %========calculating the ""local community edges""===========
for i=1:size(local_Community,1)
    Commun_size(i,1)=size(local_Community,2)-sum(local_Community(i,:)==0);
end
No=1;
for i=1:size(local_Community,1)
    for j=1:Commun_size(i,1)
        for k=1:Commun_size(i,1) 
            if Adja_Mat(local_Community(i,j),local_Community(i,k))==1
                for l=1:size(Input,1)
                    if local_Community(i,j)==Input(l,1) && local_Community(i,k)==Input(l,2)
                        Loc_edge_index(i,No)=Input(l,4); % edges of each local community
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
    %=============== calculate the Subscription of the nodes===========
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
        Subscript_node(i,j)=sum(Subscr_node(i,j,:)); % nodes eshteraki beyne har local community
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
        Subscript_edge(i,j)=sum(Subscr_edge(i,j,:)); % yalhaye eshteraki beyne har local community
    end
end
            %============="Weighted Overlapping Score" and finding the final "Communties"=============
Alfa=0.5;
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        if Commun_size(i,1)<=Commun_size(j,1)
            M_p1_down=Commun_size(i,1);
        else
            M_p1_down=Commun_size(j,1);
        end
        if  Comm_loc_edge_index(i,1)<= Comm_loc_edge_index(j,1)
            M_p2_down= Comm_loc_edge_index(i,1);
        else
            M_p2_down= Comm_loc_edge_index(j,1);
        end
        WOS_p1(i,j)=Alfa*((Subscript_node(i,j))/M_p1_down);
        WOS_p2(i,j)=(1-Alfa)*((Subscript_edge(i,j))/M_p2_down);
        WOS(i,j)= WOS_p1(i,j)+ WOS_p2(i,j); % Weighted Overlapping Score
    end
end
Beta=0.5;
k_in=1;
l_in=1;
T=1;
tabu_Index=zeros(1,1);
lenght=1;
q=1;
for i=1:size(WOS,1)
    if ismember(i,tabu_Index)==0
        Community(q,l_in)=i;          %=========final Communities=======
        l_in=l_in+1;
    for j=1:size(WOS,2)
        if ismember(j,tabu_Index)==0
            if  j>=i && WOS(i,j)>=Beta  
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
                Final_Comm_set(i,w)=local_Community(Community(i,j),k);
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
            Final_Community(i,j)= aa; % Final Comminity
    end
end
for i=1:size(Final_Community,1)
    Final_Community_size(i,1)=size(Final_Community,2)-sum(Final_Community(i,:)==0); % the size of the Fival Community
end

    %======================part3: refining=========================
    %========================"Nodes in coms" neighbours=======================
k=0;
N_in=0;
for i=1:size(Final_Community,1)
    for j=1:Final_Community_size(i,1)
        for l=1:Nodes
            if Adja_Mat(Final_Community(i,j),l)==1
                k=k+1;
                N_in(i,k,j)=l; % N=hamsaye node i (node in coms)
            end
        end
        k=0;
    end
    k=0;
end
  %==========size of the "neighbour nodes in com matrix"=======
for i=1:size(N_in,1)
    for j=1:size(N_in,3)
        N_in_size(i,j)=size(N_in,2)-sum(N_in(i,:,j)==0); 
    end
end
delta=zeros(1,1); %if both two nodes i,j where in the Community =1 , O.W=0;
landa=zeros(1,1); %if only one of the nodes was in the Community=1 , O.W=0;
for i=1:size(Final_Community,1) %i=community
    for j=1:Final_Community_size(i,1) % j=node in community
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
Final_community_Nodes=unique(FCN,'sorted');

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
    for j=1:size(Final_Community,1)
        for k=1:Final_Community_size(j,1)
            if Adja_Mat(Final_Community(j,k),nodes_not_in_com(i,1))==1 && ismember(j,Neighbourcomun_node_notin(i,:))==0
                Neighbourcomun_node_notin(i,w)=j; % Neoghbour community of "nodes not incommuntiy"
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
    for z=1:size(Final_Community2,1)
        Final_Community2_size(z,1)=size(Final_Community2,2)-sum(Final_Community2(z,:)==0);
    end
    if Neighbourcomun_node_notin_size(p,1)~=0 
        for s=1:Neighbourcomun_node_notin_size(p,1) % s= tedade community moshtarak baraye har node namojod dar communityha
            FF2=nodes_not_in_com(p,1);
            Final_Community2(Neighbourcomun_node_notin(p,s),Final_Community2_size(Neighbourcomun_node_notin(p,s),1)+1)=FF2;
            
            Final_Community2_size(Neighbourcomun_node_notin(p,s),1)=size(Final_Community2,2)-sum(Final_Community2(Neighbourcomun_node_notin(p,s),:)==0);
   %============"assigning Node in communtiy" neighbours=========
            k=0;
            N_in2=0;
%             for i=1:size(Final_Community2,1)
                for j=1:Final_Community2_size(Neighbourcomun_node_notin(p,s),1)
                    for l=1:Nodes
                        if Adja_Mat(Final_Community2(Neighbourcomun_node_notin(p,s),j),l)==1
                            k=k+1;
                            N_in2(j,k)=l; % N=hamsaye node i (node in coms)
                        end
                    end
                    k=0;
                end
                k=0;
%             end
    %==========size of the "neighbour nodes in com matrix"=======
            for i=1:size(N_in2,1)
                N_in2_size(i,1)=size(N_in2,2)-sum(N_in2(i,:)==0);
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
%=========================== M measure======================

