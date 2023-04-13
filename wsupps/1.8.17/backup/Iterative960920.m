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
for w=1:100
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
        M1_in(i)=sum(MM_in(i,:));
        M1_out(i)=sum(MM_out(i,:));
    end
    M_in=0.5*sum(M1_in(:));
    M_out=0.5*sum(M1_out(:));
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
                        N_temp(i,k,w)=j; % N_temp=hamsaye C(i,j)
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
                        Loc_edge_index(i,No)=Input(l,4);
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
    Comm_loc_edge_index(i,1)=size(Loc_edge_index,2)-sum(Loc_edge_index(i,:)==0);
end
% for i=1:size(Loc_edge_index,1)
%     for j=1:size(Loc_edge_index,1)
%         if i~=j
%             A(i,j)=sum(ismember(Loc_edge_index(i,:),Loc_edge_index(j,:)))-(size(Loc_edge_index,2)-Comm_loc_edge_index(i));
%         end
%     end
% end
    %=============== calculate the Sucscription of the nodes===========
Sub=0;
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        for k=1:Commun_size(i,1)
            for l=1:Commun_size(k,1)
                if i~=j && local_Community(i,k)==local_Community(j,l)
                    Sub=Sub+1;
                end
            end
            Subscr_node(i,j,k)=Sub;
            Sub=0;
        end
        Subscript_node(i,j)=sum(Subscr_node(i,j,:));
    end
end
    %=============== calculate the Sucscription of the edges===========
Sub=0;
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        for k=1:Commun_size(i,1)
            for l=1:Commun_size(k,1)
                if i~=j && Loc_edge_index(i,k)==Loc_edge_index(j,l)
                    Sub=Sub+1;
                end
            end
            Subscr_edge(i,j,k)=Sub;
            Sub=0;
        end
        Subscript_edge(i,j)=sum(Subscr_edge(i,j,:));
    end
end
            %=============M measure=============
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
        M_p1(i,j)=Alfa*((Subscript_node(i,j))/M_p1_down);
        M_p2(i,j)=(1-Alfa)*((Subscript_edge(i,j))/M_p2_down);
        M_measure(i,j)= M_p1(i,j)+ M_p2(i,j);
    end
end
Beta=0.5;











