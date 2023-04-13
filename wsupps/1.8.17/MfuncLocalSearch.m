clear all
clc
Input=dlmread('KZEdges.txt');
Nodes=max(max(Input));
[InRow, InClm]=size(Input);
for i=1:InRow
   Input(i,3)=1; % weight of the edges
end
%==============calculation the adjacency matrix==================
Adja_Mat = zeros(Nodes,Nodes); % Adjacency Matrix
for i=1:InRow
%     for j=1:InClm
        Adja_Mat(Input(i,1),Input(i,2))= Input(i,3);
        Adja_Mat(Input(i,2),Input(i,1))= Input(i,3); % if Graf is not Directed
%     end
end
           %===========Algorithm_part1============
Init_node=floor(1+rand*InRow-1);
LC(1,1)= Input(Init_node,1);
LC(1,2)= Input(Init_node,2)
k=0;
for i=1:size(LC,2)
    for j=1:Nodes
        if Adja_Mat(LC(1,i),j)==1
            k=k+1;
            N(i,k)=j; % N=hamsaye C(i,j)
        end
    end
    k=0;
end
        %==================Find NC=====================
z=0;
NC=0;
for i=1:size(N,2)
    for j=1:size(N,2)
        if N(1,i)==N(2,j)
            z=z+1;
            NC(1,z)=N(1,i);
        end
    end
end
           %==============find M measure===============
Delta=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
landa=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
% 
for i=1:size(LC,2)
    [C,D]=ismember(N(i,:),LC(1,:));
    CC(i,:)=C;
    DD(i,:)=D;
end
for i=1:size(CC,1)
    for j=1:size(CC,2)
        if CC(i,j)==1 
            Delta(LC(1,i),LC(1,DD(i,j)))=1;
            landa(LC(1,i),LC(1,DD(i,j)))=0;
        elseif N(i,j)~=0
            Delta(LC(1,i),N(i,j))=0;
            Delta(N(i,j),LC(1,i))=0;
            landa(LC(1,i),N(i,j))=1;
            landa(N(i,j),LC(1,i))=1;
        end
    end
end
           %============== the M measure==============
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
M=M_in/M_out
%
if size(NC,2)>=2
    for m=1:size(NC,2)
        LC_temp=LC;
        LC_temp(1,size(LC_temp,2)+1)=NC(1,m);
        k=0;
        for i=1:size(LC_temp,2)
            for j=1:Nodes
                if Adja_Mat(LC_temp(1,i),j)==1
                    k=k+1;
                    N_temp(i,k)=j; % N_temp=hamsaye C(i,j)
                end
            end
            k=0;
        end
%=========================find M measure============================
        Delta_temp=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
        landa_temp=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
        CC_temp=zeros(size(LC_temp,2),size(N_temp,2));
        DD_temp=zeros(size(LC_temp,2),size(N_temp,2));
        for l=1:size(LC_temp,2)
            [C,D]=ismember(N_temp(l,:),LC_temp(1,:));
            CC_temp(l,:)=C;
            DD_temp(l,:)=D;
        end
        for n=1:size(CC_temp,1)
            for o=1:size(CC_temp,2)
                if CC_temp(n,o)==1 
                    Delta_temp(LC_temp(1,n),LC_temp(1,DD_temp(n,o)))=1;
                    landa_temp(LC_temp(1,n),LC_temp(1,DD_temp(n,o)))=0;
                elseif N_temp(n,o)~=0
                    Delta_temp(LC_temp(1,n),N_temp(n,o))=0;
                    Delta_temp(N_temp(n,o),LC_temp(1,n))=0;
                    landa_temp(LC_temp(1,n),N_temp(n,o))=1;
                    landa_temp(N_temp(n,o),LC_temp(1,n))=1;
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
        M_temp(m,1)=M_in_temp/M_out_temp;
        %==========assigning node to LC===========
        if  M_temp(m,1)>M
            LC=LC_temp 
            M=M_temp(m,1) 
        end
    end
end

















