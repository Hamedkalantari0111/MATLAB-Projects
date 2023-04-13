clear all
clc

input= xlsread('Input1.xlsx');
[Rows,Clmn] =size(input);

input1=zeros(Rows,Clmn+2);
Sources(:,1)=input(:,1);
Targets(:,1)=input(:,2);
input1(:,1:3)=input(:,1:3);


All=zeros(2*Rows,1);
All(1:Rows,1)=Sources(:,1);
All(Rows+1:2*Rows,1)=Targets(:,1);
AllNodes=unique(All);

AdjacencyMatrix= zeros(size(AllNodes,1),size(AllNodes,1));
for i=1:size(input,1);
    k=input(i,1);
    m=input(i,2);
    AdjacencyMatrix (k,m)=1;
end

l=1;
n=1;
Average=zeros(size(AllNodes,1),1);
for i=1:size(AllNodes,1)
    Average(i,1)=(sum(AdjacencyMatrix(i,:))+sum(AdjacencyMatrix(:,i)))/2;
    if Average(i,1)<1;
        LeavesList(l,1)=i;
        l=l+1;
    elseif Average(i,1)>=1;
        NonLeavesList(n,1)=i;
        n=n+1;
    end
end

D=zeros(size(AdjacencyMatrix,1),3); % D(:,1)= node ID , D(:,2)= node distance , D(:,3)= node weight

for i=1:size(AdjacencyMatrix,1);
    for j=1:size(AdjacencyMatrix,1);
        if i==1 && j==1;
            D(j,2)=0;
            D(j,3)=1;
            D(j,1)=j;
        end
        if AdjacencyMatrix(i,j)==1;
            if D(j,2)==0;
                D(j,2)=D(i,2)+1;
                D(j,3)=D(i,3);
                D(j,1)=j;
            elseif D(j,2)>0 && D(j,2)==D(i,2)+1;
                D(j,3)=D(j,3)+D(i,3);
                D(j,1)=j;
            elseif D(j,2)>0 && D(j,2)<D(i,2)+1;
                D(j,2)=D(j,2);
                D(j,3)=D(j,3);
                D(j,1)=j;
            end
        end
    end
end

for o=1:size(input1,1);
    input1(o,4)= D(input1(o,1),3)/D(input1(o,2),3);
end

[D_Reverced , D_Reverced_Id]=sort(D(:,2),'descend');
input2=zeros(size(input1,1),1);
for o=1:size(D_Reverced_Id,1);
    for p=1:size(input1,1);
        if input1(p,2)==D_Reverced_Id(o,1);
            if ismember(LeavesList,D_Reverced_Id(o,1))==1;
                input1(p,5)= input1(p,4);
            elseif ismember(LeavesList,D_Reverced_Id(o,1))==0;
                for l=1:size(input1,1);
                    if input1(l,1)==D_Reverced_Id(o,1);
                        input2(p,1)= input2(p,1)+ input1(l,5);
                    else
                        input2(p,1)= input2(p,1);
                    end
                end
                 input1(p,5)=input1(p,4)*(1+ input2(p,1));
            end
        end
    end
end

Score= zeros(1,size(input1,1));
for i=1:size(input1,1);
    Score(1,i)=input1(i,5);
end
[ScoreValue ScoreID]= sort(Score(1,:));
value=ScoreValue(1,size(ScoreValue,2));
ID=ScoreID(1,size(ScoreID,2));















