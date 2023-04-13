% %=====================================================================
% %       Preparing data format for Community Detection Algorithm
% %=====================================================================
% clc
% clear all
% InputData=dlmread('edgesetsnotuniqe.txt');
% [Row,Column]=size(InputData);
% Source=InputData(:,1);
% Target=InputData(:,2);
% layer=InputData(:,3);
% %RankedWeight=InputData(:,4);
% for i=1:Row
%     InputData(i,4)=i; % InputData(i,5)= shenase edge dataset vorodi
% end
% save InputData.mat 
% %%
% clear all
% load ('InputData.mat','InputData','Row','Column')
% DNet=InputData;
% A=unique(InputData(:,1),'rows');
% B=unique(InputData(:,2),'rows');
% Temp=unique([InputData(:,1);InputData(:,2)]);
% for j=1:2
%     for k=1:size(DNet,1)
%         [row, col]=find(Temp==DNet(k,j));
%         DNet(k,j+4)=row;  %==============================================DNet(:,5) means DNet first node Index and %DNet(:,6) means DNet Second node Index
%     end
% end
% save InputNetwork.mat
clear all
load ('InputNetwork.mat','DNet','A','Temp')
NumofNodes=size(Temp,1); %%% number of the nodes envolved in the Layer 
[row, col]=find(DNet(:,5)==NumofNodes);
DNet_Senario1= DNet(1:max(row(:,1)),:);
AA=DNet_Senario1(:,1);
AA(size(DNet_Senario1,1)+1:2*size(DNet_Senario1,1))=DNet_Senario1(:,2);
[BB DD]=sort(AA);
j=1;
BB(1,2)=1;
for i=2:size(AA,1)
    if BB(i,1)>BB(i-1,1)
        j=j+1;
        BB(i,2)=j;
    else
        BB(i,2)=j;
    end
end
BB(:,3)=DD(:,1);
CC=unique(BB(:,1:2),'rows');
for i=1:size(DNet_Senario1,1)
    DNet_Senario1(i,7)= find(CC(:,1)==DNet_Senario1(i,1));
    DNet_Senario1(i,8)= find(CC(:,1)==DNet_Senario1(i,2));
end
DNetSenario1=[(DNet_Senario1(:,1:6)),(DNet_Senario1(:,7:8))]; % Final MAPPED Inputs for Driver Network in Modeling (MY MAPPED CASE STUDY DATA)!!!!!!!!!!!!!!
k=1;
for i=1:size(DNetSenario1,1)
    for j=1:size(DNetSenario1,1)
        if i~=j
            if DNetSenario1(i,7)== DNetSenario1(j,8) && DNetSenario1(i,8)==DNetSenario1(j,7) % find the Duplicate records in the Networc (because of Undircted Graph)
                FF(k,1)=i; % the index of the first duplicate record
                FF(k,2)=j; % the index of the Second duplicate record
                FF(k,3)=DNetSenario1(i,3);
                FF(k,4)=DNetSenario1(j,3);
                FF(k,5)=0.5*(FF(k,3)+FF(k,4));
                k=k+1;
            end
        end
    end
end
Duplicate=0;
w=1;
s=1;
for i=1:size(FF,1)
    if ismember(FF(i,1),Duplicate(:,1))==0
        Uniq_FF(s,1)=FF(i,1);
        Uniq_FF(s,2)=FF(i,5);
        Duplicate(w,1)=FF(i,2);
        Duplicate(w,2)=FF(i,5);
        w=w+1;
        s=s+1;
    end
end
for i=1:size(Uniq_FF,1)
    DNetSenario1(Uniq_FF(i,1),3)=Uniq_FF(i,2); 
    end
k=0;
for i=1:size(DNetSenario1,1)
    if ismember(i,Duplicate(:,1))==0
        k=k+1;
        DNetSenario2(k,:)=DNetSenario1(i,:); 
    else
    end
end
finaldataset=[DNetSenario2(:,1:4),DNetSenario2(:,7:8)];
xlswrite('FinalDataset.xls',finaldataset)
Length= size(CC,1);
save Finaldataset.mat
%%
% in finaldataset and excel file "FinalDataset" we have: 
% column 1= original source ID
% column 2= original target ID
% column 5= final extracted  source ID
% column 6= final extracted  target ID
% column 3= layer
% column 4= row counter ID