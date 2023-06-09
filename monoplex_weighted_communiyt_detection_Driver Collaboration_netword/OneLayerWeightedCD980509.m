% clear all
% clc
% % ================= Preparing ID of source and targets===========
% load ('Finaldataset.mat','finaldataset','NumofNodes')
% % in finaldataset and excel file "FinalDataset" we have: 
% % column 1= original source ID
% % column 2= original target ID
% % column 8= final extracted  source ID
% % column 9= final extracted  target ID
% % column 3= weight value
% % column 4= weight rank value in which "1" is low, "2" is medium and "3" is high
% % column 5= row counter ID
% AllData=[finaldataset(:,6);finaldataset(:,7)];
% DataUniq(:,2)=unique(AllData(:,1),'rows');
% DataUniq(1:size(DataUniq,1),1)=1:size(DataUniq,1);
% for i=1:size(finaldataset,1)
%     finaldataset(i,8)=DataUniq(find(DataUniq(:,2)==finaldataset(i,6)),1);
%     finaldataset(i,9)=DataUniq(find(DataUniq(:,2)==finaldataset(i,7)),1);
% end
% finaldataset(1:size(finaldataset,1),10)=1; % is the layer one lable
% %%=================== Creatind the CaseStusy Input===============
% Input=[finaldataset(:,8:10),finaldataset(:,3),finaldataset(:,5)]; %colm1=sourc, colmn2=target, colmn3=layer number, colmn4=weight, colmn5=ID
% %===================Examples================================
% % Input=dlmread('C:\Users\kalantari.h\Dropbox\Social Network\1- Paper WorkGroup\1- Monoplex CommunityDetectionandDiffusion\1-Codes\1-MonoplexDCN980429\datasets\KZOneLayer.txt');
% Nodes=NumofNodes;
% [InRow, InClm]=size(Input);
% layerNo=Input(InRow,3); % Number of the layers
% NN=1;
% for i=1:InRow
%     if Input(i,3)==1
%         Input(i,InClm+1)=NN; %colmn6= new ID
%         NN=NN+1;
%     end
% end
% LayerInputNo(1,1)=NN-1;
% %==============calculation the adjacency matrix Layer1=================
% %Adja_Mat = zeros(Nodes,Nodes); % Adjacency Matrix of the "Layer #1"
% for i=1:LayerInputNo(1,1)
%     Adja_Mat(Input(i,1),Input(i,2))= Input(i,4);  % Adjacency Matrix of the "Layer #1" and 4 is the weight column number!!!!
%     s.L1(i,1)=Input(i,1);
%     t.L1(i,1)=Input(i,2);
%     w.L1(i,1)=Input(i,4);
% end
% %%===========Calculating the "Centrality measure" for each Node============
% G = graph(s.L1,t.L1,w.L1);
% % % Adja_Mat=adjacency(G,'weighted');
% % % A=adjacency(G);
% %CC.unweighted = centrality(G,'betweenness'); % unweighted centrality measure
% CC.weighted = centrality(G,'betweenness','Cost',G.Edges.Weight); % weighted centrality measure and could use other related measures such as 'eigenvector' or 'pagerank' 
% [a.C,a.D]= sort(CC.weighted,'desc'); % Nodes sorted as Descending, D is Node name!
% for i=1:LayerInputNo(1,1)
%     Adja_Mat(Input(i,2),Input(i,1))= Input(i,4); % if Graf is not Directed
% end
% save preparing_Dataset.mat 
% %%
% %%===============Phase 0: "Edge probing"===========
% clear all
% clc
% load ('preparing_Dataset.mat','Adja_Mat','a')
% %--------------% sorting edges for creating local communities--------------
% e_phase0=cputime;
% a1=a.D;
% e_phase0_1=cputime;
% [edge_Sorted,NodeEdgesL1] = NodeProbing_Func(a1,Adja_Mat); % Finding the Sorted Edges in Layer 1
% Time.Phase0=cputime-e_phase0_1;
% save LocalCommunities_Phase0_CaseStudy.mat 
% %%
% %=================== Pahse 1: Local Community Detection===========
% clear all
% clc
% load ('LocalCommunities_Phase0_CaseStudy.mat','edge_Sorted','Time')
% load ('preparing_Dataset.mat','Nodes','Adja_Mat','LayerInputNo','Input')
% Resolution=4; %the greater the Resolution the lower the comms
% e_phase1=cputime;
% [local_Community,local_Community_size]=LocalCommunity_Phase1_Func(edge_Sorted,Adja_Mat,Nodes,Resolution); % Calculate the Lcoal ommunity of the Layer #1!!!!!!!!!!!!!!!!!!!!!!
% Time.Phase1=cputime-e_phase1;
% L1_Input=Input(1:LayerInputNo(1,1),:);
% [EQ_phase1]=EQ_Measure_Func(L1_Input,Nodes,local_Community,local_Community_size,Adja_Mat);% calculae the Q measure for the Local communities in Layer#1
% EQ.Phase1=EQ_phase1;
% save LocalCommunities_Phase1_CaseStudy.mat 
% %%==============finding the overlapping nodes proportion in comms=====
% clear all
% clc
% load ('LocalCommunities_Phase1_CaseStudy.mat','local_Community','local_Community_size')
% load ('preparing_Dataset.mat','InClm','Nodes','Adja_Mat','LayerInputNo','Input')
% %======================finding outliers in comms=====================
% [Outlayer_L1,NodesExtincom_L1]=OutlierDetectionFunc(local_Community,local_Community_size,Nodes);
% %===========finding the overlapping nodes=================
% [OverlappingNodes_L1,OverlappingNodes_Size_L1,OverlappingNodes_proportion_L1,proportion_L1,comm_count_node_L1,comm_proportion_node_L1]=overlappingDetectionFunc(local_Community,local_Community_size,NodesExtincom_L1);
% save LocalCommunities_Phase1_CaseStudy_Overlaps_infs.mat 
% %%
% %===============Phase 2: Intra-Layer Community Merging===========================
clear all
clc
load ('LocalCommunities_Phase1_CaseStudy.mat','local_Community','Input','LayerInputNo','EQ','Time')
load ('preparing_Dataset.mat','Adja_Mat','Nodes','InClm')
% Final_Community=local_Community;
% Final_Community_size=size(local_Community,2);
%FF=zeros(10,size(local_Community,1),Nodes);
m=1;
for i=1:size(local_Community,1)
    AA(i,1)=size(local_Community,2)-sum(local_Community(i,:)==0);
    if AA(i,1)<Nodes
        CC(m,1:size(local_Community,2)-sum(local_Community(i,:)==0))=local_Community(i,1:size(local_Community,2)-sum(local_Community(i,:)==0));
        m=m+1;
    end
end
F=0;
Alfa=0.8;
Gama=0.85;
Beta=0.8;
%%================================================================
e_phase2=cputime;
[Final_Community,Final_Community_size]=CommunityMerging_Phase2_Func(CC,Adja_Mat,Nodes,LayerInputNo(1,1),Input,InClm,F,Alfa,Gama,Beta); 
Time.Phase2=cputime-e_phase2;
L1_Input=Input(1:LayerInputNo(1,1),:);
[EQ_phase2]=EQ_Measure_Func(L1_Input,Nodes,Final_Community,Final_Community_size,Adja_Mat);% calculae the Q measure for the communities in Layer#1
EQ.Phase2=EQ_phase2;
save LocalCommunities_Phase2_CaseStudy.mat 
%%
%%==============finding the overlapping nodes proportion in comms=====
clear all
clc
load ('LocalCommunities_Phase2_CaseStudy.mat','Final_Community','Final_Community_size')
load ('preparing_Dataset.mat','InClm','Nodes','Adja_Mat','LayerInputNo','Input')
%======================finding outliers in comms=====================
[Outlayer_L1,NodesExtincom_L1]=OutlierDetectionFunc(Final_Community,Final_Community_size,Nodes);
%===========finding the overlapping nodes=================
[OverlappingNodes_L1,OverlappingNodes_Size_L1,OverlappingNodes_proportion_L1,proportion_L1,comm_count_node_L1,comm_proportion_node_L1]=overlappingDetectionFunc(Final_Community,Final_Community_size,NodesExtincom_L1);
save LocalCommunities_Phase2_CaseStudy_Overlaps_infs.mat 
clc
%%
% % %========================Phase 3: Refining communities ================================
clear all
clc
load ('LocalCommunities_Phase2_CaseStudy.mat','Final_Community','Final_Community_size','EQ','Time')
load ('preparing_Dataset.mat','Adja_Mat','Nodes','LayerInputNo','Input')
e_phase3=cputime;
[FinalCommunity_phase3,FinalCommunity_size_phase3,NodesExistinAllComs,Outlayer]=CommunityRefining_phase4_Func(Final_Community,Final_Community_size,Adja_Mat,Nodes); % Final Comm Phase 3 for layer1
Time.Phase3=cputime-e_phase3;
%-------------------Calculating the EQ measure-------------------
L1_Input=Input(1:LayerInputNo(1,1),:);
[EQ_phase3]=EQ_Measure_Func(L1_Input,Nodes,FinalCommunity_phase3,FinalCommunity_size_phase3,Adja_Mat);
EQ.Phase3=EQ_phase3;
save LocalCommunities_Phase3_CaseStudy.mat 
%%
%%==============finding the overlapping nodes proportion in comms=====
clear all
clc
load ('LocalCommunities_Phase3_CaseStudy.mat','FinalCommunity_phase3','FinalCommunity_size_phase3')
load ('preparing_Dataset.mat','InClm','Nodes','Adja_Mat','LayerInputNo','Input')
%======================finding outliers in comms=====================
[Outlayer_L1,NodesExtincom_L1]=OutlierDetectionFunc(FinalCommunity_phase3,FinalCommunity_size_phase3,Nodes);
%===========finding the overlapping nodes=================
[OverlappingNodes_L1,OverlappingNodes_Size_L1,OverlappingNodes_proportion_L1,proportion_L1,comm_count_node_L1,comm_proportion_node_L1]=overlappingDetectionFunc(FinalCommunity_phase3,FinalCommunity_size_phase3,NodesExtincom_L1);
save LocalCommunities_Phase3_CaseStudy_Overlaps_infs.mat 
clc
%%
% %======================        resutls       ===========================
load ('LocalCommunities_Phase0_CaseStudy.mat')
load ('LocalCommunities_Phase1_CaseStudy.mat')
load ('LocalCommunities_Phase2_CaseStudy.mat')
load ('LocalCommunities_Phase3_CaseStudy.mat','EQ','Time')
disp('EQ measure for phase One is:'),EQ.Phase1
disp('EQ measure for phase Two is:'),EQ.Phase2
disp('EQ measure for phase Three is:'),EQ.Phase3
Time.TotalCpuTime=Time.Phase0+Time.Phase1+Time.Phase2+Time.Phase3
save Results_CaseStudy.mat
%%
% %=========================== ploting the graph=============================
load ('preparing_Dataset.mat','InClm','Nodes','Adja_Mat','LayerInputNo','Input')
load ('LocalCommunities_Phase3_CaseStudy.mat','FinalCommunity_phase3','FinalCommunity_size_phase3')
for i=1:size(FinalCommunity_phase3,1)
    k=1;
    for l=1:size(FinalCommunity_phase3,1)
        for j=1:FinalCommunity_size_phase3(i,1)
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

s = Input(:,1); 
t = Input(:,2); 
G = graph(s,t);
D = degree(G); % Degree of the nodes!
h = plot(G);
aa=1;
bb=2;
A=FinalCommunity_phase3(aa,1:FinalCommunity_size_phase3(aa,1)); % Frist Community Nodes
B=FinalCommunity_phase3(bb,1:FinalCommunity_size_phase3(bb,1)); % Second Community Nodes
%C= Outlayer(:,1);                                             % is the Outlayer Nodes
highlight(h,A,'NodeColor','g')
highlight(h,B,'NodeColor','r')
%highlight(h,C,'NodeColor','y')
for j=1:size(OverlappingNodes_Size,1)
    X=OverlappingNodes(j,1:OverlappingNodes_Size(j,1));
    highlight(h,X)
end
save LocalCommunities_Phase4_CaseStudy.mat 
% %=============================================
% % for i=1:size(FinalCommunity_size_phase3,1)
% %     Final_comm(i,1)=i;
% %     Final_comm(i,2)=FinalCommunity_size_phase3(i,1);
% % end
% % FinalCommSorted= sort(Final_comm(:,2),'desc')
% % for i=1:size(OverlappingNodes_Size,1)
% %     Final_overlappingcomm(i,1)=i;
% %     Final_overlappingcomm(i,2)=OverlappingNodes_Size(i,1);
% % end
% % FinalOverlappingNodeCommSorted= sort(Final_overlappingcomm(:,2),'desc')
%===========================================
%%
clear all
clc
load ('preparing_Dataset.mat','finaldataset','Nodes','Input')
load ('LocalCommunities_Phase3_CaseStudy.mat','FinalCommunity_phase3','FinalCommunity_size_phase3')
load('LocalCommunities_Phase3_CaseStudy_Overlaps_infs.mat','Outlayer_L1','OverlappingNodes_L1','OverlappingNodes_Size_L1')
% ===============conevrting final communities node ID
for i=1:size(FinalCommunity_phase3,1)
    for j=1:FinalCommunity_size_phase3(i,1)
        [AA DD]=ismember(FinalCommunity_phase3(i,j),finaldataset(:,8));
        if AA==1
            FinalCommunity_mapped(i,j)=finaldataset(DD,1);
        else
            [EE FF]=ismember(FinalCommunity_phase3(i,j),finaldataset(:,9));
            FinalCommunity_mapped(i,j)=finaldataset(FF,2);
        end
    end
end
% ===============conevrting outliers node ID
for i=1:size(Outlayer_L1,1)-1
    [AA DD]=ismember(Outlayer_L1(i,1),finaldataset(:,8));
    if AA==1
        Outlayer_mapped(i,1)=finaldataset(DD,1);
    else
        [EE FF]=ismember(Outlayer_L1(i,1),finaldataset(:,9));
        Outlayer_mapped(i,1)=finaldataset(FF,2);
    end
end
% ===============conevrting OverlappingNodes ID
for i=1:size(OverlappingNodes_L1,1)
    for j=1:OverlappingNodes_Size_L1(i,1)
        [AA DD]=ismember(OverlappingNodes_L1(i,j),finaldataset(:,8));
        if AA==1
            OverlappingNodes_mapped(i,j)=finaldataset(DD,1);
        else
            [EE FF]=ismember(OverlappingNodes_L1(i,j),finaldataset(:,9));
            OverlappingNodes_mapped(i,j)=finaldataset(FF,2);
        end
    end
end
save FinalResultConverted.mat 
clear all
clc
load('FinalResultConverted.mat','FinalCommunity_mapped','Outlayer_mapped','OverlappingNodes_mapped')
save ResultforGHODSI.mat