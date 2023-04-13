% clear all
% clc
% % ================= Preparing ID of source and targets===========
% load ('DriverNetLayerOneFinal.mat','ResultDataLayerOne')
% load ('DriverNetLayerTwoFinalOrigDescNavi.mat','ResultDataLayerTwo')
% AllData=[ResultDataLayerOne(:,3);ResultDataLayerOne(:,4);ResultDataLayerTwo(:,3);ResultDataLayerTwo(:,4)];
% DataUniq(:,2)=unique(AllData(:,1),'rows');
% DataUniq(1:size(DataUniq,1),1)=1:size(DataUniq,1);
% for i=1:size(ResultDataLayerOne,1)
%     ResultDataLayerOne(i,5)=DataUniq(find(DataUniq(:,2)==ResultDataLayerOne(i,3)),1);
%     ResultDataLayerOne(i,6)=DataUniq(find(DataUniq(:,2)==ResultDataLayerOne(i,4)),1);
% end
% ResultDataLayerOne(1:size(ResultDataLayerOne,1),7)=1; % is the layer one lable
% for i=1:size(ResultDataLayerTwo,1)
%     ResultDataLayerTwo(i,5)=DataUniq(find(DataUniq(:,2)==ResultDataLayerTwo(i,3)),1);
%     ResultDataLayerTwo(i,6)=DataUniq(find(DataUniq(:,2)==ResultDataLayerTwo(i,4)),1);
% end
% ResultDataLayerTwo(1:size(ResultDataLayerTwo,1),7)=2; % is ther layer two lable;
% %%=================== Creatind the CaseStusy Input===============
% Input=[ResultDataLayerOne(:,5:7);ResultDataLayerTwo(:,5:7)];
% %===================Examples================================
% % Input=dlmread('I:\1- University Courses\2-PhD\3-Proposal\1-Social Network\Report\Coding\971007-heterogenousCDonDriverNet-TwoLayer\2-Coding_TWO_Layer_971020\980205_TWOLAYERCaseStudy\DataSets\CaseStudyDataSet.txt');
% Nodes=max(max(Input));
% [InRow, InClm]=size(Input);
% layerNo=Input(InRow,3); % Number of the layers
% NN=1;
% for i=1:InRow
%     if Input(i,3)==1
%         Input(i,InClm+1)=1; % weight of the edges
%         Input(i,InClm+2)=NN; % Index of the edges
%         NN=NN+1;
%     end
% end
% MM=1;
% for i=NN:InRow
%     if Input(InRow,3)==2
%         Input(i,InClm+1)=1; % weight of the edges
%         Input(i,InClm+2)=MM; % Index of the edges
%         MM=MM+1;
%     end
% end
% LayerInputNo(1,1)=NN-1;
% LayerInputNo(1,2)=(MM-1);
% tag=1;
% Input(i,6)=0;
% for i=1:LayerInputNo(1,1)
%     for j=1+LayerInputNo(1,1):LayerInputNo(1,1)+LayerInputNo(1,2)
%         if Input(i,1)==Input(j,1) && Input(i,2)==Input(j,2)
%             Input(i,6)=tag;
%             Input(j,6)=tag;
%             tag=tag+1;
%         else
%             
%         end
%     end
% end
% %==============calculation the adjacency matrix=================
% %Layer1
% Adja_Mat = zeros(Nodes,Nodes); % Adjacency Matrix of the "Layer #1"
% for i=1:LayerInputNo(1,1)
%     Adja_Mat(Input(i,1),Input(i,2))= Input(i,InClm+1);  % Adjacency Matrix of the "Layer #1"
%     %Adja_Mat(Input(i,2),Input(i,1))= Input(i,InClm+1); % if Graf is not Directed
%     s.L1(i,1)=Input(i,1);
%     t.L1(i,1)=Input(i,2);
% end
% %layer2
% Adja_Mat_L2 = zeros(Nodes,Nodes); % Adjacency Matrix Of the "layer #2"
% x=1;
% for i=LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2)
%     Adja_Mat_L2(Input(i,1),Input(i,2))= Input(i,InClm+1); % Adjacency Matrix Of the "layer #2"
%     %Adja_Mat_L2(Input(i,2),Input(i,1))= Input(i,InClm+1); % if Graf is not Directed
%     s.L2(x,1)=Input(i,1);
%     t.L2(x,1)=Input(i,2);
%     x=x+1;
% end
% %%======================Creating the "totoal Adj_Mat"======================
% % AA=Adja_Mat;
% % BB=Adja_Mat_L2;
% % Both_Mat_weighted=AA+BB;
% % Adja_Mat_Both0=sign(Both_Mat_weighted);
% % [Neww_Inpt,edgeCount]=Adjamat2Edge(Adja_Mat_Both0);
% R=0;
% New=zeros(size(Input,2));
% for i=1:LayerInputNo(1,1)
%     for j=LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2)
%         [AA BB]=ismember(Input(j,1:2),Input(i,1:2));
%         if sum(AA(1,:))==2
%             R=R+1;
%             New(R,:)=Input(j,:);
%             break
%         end
%     end
% end
% if New ~=0
%     Inters_in_L2=unique(New,'rows'); % mojod dar layer 2
%     [V W]=ismember(Input(LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2),5),Inters_in_L2(:,5));
% else
%     V= zeros(LayerInputNo(1,2),1);
% end
% S=0;
% for i=LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2)
%     if V(i-LayerInputNo(1,1),1)==0
%         S=S+1;
%         New_in_L2(S,:)= Input(i,1:5); % Naa mojood dar L2
%     end
% end
% x=LayerInputNo(1,1)+1;
% for i=1: size(New_in_L2,1)
%     New_in_L2(i,6)=x;
%     x=x+1;
% end
% Neww_Input=[Input(1:LayerInputNo,1:5);New_in_L2(:,1:4),New_in_L2(:,6) ];
% New_Input=sortrows(Neww_Input); % new s and t from two layers (edges from two layers)
% % adja_mat of both two layers
% Adja_Mat_both = zeros(Nodes,Nodes); % Adjacency Matrix of the both two layers
% for i=1:size(New_Input,1)
%     Adja_Mat_both(New_Input(i,1),New_Input(i,2))= New_Input(i,InClm+1);  % Adjacency Matrix of the both two layers
%     %Adja_Mat_both(New_Input(i,2),New_Input(i,1))= New_Input(i,InClm+1); % if Graf is not Directed
%     s.Lb(i,1)=New_Input(i,1);
%     t.Lb(i,1)=New_Input(i,2);
% end
% %%===========Calculating the "Centrality measure" for each Node============
% G.L1 = graph(s.L1,t.L1);
% CC.L1 = centrality(G.L1,'betweenness'); % could use other related measures such as 'eigenvector' or 'pagerank' 
% [a.C,a.D]= sort(CC.L1,'desc'); % Nodes sorted as Descending, D is Node name!
% 
% G.L2 = graph(s.L2,t.L2);
% CC.L2 = centrality(G.L2,'betweenness'); % could use other related measures such as 'eigenvector' or 'pagerank' 
% [a.C2,a.D2]= sort(CC.L2,'desc'); % Nodes sorted as Descending, D is Node name!
% 
% for i=1:LayerInputNo(1,1)
%     Adja_Mat(Input(i,2),Input(i,1))= Input(i,InClm+1); % if Graf is not Directed
% end
% %layer2
% for i=LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2)
%     Adja_Mat_L2(Input(i,2),Input(i,1))= Input(i,InClm+1); % if Graf is not Directed
% end
% %for both two layers
% for i=1:size(New_Input,1)
%     Adja_Mat_both(New_Input(i,2),New_Input(i,1))= New_Input(i,InClm+1); % if Graf is not Directed
% end
% % for i=1:size(Neww_Inpt,1)
% %     Adja_Mat_Both0(Neww_Inpt(i,2),Neww_Inpt(i,1))= Neww_Inpt(i,InClm+1); % if Graf is not Directed
% % end
% save preparing_CaseStudy.mat 
% %
% %===============Phase 0: "Edge probing"===========
% clear all
% clc
% load ('preparing_CaseStudy.mat','Adja_Mat','Adja_Mat_L2','a')
% %--------------% sorting edges for creating local communities--------------
% e_phase0=cputime;
% a1=a.D;
% a2=a.D2;
% e_phase0_1=cputime;
% [edge_Sorted,NodeEdgesL1] = NodeProbing_Func(a1,Adja_Mat); % Finding the Sorted Edges in Layer 1
% Time_phase0_1=cputime-e_phase0_1;
% e_phase1_2=cputime;
% [edge_Sorted_L2 NodeEdgesL2]= NodeProbing_Func(a2,Adja_Mat_L2); % Finding the Sorted Edges in Layer 2
% Time_phase0_2=cputime-e_phase1_2;
% Time_p0=Time_phase0_1+Time_phase0_2;
% save LocalCommunities_Phase0_CaseStudy.mat 
% %
% %=================== Pahse 1: Local Community Detection===========
% clear all
% clc
% load ('LocalCommunities_Phase0_CaseStudy.mat','edge_Sorted')
% load ('preparing_CaseStudy.mat','Nodes','Adja_Mat','LayerInputNo','Input')
% Resolution=4; %the greater the Resolution the lower the comms
% e_phase1_1=cputime;
% [local_Community,local_Community_size]=LocalCommunity_Phase1_Func(edge_Sorted,Adja_Mat,Nodes,Resolution); % Calculate the Lcoal ommunity of the Layer #1!!!!!!!!!!!!!!!!!!!!!!
% Time_phase1_1=cputime-e_phase1_1;
% L1_Input=Input(1:LayerInputNo(1,1),:);
% [EQ_Local_L1]=EQ_Measure_Func(L1_Input,Nodes,local_Community,local_Community_size,Adja_Mat);% calculae the Q measure for the Local communities in Layer#1
% save LocalCommunities_Phase1_1_CaseStudy.mat 
% clear all
% clc
% load ('LocalCommunities_Phase0_CaseStudy.mat','edge_Sorted_L2')
% load ('preparing_CaseStudy.mat','Nodes','Adja_Mat_L2','LayerInputNo','Input')
% load ('LocalCommunities_Phase1_1_CaseStudy.mat','Time_phase1_1')
% Resolution=4;
% e_phase1_2=cputime;
% [local_Community2,local_Community_size_L2]=LocalCommunity_Phase1_Func(edge_Sorted_L2,Adja_Mat_L2,Nodes,Resolution); % Calculate the Lcoal ommunity of the Layer #1!!!!!!!!!!!!!!!!!!!!!!
% Time_phase1_2=cputime-e_phase1_2;
% L2_Input=Input(LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2),:);
% [EQ_Local_L2]=EQ_Measure_Func(L2_Input,Nodes,local_Community2,local_Community_size_L2,Adja_Mat_L2);% calculae the Q measure for the Local communities in Layer#2
% Time_p1=Time_phase1_1+Time_phase1_2;
% save LocalCommunities_Phase1_2_CaseStudy.mat 
% %
% %==============finding the overlapping nodes proportion in comms=====
% clear all
% clc
% load ('LocalCommunities_Phase1_1_CaseStudy.mat','local_Community','local_Community_size')
% load ('LocalCommunities_Phase1_2_CaseStudy.mat','local_Community2','local_Community_size_L2')
% load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input')
% %======================finding outliers in comms=====================
% [Outlayer_L1,NodesExtincom_L1]=OutlierDetectionFunc(local_Community,local_Community_size,Nodes);
% [Outlayer_L2,NodesExtincom_L2]=OutlierDetectionFunc(local_Community2,local_Community_size_L2,Nodes);
% %===========finding the overlapping nodes=================
% [OverlappingNodes_L1,OverlappingNodes_Size_L1,OverlappingNodes_proportion_L1,proportion_L1,comm_count_node_L1,comm_proportion_node_L1]=overlappingDetectionFunc(local_Community,local_Community_size,NodesExtincom_L1);
% [OverlappingNodes_L2,OverlappingNodes_Size_L2,OverlappingNodes_proportion_L2,proportion_L2,comm_count_node_L2,comm_proportion_node_L2]=overlappingDetectionFunc(local_Community2,local_Community_size_L2,NodesExtincom_L2);
% save LocalCommunities_Phase1_CaseStudy_Overlaps_infs.mat 
% %
% %===============Phase 2: Intra-Layer Community Merging===========================
% clear all
% clc
% load ('LocalCommunities_Phase1_1_CaseStudy.mat','local_Community','Input','LayerInputNo')
% load ('preparing_CaseStudy.mat','Adja_Mat','Nodes','InClm')
% % Final_Community=local_Community;
% % Final_Community_size=size(local_Community,2);
% %FF=zeros(10,size(local_Community,1),Nodes);
% m=1;
% for i=1:size(local_Community,1)
%     AA(i,1)=size(local_Community,2)-sum(local_Community(i,:)==0);
%     if AA(i,1)<Nodes
%         CC(m,1:size(local_Community,2)-sum(local_Community(i,:)==0))=local_Community(i,1:size(local_Community,2)-sum(local_Community(i,:)==0));
%         m=m+1;
%     end
% end
% F=0;
% Alfa=0.6;
% Gama=0.9;
% Beta=0.35;
% %%================================================================
% e_phase2_1=cputime;
% [Final_Community,Final_Community_size,MergedCommunities]=CommunityMerging_Phase2_Func(CC,Adja_Mat,Nodes,LayerInputNo(1,1),Input,InClm,F,Alfa,Gama,Beta); % Id od merged communities
% Time_phase2_1=cputime-e_phase2_1;
% L1_Input=Input(1:LayerInputNo(1,1),:);
% [EQ_l1]=EQ_Measure_Func(L1_Input,Nodes,Final_Community,Final_Community_size,Adja_Mat);% calculae the Q measure for the communities in Layer#1
% save LocalCommunities_Phase2_1_Casestudy.mat 
% clear all
% clc
% load ('LocalCommunities_Phase1_2_CaseStudy.mat','local_Community2','Input','LayerInputNo')
% load ('preparing_CaseStudy.mat','Adja_Mat_L2','Nodes','InClm')
% load ('LocalCommunities_Phase2_1_Casestudy.mat','Time_phase2_1')
% m=1;
% for i=1:size(local_Community2,1)
%     DD(i,1)=size(local_Community2,2)-sum(local_Community2(i,:)==0);
%     if DD(i,1)<Nodes
%         BB(m,1:size(local_Community2,2)-sum(local_Community2(i,:)==0))=local_Community2(i,1:size(local_Community2,2)-sum(local_Community2(i,:)==0));
%         m=m+1;
%     end
% end
% Alfa=0.6;
% Gama=0.9;
% Beta=0.35;
% F=LayerInputNo(1,1);
% e_phase2_2=cputime;
% [Final_Community_L2_out,Final_Community_size_L2_out,MergedCommunities_L2]=CommunityMerging_Phase2_Func(BB,Adja_Mat_L2,Nodes,LayerInputNo(1,2),Input,InClm,F,Alfa,Gama,Beta);
% Time_phase2_2=cputime-e_phase2_2;
% Final_Community_L2=Final_Community_L2_out;
% Final_Community_size_L2=Final_Community_size_L2_out;
% L2_Input=Input(LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2),:);
% [EQ_l2]=EQ_Measure_Func(L2_Input,Nodes,Final_Community_L2,Final_Community_size_L2,Adja_Mat_L2);% calculae the Q measure for the communities in Layer#2
% Time_p2=Time_phase2_2+Time_phase2_1%
% save LocalCommunities_Phase2_2_Casestudy.mat 
% %
% %==============finding the overlapping nodes proportion in comms=====
% clear all
% clc
% load ('LocalCommunities_Phase2_1_CaseStudy.mat','Final_Community','Final_Community_size')
% load ('LocalCommunities_Phase2_2_CaseStudy.mat','Final_Community_L2','Final_Community_size_L2')
% load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input')
% %======================finding outliers in comms=====================
% [Outlayer_L1,NodesExtincom_L1]=OutlierDetectionFunc(Final_Community,Final_Community_size,Nodes);
% [Outlayer_L2,NodesExtincom_L2]=OutlierDetectionFunc(Final_Community_L2,Final_Community_size_L2,Nodes);
% %===========finding the overlapping nodes=================
% [OverlappingNodes_L1,OverlappingNodes_Size_L1,OverlappingNodes_proportion_L1,proportion_L1,comm_count_node_L1,comm_proportion_node_L1]=overlappingDetectionFunc(Final_Community,Final_Community_size,NodesExtincom_L1);
% [OverlappingNodes_L2,OverlappingNodes_Size_L2,OverlappingNodes_proportion_L2,proportion_L2,comm_count_node_L2,comm_proportion_node_L2]=overlappingDetectionFunc(Final_Community_L2,Final_Community_size_L2,NodesExtincom_L2);
% save LocalCommunities_Phase2_CaseStudy_Overlaps_infs.mat 
% clc
% %
% %========================Phase 3: Refining communities ================================
% % Final_Community_Phase3=Final_Community;
% % Final_Community_size_Phase3=Final_Community_size;
% % Adja_Mat_both=Adja_Mat;
% % New_Input=Input(1:LayerInputNo(1,1),:);
% clear all
% clc
% load ('LocalCommunities_Phase2_1_Casestudy.mat','Final_Community','Final_Community_size')
% load ('preparing_CaseStudy.mat','Adja_Mat','Nodes','LayerInputNo','Input')
% e_phase3_1_1=cputime;
% [FinalCommunity_phase3_1,FinalCommunity_size_phase3_1,NodesExistinAllComs,Outlayer]=CommunityRefining_phase4_Func(Final_Community,Final_Community_size,Adja_Mat,Nodes); % Final Comm Phase 3 for layer1
% Time_p3_1=cputime-e_phase3_1_1;
% %-------------------Calculating the EQ measure-------------------
% L1_Input=Input(1:LayerInputNo(1,1),:);
% [EQ_Phase3_1]=EQ_Measure_Func(L1_Input,Nodes,FinalCommunity_phase3_1,FinalCommunity_size_phase3_1,Adja_Mat);
% save LocalCommunities_Phase3_1_Casestudy.mat 
% clear all
% clc
% load ('LocalCommunities_Phase2_2_Casestudy.mat','Final_Community_L2','Final_Community_size_L2')
% load ('preparing_CaseStudy.mat','Adja_Mat_L2','Nodes','LayerInputNo','Input')
% load ('LocalCommunities_Phase3_1_Casestudy.mat','Time_p3_1')
% e_phase3_1_2=cputime;
% [FinalCommunity_phase3_2,FinalCommunity_size_phase3_2,NodesExistinAllComs,Outlayer]=CommunityRefining_phase4_Func(Final_Community_L2,Final_Community_size_L2,Adja_Mat_L2,Nodes); % Final Comm Phase 3 for layer2
% Time_p3_2=cputime-e_phase3_1_2;
% Time_p3=Time_p3_1+Time_p3_2;
% %-------------------Calculating the EQ measure-------------------
% L2_Input=Input(LayerInputNo(1,1)+1:LayerInputNo(1,1)+LayerInputNo(1,2),:);
% [EQ_Phase3_2]=EQ_Measure_Func(L2_Input,Nodes,FinalCommunity_phase3_2,FinalCommunity_size_phase3_2,Adja_Mat_L2);
% save LocalCommunities_Phase3_2_Casestudy.mat 
% %
% %==============finding the overlapping nodes proportion in comms=====
% clear all
% clc
% load ('LocalCommunities_Phase3_1_CaseStudy.mat','FinalCommunity_phase3_1','FinalCommunity_size_phase3_1')
% load ('LocalCommunities_Phase3_2_CaseStudy.mat','FinalCommunity_phase3_2','FinalCommunity_size_phase3_2')
% load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input')
% %======================finding outliers in comms=====================
% [Outlayer_L1,NodesExtincom_L1]=OutlierDetectionFunc(FinalCommunity_phase3_1,FinalCommunity_size_phase3_1,Nodes);
% [Outlayer_L2,NodesExtincom_L2]=OutlierDetectionFunc(FinalCommunity_phase3_2,FinalCommunity_size_phase3_2,Nodes);
% %===========finding the overlapping nodes=================
% [OverlappingNodes_L1,OverlappingNodes_Size_L1,OverlappingNodes_proportion_L1,proportion_L1,comm_count_node_L1,comm_proportion_node_L1]=overlappingDetectionFunc(FinalCommunity_phase3_1,FinalCommunity_size_phase3_1,NodesExtincom_L1);
% [OverlappingNodes_L2,OverlappingNodes_Size_L2,OverlappingNodes_proportion_L2,proportion_L2,comm_count_node_L2,comm_proportion_node_L2]=overlappingDetectionFunc(FinalCommunity_phase3_2,FinalCommunity_size_phase3_2,NodesExtincom_L2);
% save LocalCommunities_Phase3_CaseStudy_Overlaps_infs.mat 
% clc
% %
% %======Phase4: "Inter-Layers community Mergging with NCWOS measure" ============
% clear all
% clc
% load ('LocalCommunities_Phase3_1_Casestudy.mat','FinalCommunity_phase3_1','FinalCommunity_size_phase3_1')
% load ('LocalCommunities_Phase3_2_Casestudy.mat','FinalCommunity_phase3_2','FinalCommunity_size_phase3_2')
% load ('preparing_CaseStudy.mat','New_Input','InClm','Nodes','Adja_Mat_both')
% load ('LocalCommunities_Phase2_2_Casestudy.mat','Alfa','Gama','Beta')
% [TwoLayersCommunity_Local,TwoLayersCommunity_Local_Size]= TwoLayers_Intersect_Func(FinalCommunity_phase3_1,FinalCommunity_size_phase3_1,FinalCommunity_phase3_2,FinalCommunity_size_phase3_2); % Intersect two layers communities
% Len=max(max(FinalCommunity_size_phase3_1(:,1)),max(FinalCommunity_size_phase3_2(:,1)));
% TwoLayerLocalCommunity=0;
% for i=1: size(FinalCommunity_size_phase3_1,1)
%     TwoLayerLocalCommunity(i,1:FinalCommunity_size_phase3_1(i,1))=FinalCommunity_phase3_1(i,1:FinalCommunity_size_phase3_1(i,1));
% end
% sizeTL=size(TwoLayerLocalCommunity,1);
% for i=1: size(FinalCommunity_size_phase3_2,1)
%     TwoLayerLocalCommunity(i+sizeTL,1:FinalCommunity_size_phase3_2(i,1))=FinalCommunity_phase3_2(i,1:FinalCommunity_size_phase3_2(i,1));
% end
% F=0;
% New_Input_Size=size(New_Input,1);
% e_phase4_1=cputime;
% [Final_Community_Phase4,Final_Community_size_Phase4,MergedCommunities]=CommunityMerging_Phase2_Func(TwoLayerLocalCommunity,Adja_Mat_both,Nodes,New_Input_Size,New_Input,InClm,F,Alfa,Gama,Beta); % dar jadvale "merged communities" ta community number 60 baraye layer 1 ast az 61 baraye layer 2 hast
% Time_phase4_1=cputime-e_phase4_1;
% Time_p4=Time_phase4_1;
% [EQ_Phase4]=EQ_Measure_Func(New_Input,Nodes,Final_Community_Phase4,Final_Community_size_Phase4,Adja_Mat_both);
% save LocalCommunities_Phase4_Casestudy.mat 
% %
% %==============finding the overlapping nodes proportion in comms=====
% clear all
% clc
% load ('LocalCommunities_Phase4_Casestudy.mat','Final_Community_Phase4','Final_Community_size_Phase4')
% load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input')
% %======================finding outliers in comms=====================
% [Outlayer,NodesExtincom]=OutlierDetectionFunc(Final_Community_Phase4,Final_Community_size_Phase4,Nodes);
% %===========finding the overlapping nodes=================
% [OverlappingNodes,OverlappingNodes_Size,OverlappingNodes_proportion,proportion,comm_count_node,comm_proportion_node]=overlappingDetectionFunc(Final_Community_Phase4,Final_Community_size_Phase4,NodesExtincom);
% save LocalCommunities_Phase4_CaseStudy_Overlaps_infs.mat 
% clc
% %
% % %=================== phase4-prime: "Inter-Layers community Mergging" with THE REDUNDANCY measure======================
% % clear all
% % clc
% % load ('LocalCommunities_Phase3_1_Casestudy.mat','FinalCommunity_phase3_1','FinalCommunity_size_phase3_1')
% % load ('LocalCommunities_Phase3_2_Casestudy.mat','FinalCommunity_phase3_2','FinalCommunity_size_phase3_2')
% % load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input','InClm')
% % [Ro_Final_Community]=CommuityMergingwithRedundancy(FinalCommunity_phase3_1,FinalCommunity_size_phase3_1,FinalCommunity_phase3_2,FinalCommunity_size_phase3_2,Adja_Mat,Adja_Mat_L2,LayerInputNo,Input,InClm);
% % save FinalCommunity_withRedundacyMeasure_Casestudy.mat 
% %
% %======================        resutls       ===========================
% clear all
% clc
% load ('LocalCommunities_Phase0_CaseStudy.mat','Time_p0')
% load ('LocalCommunities_Phase1_1_CaseStudy.mat','EQ_Local_L1')
% load ('LocalCommunities_Phase1_2_CaseStudy.mat','EQ_Local_L2','Time_p1')
% load ('LocalCommunities_Phase2_1_Casestudy.mat','EQ_l1')
% load ('LocalCommunities_Phase2_2_Casestudy.mat','EQ_l2','Time_p2')
% load ('LocalCommunities_Phase3_1_Casestudy.mat','EQ_Phase3_1')
% load ('LocalCommunities_Phase3_2_Casestudy.mat','EQ_Phase3_2','Time_p3')
% load ('LocalCommunities_Phase4_Casestudy.mat','EQ_Phase4','Time_p4')
% EQ_Local_L1
% EQ_Local_L2
% EQ_l1
% EQ_l2
% EQ_Phase3_1
% EQ_Phase3_2
% EQ_Phase4
% Time_p0
% Time_p1
% Time_p2
% Time_p3
% Time_p4
% Total_CPUT_TIME=Time_p0+Time_p1+Time_p2+Time_p3+Time_p4
% save Results_Casestudy.mat
% % %
% % %=================== CALCULATING THE REDUNDANCY INDEX======================
% clear all
% clc
% load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input')
% load ('LocalCommunities_Phase4_Casestudy.mat','Final_Community_Phase4','Final_Community_size_Phase4')
% [ro_Total,ro,pc,pc_2hat,proportion]=RedundancyFunc(Final_Community_Phase4,Final_Community_size_Phase4,Adja_Mat,Adja_Mat_L2,LayerInputNo(1,1),LayerInputNo(1,2),Input,InClm);
% save Redundacy_Casestudy.mat 
% %
% %=========================== ploting the graph=============================
% clear all
% clc
% load ('preparing_CaseStudy.mat','InClm','Nodes','Adja_Mat','Adja_Mat_L2','LayerInputNo','Input','New_Input')
% load ('LocalCommunities_Phase4_Casestudy.mat','Final_Community_Phase4','Final_Community_size_Phase4')
% for i=1:size(Final_Community_Phase4,1)
%     k=1;
%     for l=1:size(Final_Community_Phase4,1)
%         for j=1:Final_Community_size_Phase4(i,1)
%             if ismember(Final_Community_Phase4(i,j),Final_Community_Phase4(l,:))==1 && i~=l
%                 OL(i,k)=Final_Community_Phase4(i,j);
%                 k=k+1;
%             end
%         end
%     end
%     
% end
% OverlappingNodes=zeros(size(Final_Community_Phase4,1),size(OL,2));
% for i=1:size(Final_Community_Phase4,1)
%       a=unique(OL(i,:));
%       b=sort(a,'desc');
%       OverlappingNodes(i,1:size(b,2))=b; % Overlapping Nodes in each Community
%       OverlappingNodes_Size(i,1)=size(OverlappingNodes,2)-sum(OverlappingNodes(i,:)==0); % Number of the overlapping nodes in each Community
% end
% 
% s = New_Input(:,1); 
% t = New_Input(:,2); 
% G = graph(s,t);
% D = degree(G); % Degree of the nodes!
% h = plot(G);
% aa=3;
% bb=7;
% A=Final_Community_Phase4(aa,1:Final_Community_size_Phase4(aa,1)); % Frist Community Nodes
% B=Final_Community_Phase4(bb,1:Final_Community_size_Phase4(bb,1)); % Second Community Nodes
% %C= Outlayer(:,1);                                             % is the Outlayer Nodes
% highlight(h,A,'NodeColor','g')
% highlight(h,B,'NodeColor','r')
% %highlight(h,C,'NodeColor','y')
% % for j=1:size(OverlappingNodes_Size,1)
% %     X=OverlappingNodes(j,1:OverlappingNodes_Size(j,1));
% %     highlight(h,X)
% % end
% % save LocalCommunities_Phase5-Casestudy.mat 
% %=============================================
% for i=1:size(Final_Community_size_Phase4,1)
%     Final_comm(i,1)=i;
%     Final_comm(i,2)=Final_Community_size_Phase4(i,1);
% end
% FinalCommSorted= sort(Final_comm(:,2),'desc')
% for i=1:size(OverlappingNodes_Size,1)
%     Final_overlappingcomm(i,1)=i;
%     Final_overlappingcomm(i,2)=OverlappingNodes_Size(i,1);
% end
% FinalOverlappingNodeCommSorted= sort(Final_overlappingcomm(:,2),'desc')
% %%
% %=========== extracting degargoonsazha==============================
% clear all
% clc
% load ('LocalCommunities_Phase4_Casestudy.mat','Final_Community_Phase4','Final_Community_size_Phase4','FinalCommunity_phase3_1','FinalCommunity_phase3_2','FinalCommunity_size_phase3_1','FinalCommunity_size_phase3_2','MergedCommunities')
% MergedCommunities_Final=MergedCommunities;
% MergedCommunities=0;
% FinalCommunity_phase4=Final_Community_Phase4;
% Final_Community_Phase4=0;
% load ('LocalCommunities_Phase2_1_CaseStudy.mat','Final_Community','Final_Community_size','MergedCommunities','local_Community')
% MergedCommunities_L1=MergedCommunities;
% MergedCommunities=0;
% FinalCommunity_phase2_1=Final_Community;
% Final_Community=0;
% FinalCommunity_size_1=Final_Community_size;
% Final_Community_size=0;
% load ('LocalCommunities_Phase2_2_CaseStudy.mat','Final_Community_L2','Final_Community_size_L2','MergedCommunities_L2','local_Community2')
% FinalCommunity_phase2_2=Final_Community_L2;
% Final_Community_L2=zeros(1,1);
% FinalCommunity_size_L2=Final_Community_size_L2;
% Final_Community_size_L2=0;
% 
% for i=1:size(MergedCommunities_Final,1)
%     MergedCommunities_Final_Size(i,1)=size(MergedCommunities_Final,2)-sum(MergedCommunities_Final(i,:)==0);
% end
% for i=1:size(MergedCommunities_L1,1)
%     MergedCommunities_L1_Size(i,1)=size(MergedCommunities_L1,2)-sum(MergedCommunities_L1(i,:)==0);
% end
% for i=1:size(MergedCommunities_L2,1)
%     MergedCommunities_L2_Size(i,1)=size(MergedCommunities_L2,2)-sum(MergedCommunities_L2(i,:)==0);
% end
% w=1;
% x=1;
% Merged_LocalCommunity_LC1=zeros(size(MergedCommunities_Final,1),1);
% Merged_LocalCommunity_LC2=zeros(size(MergedCommunities_Final,1),1);
% for i=1:size(MergedCommunities_Final,1)
%     for j=1:MergedCommunities_Final_Size(i,1)
%         if MergedCommunities_Final(i,j)<=size(FinalCommunity_size_phase3_1,1)
%             for k=1:MergedCommunities_L1_Size(MergedCommunities_Final(i,j),1)
%                 Merged_LocalCommunity_LC1(i,w)= MergedCommunities_L1(MergedCommunities_Final(i,j),k);% local communities' ID from Layer #1 in final comms
%                 w=w+1;
%             end
%         else
%             for l=1:MergedCommunities_L2_Size(MergedCommunities_Final(i,j)-size(FinalCommunity_size_phase3_1,1),1) 
%                 Merged_LocalCommunity_LC2(i,x)= MergedCommunities_L2(MergedCommunities_Final(i,j)-size(FinalCommunity_size_phase3_1,1),l);% local communities' ID from Layer #2 in final comms
%                 x=x+1;
%             end
%         end
%     end
%     w=1;
%     x=1;
% end
% 
% for i=1:size(Merged_LocalCommunity_LC1,1)
%     Merged_LocalCommunity_LC1_Size(i,1)=size(Merged_LocalCommunity_LC1,2)-sum(Merged_LocalCommunity_LC1(i,:)==0);
% end
% for i=1:size(Merged_LocalCommunity_LC2,1)
%     Merged_LocalCommunity_LC2_Size(i,1)=size(Merged_LocalCommunity_LC2,2)-sum(Merged_LocalCommunity_LC2(i,:)==0);
% end
% 
% FinalCommunity_phase4_LC1_Count=zeros(size(FinalCommunity_phase4,1),size(FinalCommunity_phase4,2));
% for i=1:size(FinalCommunity_phase4,1)
%     for j=1:Final_Community_size_Phase4(i,1)
%         w=0;
%         for k=1:Merged_LocalCommunity_LC1_Size(i,1)
%             if ismember(FinalCommunity_phase4(i,j),local_Community(Merged_LocalCommunity_LC1(i,k),:))==1 
%                 w=w+1;
%                 FinalCommunity_phase4_LC1_Count(i,j)=w;
%             end
%         end
%     end
% end
% 
% FinalCommunity_phase4_LC2_Count=zeros(size(FinalCommunity_phase4,1),size(FinalCommunity_phase4,2));
% for i=1:size(FinalCommunity_phase4,1)
%     for j=1:Final_Community_size_Phase4(i,1)
%         w=0;
%         for k=1:Merged_LocalCommunity_LC2_Size(i,1)
%             if ismember(FinalCommunity_phase4(i,j),local_Community2(Merged_LocalCommunity_LC2(i,k),:))==1 
%                 w=w+1;
%                 FinalCommunity_phase4_LC2_Count(i,j)=w;
%             end
%         end
%     end
% end
% 
% FinalCommunity_phase4_Total_Nodes_Count=FinalCommunity_phase4_LC1_Count +FinalCommunity_phase4_LC2_Count; %%  the local comms node_counter in each final comm
% 
% %------- sort the nodes based on their count or proportion
% for u=1:size(FinalCommunity_phase4_Total_Nodes_Count,1)
%     [S1 I1]=sort(FinalCommunity_phase4_Total_Nodes_Count(u,:),'desc');
%     ISortednodes(u,1:size(S1,2)-sum(S1(1,:)==0))=I1(1,1:size(S1,2)-sum(S1(1,:)==0)); % index of the sorted nodes 
%     Sortednodes(u,1:size(S1,2)-sum(S1(1,:)==0))=S1(1,1:size(S1,2)-sum(S1(1,:)==0)); % sorted nodes counters
%     Sum(u,1)=sum(Sortednodes(u,:));
% end
% clc
% param=0.9; %((counternode/Sum)*100)percentage of the Nodes are used to create the (param) percentage  of local communities. 
% counter=zeros(size(FinalCommunity_phase4_Total_Nodes_Count,1),1);
% for i=1:size(FinalCommunity_phase4_Total_Nodes_Count,1)
%     for j=1:size(FinalCommunity_phase4_Total_Nodes_Count,2)
%         counter(i,1)=counter(i,1)+Sortednodes(i,j);
%         if counter(i,1)<=param*Sum(i,1)
%             counternode(i,1)=j;
%         else
%             break
%         end
%     end
% end
% for i=1:size(counternode,1)
%     for j=1:counternode(i,1)
%         Volunteer_Nodes(i,j)=FinalCommunity_phase4(i,ISortednodes(i,j)); % nodes in each cummunity which can be used in deffusion model!!!!!!!!!!!!!
%     end
% end
% for i=1:size(Volunteer_Nodes,1)
%     Volunteer_Nodes_Size(i,1)=size(Volunteer_Nodes,2)-sum(Volunteer_Nodes(i,:)==0);
% end
% save Volunteer_Nodes_Extraction_Casestudy.mat 

%%%======================= extracting edges of the volunteer nodes and calculating the alpha cut value for each volunteer nodes========================
clear all
clc
load ('Volunteer_Nodes_Extraction_Casestudy.mat','Volunteer_Nodes','Volunteer_Nodes_Size')
load ('LocalCommunities_Phase4_Casestudy.mat','Final_Community_Phase4','Final_Community_size_Phase4')
load ('TwoLayersFuzzyRalationValue.mat','FuzzyDataLayerOne','FuzzyDataLayerTwo')
load ('LocalCommunities_Phase0_CaseStudy.mat','Adja_Mat','Adja_Mat_L2')
load ('preparing_CaseStudy.mat','ResultDataLayerOne','ResultDataLayerTwo')
w=1;
for i=1:size(FuzzyDataLayerOne,1)
    for j=1:size(ResultDataLayerOne,1)
        if FuzzyDataLayerOne(i,1)== ResultDataLayerOne(j,1) && FuzzyDataLayerOne(i,2)==ResultDataLayerOne(j,2)
            FuzzyDataLayerOne(i,9)=ResultDataLayerOne(j,5); %extracting the used Index of the source (used in coding)
            FuzzyDataLayerOne(i,10)=ResultDataLayerOne(j,6); %extracting the used Index of the target (used in coding)
        end
    end
end

for i=1:size(Volunteer_Nodes_Size,1)
    for j=1:Volunteer_Nodes_Size(i,1)
        A=Volunteer_Nodes(i,j);
        for k=1:Final_Community_size_Phase4(i,1)
            B=Final_Community_Phase4(i,k);
            for m=1:size(FuzzyDataLayerOne,1)
                if A== FuzzyDataLayerOne(m,9) && B==FuzzyDataLayerOne(m,10)
                    Volunteer_edges_L1(w,1)=A;
                    Volunteer_edges_L1(w,2)=B;
                    Volunteer_edges_L1(w,3)=FuzzyDataLayerOne(m,5);
                    Volunteer_edges_L1(w,4)=FuzzyDataLayerOne(m,6);
                    w=w+1;
                    break;
                elseif A== FuzzyDataLayerOne(m,10) && B==FuzzyDataLayerOne(m,9)
                    Volunteer_edges_L1(w,1)=A;
                    Volunteer_edges_L1(w,2)=B;
                    Volunteer_edges_L1(w,3)=FuzzyDataLayerOne(m,7);
                    Volunteer_edges_L1(w,4)=FuzzyDataLayerOne(m,8);
                    w=w+1;
                    break;
                end
            end
        end
    end
end
for i=1:size(FuzzyDataLayerTwo,1)
    for j=1:size(ResultDataLayerTwo,1)
        if FuzzyDataLayerTwo(i,1)== ResultDataLayerTwo(j,1) && FuzzyDataLayerTwo(i,2)==ResultDataLayerTwo(j,2)
            FuzzyDataLayerTwo(i,9)=ResultDataLayerTwo(j,5); %extracting the used Index of the source (used in coding)
            FuzzyDataLayerTwo(i,10)=ResultDataLayerTwo(j,6);%extracting the used Index of the target (used in coding)
        end
    end
end
x=1;
for i=1:size(Volunteer_Nodes_Size,1)
    for j=1:Volunteer_Nodes_Size(i,1)
         A=Volunteer_Nodes(i,j);
        for k=1:Final_Community_size_Phase4(i,1)       
            B=Final_Community_Phase4(i,k);
            for m=1:size(FuzzyDataLayerTwo,1)
                if A== FuzzyDataLayerTwo(m,9) && B==FuzzyDataLayerTwo(m,10)
                    Volunteer_edges_L2(x,1)=A;
                    Volunteer_edges_L2(x,2)=B;
                    Volunteer_edges_L2(x,3)=FuzzyDataLayerTwo(m,5);
                    Volunteer_edges_L2(x,4)=FuzzyDataLayerTwo(m,6);
                    x=x+1;
                elseif A== FuzzyDataLayerTwo(m,10) && B==FuzzyDataLayerTwo(m,9)
                    Volunteer_edges_L2(x,1)=A;
                    Volunteer_edges_L2(x,2)=B;
                    Volunteer_edges_L2(x,3)=FuzzyDataLayerTwo(m,7);
                    Volunteer_edges_L2(x,4)=FuzzyDataLayerTwo(m,8);
                    x=x+1;
                end
            end
        end
    end
end
%=============== calculating the alpha cut value for fuzzy relation=====
alpha=0.6;
for i=1:size(Volunteer_edges_L1,1)
    Volunteer_edges_L1(i,5)= Volunteer_edges_L1(i,3)+alpha*(Volunteer_edges_L1(i,4)-Volunteer_edges_L1(i,3)); % clmn 5=alpha cut value
end
for i=1:size(Volunteer_edges_L2,1)
    Volunteer_edges_L2(i,5)= Volunteer_edges_L2(i,3)+alpha*(Volunteer_edges_L2(i,4)-Volunteer_edges_L2(i,3)); % clmn 5=alpha cut value
end
save VolunteerNodes_Edges_AlphaCutValues_Casestudy.mat 
clc
clear all
load ('VolunteerNodes_Edges_AlphaCutValues_Casestudy.mat','Volunteer_Nodes','Volunteer_Nodes_Size','Volunteer_edges_L1','Volunteer_edges_L2','FuzzyDataLayerOne','FuzzyDataLayerTwo','Adja_Mat','Adja_Mat_L2','Final_Community_Phase4','Final_Community_size_Phase4')
save FinalResult_for_OR_Modeling_Casestudy.mat 
