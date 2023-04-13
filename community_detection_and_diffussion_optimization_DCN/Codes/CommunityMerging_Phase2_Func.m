function [Final_Community,Final_Community_size,Community]=CommunityMerging_Phase2_Func(LC,AdjaMat,Nodes,LayerInputNo,Input,InClm,F,Alfa,Gama,Beta) %LC=local_Community
%==================== WOS_NCA Algorithm

[WOS,NCA,NCA_p1,NCA_p2]=WOS_NCA_Function(LC,AdjaMat,Input,InClm,LayerInputNo,F,Nodes,Alfa,Gama); % my new WOS function!!!!!!!!!

%===================== Final Community of Pahse 2==========================
%[Final_Community,Final_Community_size]=FinalCommunity_Phase2(LC,WOS,Beta);
%% "Calculatingthe Final Community of Phase2" with WOS_NCA
[Final_Community,Final_Community_size,Community]=FinalCommunity_Phase2_WOS_NCA_Itrtv(LC,WOS,Beta,AdjaMat,Input,InClm,LayerInputNo,F,Nodes,Alfa,Gama); % "Calculatingthe Final Community of Phase2" with WOS_NCA_Itrtv