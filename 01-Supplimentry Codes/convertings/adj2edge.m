clc
clear all
Input=dlmread('I:\1- University Courses\2-PhD\3-Proposal\1-Coding\01-Supplimentry Codes\convertings\TwoLayerDataSets\DataSet\KAPFM.txt');
adjl1=Input(1:(size(Input,1)/2),1:size(Input,2));
adjl2=Input((size(Input,1)/2)+1:size(Input,1),1:size(Input,2));
n=length(adjl1); % number of nodes in L1
g=length(adjl2); % number of nodes in L2
[EL1a EL1b]=find(adjl1>0)
edges_L1=[EL1a EL1b]
edges_L1(:,3)=1
edges_L1(:,4)=EL1a+EL1b;
[EL2a EL2b]=find(adjl2>0)
edges_L2=[EL2a EL2b]
edges_L2(:,3)=2
edges_L2(:,4)=EL2a+EL2b;
edgeCount_L1=length(edges_L1)
edgeCount_L2=length(edges_L2)
Total_Edge_Count=edgeCount_L1+edgeCount_L2
edgeSet=[edges_L1;edges_L2];
% edgeSet(:,3)=edgeSet(:,1)+edgeSet(:,2);