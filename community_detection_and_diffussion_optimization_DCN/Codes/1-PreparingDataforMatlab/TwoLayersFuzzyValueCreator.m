% clear all
% clc
% InputL1=dlmread('I:\1- University Courses\2-PhD\3-Proposal\1-Coding\6-CodingDCNTLDiffusionOpt\Codes\1-PreparingDataforMatlab\DNetWeightedL1.txt');
% InputL2=dlmread('I:\1- University Courses\2-PhD\3-Proposal\1-Coding\6-CodingDCNTLDiffusionOpt\Codes\1-PreparingDataforMatlab\DNetWeightedL2.txt');
% Totalrows=[InputL1(:,1:2);InputL2(:,1:2)];
% TotalRowUniq=unique(Totalrows(:,1:2),'rows');
% TotalRowUniq(1:size(TotalRowUniq,1),3:4)=zeros(size(TotalRowUniq,1),2);
% for i=1:size(TotalRowUniq,1)
%     for j=1:size(InputL1,1)
%         if TotalRowUniq(i,1)== InputL1(j,1) && TotalRowUniq(i,2)==InputL1(j,2)
%             TotalRowUniq(i,3)=InputL1(j,3); % Relation value from L1
%             break;
%         end
%     end
% end
% for i=1:size(TotalRowUniq,1)
%     for j=1:size(InputL2,1)
%         if TotalRowUniq(i,1)== InputL2(j,1) && TotalRowUniq(i,2)==InputL2(j,2)
%             TotalRowUniq(i,4)=InputL2(j,3); % Relation value from L2
%             break;
%         end
%     end
% end
% save TotalNetworkRelationValue.mat
%%
% =================calculating the source and target relation values=====
clear all
clc
load ('DriverNetLayerOneFinal.mat','ResultDataLayerOne') % clmn3=IDfrstDriver,clmn4=IDScndDriver 
load ('TotalNetworkRelationValue.mat','TotalRowUniq')
%calculating the fuzzy relation value for the source node for the layer one graph nodes
for i=1:size(ResultDataLayerOne,1)
    for j=1:size(TotalRowUniq,1)
        if ResultDataLayerOne(i,3)==TotalRowUniq(j,1) && ResultDataLayerOne(i,4)==TotalRowUniq(j,2)
            ResultDataLayerOne(i,5:6)=TotalRowUniq(j,3:4); % clmn3=IDfrstDriver,clmn4=IDScndDriver and
            %  clmn5=FirstDriver relation value of layer one
            % ,clmn6=FirstDriver relation value of layer two
            break;
        end
    end
end
%calculating the relation value for the target node for the layer one graph nodes
for i=1:size(ResultDataLayerOne,1)
    for j=1:size(TotalRowUniq,1)
        if ResultDataLayerOne(i,3)==TotalRowUniq(j,2) && ResultDataLayerOne(i,4)==TotalRowUniq(j,1)
            ResultDataLayerOne(i,7:8)=TotalRowUniq(j,3:4); % clmn3=IDfrstDriver,clmn4=IDScndDriver and 
            %clmn 7 = Second Driver relation value of layer one
            %clmn 8 = Second Driver relation value of layer two
            break;
        end
    end
end
load ('DriverNetLayerTwoFinalOrigDescNavi.mat','ResultDataLayerTwo')
%calculating the lower and upper bound of fuzzy relation value for the source node for the layer two graph nodes
for i=1:size(ResultDataLayerTwo,1)
    for j=1:size(TotalRowUniq,1)
        if ResultDataLayerTwo(i,3)==TotalRowUniq(j,1) && ResultDataLayerTwo(i,4)==TotalRowUniq(j,2)
            ResultDataLayerTwo(i,5:6)=TotalRowUniq(j,3:4); % clmn3=IDfrstDriver,clmn4=IDScndDriver and
            %  clmn5=FirstDriver relation value of layer one
            % ,clmn6=FirstDriver relation value of layer two
            break;
        end
    end
end
%calculating the lower and upper bound of fuzzy relation value for the target node for the layer two graph nodes
for i=1:size(ResultDataLayerTwo,1)
    for j=1:size(TotalRowUniq,1)
        if ResultDataLayerTwo(i,3)==TotalRowUniq(j,2) && ResultDataLayerTwo(i,4)==TotalRowUniq(j,1)
            ResultDataLayerTwo(i,7:8)=TotalRowUniq(j,3:4); % clmn3=IDfrstDriver,clmn4=IDScndDriver and 
            %clmn 7 = Second Driver relation value of layer one
            %clmn 8 = Second Driver relation value of layer two
            break;
        end
    end
end
%%
%============creating the fuzzy number for node in the layer #1==========================
W_L1=2;
W_L2=1;
A=0;
B=0;
for i=1:size(ResultDataLayerOne,1)
    ResultDataLayerOne(i,9) = W_L1 * ResultDataLayerOne(i,5);
    ResultDataLayerOne(i,10)= W_L2 * ResultDataLayerOne(i,6);
    ResultDataLayerOne(i,11)= W_L1 * ResultDataLayerOne(i,7);
    ResultDataLayerOne(i,12)= W_L2 * ResultDataLayerOne(i,8);
end
FDataLayerOne(1:size(ResultDataLayerOne,1),1:4)=ResultDataLayerOne(1:size(ResultDataLayerOne,1),1:4);
for i=1:size(ResultDataLayerOne,1)
    for j=1:2
        FDataLayerOne(i,4+j)=ResultDataLayerOne(i,8+j)/max( ResultDataLayerOne(i,9), ResultDataLayerOne(i,10));
        FDataLayerOne(i,6+j)=ResultDataLayerOne(i,10+j)/max( ResultDataLayerOne(i,11), ResultDataLayerOne(i,12));
    end
end
%calculating the fuzzy relation value for the source node
FuzzyDataLayerOne(1:size(FDataLayerOne,1),1:4)=FDataLayerOne(1:size(FDataLayerOne,1),1:4);
for i=1:size(FDataLayerOne,1)
    if FDataLayerOne(i,5)>=FDataLayerOne(i,6)
        Max=FDataLayerOne(i,5);
        Min=FDataLayerOne(i,6);
    else
        Max=FDataLayerOne(i,6);
        Min=FDataLayerOne(i,5);
    end
    FuzzyDataLayerOne(i,5)=Min; % Lower Bound of the fuzzy relation value for source node in layer one
    FuzzyDataLayerOne(i,6)=Max; % Upper Bound of the fuzzy relation value for source node  in layer one
end
%calculating the fuzzy relation value for the target node
for i=1:size(FDataLayerOne,1)
    if FDataLayerOne(i,7)>=FDataLayerOne(i,8)
        Max=FDataLayerOne(i,7);
        Min=FDataLayerOne(i,8);
    else
        Max=FDataLayerOne(i,8);
        Min=FDataLayerOne(i,7);
    end
    FuzzyDataLayerOne(i,7)=Min; % Lower Bound of the fuzzy relation value for Target node  in layer one
    FuzzyDataLayerOne(i,8)=Max; % Upper Bound of the fuzzy relation value for Target node  in layer one
end
%%
%============creating the fuzzy number for node in the layer #2==========================
for i=1:size(ResultDataLayerTwo,1)
    ResultDataLayerTwo(i,9) = W_L1 * ResultDataLayerTwo(i,5);
    ResultDataLayerTwo(i,10)= W_L2 * ResultDataLayerTwo(i,6);
    ResultDataLayerTwo(i,11)= W_L1 * ResultDataLayerTwo(i,7);
    ResultDataLayerTwo(i,12)= W_L2 * ResultDataLayerTwo(i,8);
end
FDataLayerTwo(1:size(ResultDataLayerTwo,1),1:4)=ResultDataLayerTwo(1:size(ResultDataLayerTwo,1),1:4);
for i=1:size(ResultDataLayerTwo,1)
    for j=1:2
        FDataLayerTwo(i,4+j)=ResultDataLayerTwo(i,8+j)/max( ResultDataLayerTwo(i,9), ResultDataLayerTwo(i,10));
        FDataLayerTwo(i,6+j)=ResultDataLayerTwo(i,10+j)/max( ResultDataLayerTwo(i,11), ResultDataLayerTwo(i,12));
    end
end
%calculating the fuzzy relation value for the source node in the layer #2
FuzzyDataLayerTwo(1:size(FDataLayerTwo,1),1:4)=FDataLayerTwo(1:size(FDataLayerTwo,1),1:4);
for i=1:size(FDataLayerTwo,1)
    if FDataLayerTwo(i,5)>=FDataLayerTwo(i,6)
        Max=FDataLayerTwo(i,5);
        Min=FDataLayerTwo(i,6);
    else
        Max=FDataLayerTwo(i,6);
        Min=FDataLayerTwo(i,5);
    end
    FuzzyDataLayerTwo(i,5)=Min; % Lower Bound of the fuzzy relation value for source node in the layer #2
    FuzzyDataLayerTwo(i,6)=Max; % Upper Bound of the fuzzy relation value for source node in the layer #2
end
%calculating the fuzzy relation value for the target node
for i=1:size(FDataLayerTwo,1)
    if FDataLayerTwo(i,7)>=FDataLayerTwo(i,8)
        Max=FDataLayerTwo(i,7);
        Min=FDataLayerTwo(i,8);
    else
        Max=FDataLayerTwo(i,8);
        Min=FDataLayerTwo(i,7);
    end
    FuzzyDataLayerTwo(i,7)=Min; % Lower Bound of the fuzzy relation value for Target node for layer # 2
    FuzzyDataLayerTwo(i,8)=Max; % Upper Bound of the fuzzy relation value for Target node for layer # 2
end
save TwoLayersFuzzyRalationValue.mat % 