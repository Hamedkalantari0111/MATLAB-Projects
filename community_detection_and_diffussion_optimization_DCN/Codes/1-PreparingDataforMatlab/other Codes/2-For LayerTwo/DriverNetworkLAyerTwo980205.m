% %=====================================================================
% %                          Preparing Second Layer Data (Edge= BarnamehID)
% %=====================================================================
% clc
% clear all
% InputData=dlmread('DNetHamsafar.txt');
% [Row,Column]=size(InputData);
% FIRSTDRIVERID_Orig=InputData(:,1);
% FIRSTDRIVERID_Dest=InputData(:,2);
% % FRIGHTERID_Orig   =InputData(:,3);
% % FRIGHTERID_Dest   =InputData(:,4);
% % COMPANYID_Orig    =InputData(:,5);
% % COMPANYID_Dest    =InputData(:,6);
% % ORIGINZONEID_Orig =InputData(:,7);
% % ORIGINZONEID_Dest =InputData(:,8);
% % DESTZONEID_Orig   =InputData(:,9);
% % DESTZONEID_Dest   =InputData(:,10);
% save Input_Layer2Hamsafar.mat % Layer 1 is the Layer of the Drivers Network of which the edge is the "frighterID"
% clear all
% clc
% load ('DriverNetLayerOneFinal.mat','ResultDataLayerOne')
% load ('Input_Layer2Hamsafar.mat','InputData','Row','Column')
% LoneN=[ResultDataLayerOne(:,3);ResultDataLayerOne(:,4)];
% LayerOneNodes=unique(LoneN,'rows');
% k=1;
% for i=1:size(InputData,1)
%     for j=1:size(LayerOneNodes,1)
%         if InputData(i,1)==LayerOneNodes(j,1)
%             for m=1:size(LayerOneNodes,1)
%                 if InputData(i,2)==LayerOneNodes(m,1)
%                     NewInputData(k,1:2)=InputData(i,1:2);
%                     k=k+1;
%                 end
%             end
%         end
%     end
% end
% 
% 
% DNet=unique(NewInputData(:,1:2),'rows');
% A=unique(NewInputData(:,1),'rows');
% B=unique(NewInputData(:,2),'rows');
% % C=[A(:,1),B(:,1)];
% % for i=1: size(A,1)
% %     if C(i,1)==C(i,2)
% %         C(i,3)=1;
% %     else
% %         C(i,3)=0;
% %     end
% % end
% Temp=unique([NewInputData(:,1);NewInputData(:,2)]);
% for j=1:2
%     for k=1:size(DNet,1)
%         [row, col]=find(Temp==DNet(k,j));
%         DNet(k,j+2)=row;  %DNet(:,3) means DNet first node Index and %DNet(:,4) means DNet Second node Index
%     end
% end
% save InputNetwork_Layer2Hamsafar.mat % Layer 1 is the Layer of the Drivers Network of which the edge is the "frighterID"
% clear all
% load ('InputNetwork_Layer2Hamsafar.mat','DNet','A')
% NumofNodes= max(DNet(:,3)); %%% number of the nodes envolved in the Layer 1 %%%
% [row, col]=find(DNet(:,3)==NumofNodes);
% DNet_Senario1= DNet(1:max(row(:,1)),:);
% AA=DNet_Senario1(:,1);
% AA(size(DNet_Senario1,1)+1:2*size(DNet_Senario1,1))=DNet_Senario1(:,2);
% [BB DD]=sort(AA);
% j=1;
% BB(1,2)=1;
% for i=2:size(AA,1)
%     if BB(i,1)>BB(i-1,1)
%         j=j+1;
%         BB(i,2)=j;
%     else
%         BB(i,2)=j;
%     end
% end
% BB(:,3)=DD(:,1);
% CC=unique(BB(:,1:2),'rows');
% for i=1:size(DNet_Senario1,1)
%     DNet_Senario1(i,5)= find(CC(:,1)==DNet_Senario1(i,1));
%     DNet_Senario1(i,6)= find(CC(:,1)==DNet_Senario1(i,2));
% end
% DNetSenario1=[(DNet_Senario1(:,1:2)),(DNet_Senario1(:,5:6))]; % Final MAPPED Inputs for Driver Network in Modeling (MY MAPPED CASE STUDY DATA)!!!!!!!!!!!!!!
% k=1;
% for i=1:size(DNetSenario1,1)
%     for j=1:size(DNetSenario1,1)
%         if i~=j
%             if DNetSenario1(i,3)== DNetSenario1(j,4) && DNetSenario1(i,4)==DNetSenario1(j,3) % find the Duplicate records in the Networc (because of Undircted Graph)
%                 FF(k,1)=i; % the index of the first duplicate record
%                 FF(k,2)=j; % the index of the Second duplicate record
%                 k=k+1;
%             end
%         end
%     end
% end
% FF(:,3)=0;
% for m=1:size(FF,1)
%     [a b]=find(FF(:,2)==FF(m,1));
%     if FF(a,3)==0
%         FF(a,3)=0;
%         FF(m,3)=1;
%     else
%     end   
% end
% k=0;
% for i=1:size(DNetSenario1,1)
%     [A B]=ismember(i,FF(:,1));
%     if A==0
%         k=k+1;
%         DNetSenario2(k,:)=DNetSenario1(i,:); 
%     elseif A==1 && FF(B,3)==1
%         k=k+1;
%         DNetSenario2(k,:)=DNetSenario1(i,:);  %Final Result
%     end
% end
% %csvwrite('MyFile11.txt',DNet_Senario1)
% xlswrite('layerTwoMatlabPrepared.xls',DNetSenario2(:,3:4)) %clmn1=IDfrstDriver,clmn2=IDScndDriver,clmn3=SourceCode(FirstDriverMappedID),clmn4=TargetCode(SecondDriverMappedID)
% Length= size(CC,1);
% save DNetlayerTwoResult.mat
% clear all
% load ('DNetlayerTwoResult.mat','DNetSenario2', 'Length')
% % % %==============calculation the adjacency matrix==================
% % % AdjaMat = zeros(Length,Length);%max(DNet_Senario1(:,3)),max(DNet_Senario1(:,4))); % Adjacency Matrix
% % % AdjaMat(sub2ind([Length,Length], DNetSenario1(:,3),DNetSenario1(:,4))) = 1;
% % % %AdjaMat(sub2ind([Length,Length], DNetSenario1(:,4),DNetSenario1(:,3))) = 1;
% % F=DNetSenario2(:,3);
% % G=DNetSenario2(:,4);
% % GG = graph(F,G);
% % h=plot(GG)
% % % bg = biograph(AdjaMat);  % make biograph object
% % % dolayout(bg);   % automatically calculate positions for nodes
% % % view(bg); % what it says on the tin
% % % % % for i=1:Length
% % %     Adja_Mat(DNetSenario1(i,1),DNetSenario1(i,2))= 1;
% % %     Adja_Mat(DNetSenario1(i,2),DNetSenario1(i,1))= 1; % if Graf is not Directed
% % % end
% % save DGraphLayer1.mat % LAYER1=DNetwork on same Navy (Edge=FreighterID or Plaque)
% %%=====================================================================
% k=1;
% ResultDataLayerTwo=dlmread('ResultFromGephiL2.txt');
% for i=1:size(ResultDataLayerTwo,1)
%     for j=1:size(DNetSenario2,1)
%         if ResultDataLayerTwo(i,1)==DNetSenario2(j,3) && ResultDataLayerTwo(i,2)==DNetSenario2(j,4)
%             ResultDataLayerTwo(i,3:4)=DNetSenario2(j,1:2); % clmn3=IDfrstDriver,clmn4=IDScndDriver
%         end
%     end
% end
% save DriverNetLayerTwoFinal.mat % 
%%=====================================================================
clear all
clc
load ('DriverNetLayerTwoFinal.mat','ResultDataLayerTwo')
load ('TotalNetFuzzyRelationValue.mat','TotalRowUniq')
for i=1:size(ResultDataLayerTwo,1)
    for j=1:size(TotalRowUniq,1)
        if ResultDataLayerTwo(i,3)==TotalRowUniq(j,1) && ResultDataLayerTwo(i,4)==TotalRowUniq(j,2)
            ResultDataLayerTwo(i,5:6)=TotalRowUniq(j,3:4); % clmn3=IDfrstDriver,clmn4=IDScndDriver 
            % and clmn5=FirstDriver lower bound of fuzzy relation value ,clmn6=FirstDriver upper bound of fuzzy relation value
            break;
        end
    end
end
save DriverNetL2withFuzzyRelationValue.mat % 
