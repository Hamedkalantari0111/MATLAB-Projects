% %%=====================================================================
% %%                          Preparing First Layer Data (Edge= CompanyID)
% %%=====================================================================
% clc
% clear all
% InputData=dlmread('test3.txt');
% [Row,Column]=size(InputData);
% % FIRSTDRIVERID_Orig=InputData(:,1);
% % FIRSTDRIVERID_Dest=InputData(:,2);
% % FRIGHTERID_Orig   =InputData(:,3);
% % FRIGHTERID_Dest   =InputData(:,4);
% % COMPANYID_Orig    =InputData(:,5);
% % COMPANYID_Dest    =InputData(:,6);
% % ORIGINZONEID_Orig =InputData(:,7);
% % ORIGINZONEID_Dest =InputData(:,8);
% % DESTZONEID_Orig   =InputData(:,9);
% % DESTZONEID_Dest   =InputData(:,10);
% save Input_Layer2.mat    % Layer 2 is the Layer of the Drivers Network of which the edge is the "CompanyID"
% clear all
% load ('Input_Layer2.mat','InputData','Row','Column')
% DNet=unique(InputData(:,2:3),'rows');
% A=unique(InputData(:,2),'rows');
% B=unique(InputData(:,3),'rows');
% C=[A(:,1),B(:,1)];
% for i=1: size(A,1)
%     if C(i,1)==C(i,2)
%         C(i,3)=1;
%     else
%         C(i,3)=0;
%     end
% end
% for j=1:2
%     for k=1:size(DNet,1)
%         [row, col]=find(A==DNet(k,j));
%         DNet(k,j+2)=row;  %DNet(:,3) means DNet first node Index and %DNet(:,4) means DNet Second node Index
%     end
% end
% save InputNetwork_layer2.mat   % Layer 2 is the Layer of the Drivers Network of which the edge is the "CompanyID"
% clear all
% load ('InputNetwork_layer2.mat','DNet','A')
% NumofRows=100;
% [row, col]=find(DNet(:,3)==NumofRows);
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
%             if DNetSenario1(i,3)== DNetSenario1(j,4) && DNetSenario1(i,4)==DNetSenario1(j,3)
%                 FF(k,1)=i;
%                 FF(k,2)=j;
%                 k=k+1;
%             end
%         end
%     end
% end
% k=0;
% for i=1:size(DNetSenario1,1)
%     if ismember(i,FF(:,2))==0
%         k=k+1;
%         DNetSenario2(k,:)=DNetSenario1(i,:); 
%     else
%     end
% end
% %csvwrite('MyFile11.txt',DNet_Senario1)
% % xlswrite('MyFile5.xls',DNetSenario2(:,3:4))
% Length= size(CC,1);
% save DNetSenario01_Layer2.mat
clear all
load ('DNetSenario01_Layer2.mat','DNetSenario2', 'Length')
% %==============calculation the adjacency matrix==================
% AdjaMat = zeros(Length,Length);%max(DNet_Senario1(:,3)),max(DNet_Senario1(:,4))); % Adjacency Matrix
% AdjaMat(sub2ind([Length,Length], DNetSenario1(:,3),DNetSenario1(:,4))) = 1;
% %AdjaMat(sub2ind([Length,Length], DNetSenario1(:,4),DNetSenario1(:,3))) = 1;
F=DNetSenario2(:,3);
G=DNetSenario2(:,4);
GG = graph(F,G);
h=plot(GG)
% bg = biograph(AdjaMat);  % make biograph object
% dolayout(bg);   % automatically calculate positions for nodes
% view(bg); % what it says on the tin
% % % for i=1:Length
%     Adja_Mat(DNetSenario1(i,1),DNetSenario1(i,2))= 1;
%     Adja_Mat(DNetSenario1(i,2),DNetSenario1(i,1))= 1; % if Graf is not Directed
% end
save DGraphLayer2.mat
%%=====================================================================
                          %Preparing Second Layer Data (Edge= CompanyID)
%%=====================================================================


