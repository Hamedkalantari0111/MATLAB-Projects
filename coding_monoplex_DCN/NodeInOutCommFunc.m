function[NodesInCommunity,NodesOutCommunity]= NodeInOutCommFunc(LC,Nodes,NeighNodesLocalCom) % LC=local_Community
 Delta_temp=zeros(Nodes,Nodes); %if both two nodes i,j were in the Community =1 , O.W=0;
 landa_temp=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
 CC_temp=zeros(size(LC,1),size(NeighNodesLocalCom,2));
 DD_temp=zeros(size(LC,1),size(NeighNodesLocalCom,2));
 for i=1:size(LC,1)
     [C,D]=ismember(NeighNodesLocalCom(i,:),LC(i,1:size(LC,2)-sum(LC(i,:)==0)));
     CC_temp(i,:)=C;
     DD_temp(i,:)=D; % the index of the local_Community where the node of the "NeighNodesLocalCom" is there  
 end
 w=0;
 y=0;
 for k=1:size(DD_temp,1)
     for l=1:size(DD_temp,2)
         w=w+1;
         y=y+1;
         if DD_temp(k,l)==0 && NeighNodesLocalCom(k,l)~=0  
             NodesOutCommunity1(k,w)=NeighNodesLocalCom(k,l);
              NodesInCommunity1(k,y)=0;
         else if DD_temp(k,l)~=0 && NeighNodesLocalCom(k,l)~=0
                 NodesInCommunity1(k,y)=NeighNodesLocalCom(k,l);
                 NodesOutCommunity1(k,w)=0;
             end
         end
     end
     w=0;
     y=0;
     xx=unique(NodesOutCommunity1(k,:));
     xxx=sort(xx,'Desc');
     NodesOutCommunity(k,1:size(xxx,2)-sum(xxx(1,:)==0))=xxx(1,1:size(xxx,2)-sum(xxx(1,:)==0));
     yy=unique(NodesInCommunity1(k,:));
     yyy=sort(yy,'Desc');
     NodesInCommunity(k,1:size(yyy,2)-sum(yyy(1,:)==0))=yyy(1,1:size(yyy,2)-sum(yyy(1,:)==0));
     xx=0;
     yy=0;
 end