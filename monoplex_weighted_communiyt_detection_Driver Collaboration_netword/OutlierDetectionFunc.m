function [Outlayer,NodesExistinAllComs]=OutlierDetectionFunc(Community,community_size,Nodes)
k=0;
for i=1: size(Community,1)
    community_size(i,1)=size(Community,2)-sum(Community(i,:)==0); 
    NodesExtincom(k+1:k+community_size(i,1),1)=Community(i,1:community_size(i,1));
    k=k+community_size(i,1);
end
NodesExistinAllComs=unique(NodesExtincom(:,1));
w=1;
for l=1:Nodes
    if ismember(l,NodesExistinAllComs(:,1))==0
        Outlayer(w,1)= l;
        w=w+1;
    else
        Outlayer(w,1)=0;
    end
end
