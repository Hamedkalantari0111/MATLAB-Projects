function [OverlappingNodes,OverlappingNodes_Size,OverlappingNodes_proportion,proportion,comm_count_node,comm_proportion_node]=overlappingDetectionFunc(Community,community_size,NodesExistincom)
OL=0;
a=0;
b=0;
OverlappingNodes_proportion=0;
OverlappingNodes_Size=0;
OverlappingNodes=0;
for i=1:size(Community,1)
    k=1;
    for l=1:size(Community,1)
        for j=1:community_size(i,1)
            if ismember(Community(i,j),Community(l,:))==1 && i~=l
                OL(i,k)=Community(i,j);
                k=k+1;
            else
                OL(i,k)=0;
            end
        end
    end
    
end
OverlappingNodes=zeros(size(Community,1),size(OL,2));
for i=1:size(Community,1)
      a=unique(OL(i,:));
      b=sort(a,'desc');
      OverlappingNodes(i,1:size(b,2))=b; % Overlapping Nodes in each Community
      OverlappingNodes_Size(i,1)=size(OverlappingNodes,2)-sum(OverlappingNodes(i,:)==0); % Number of the overlapping nodes in each Community
end
for i=1:size(community_size,1)
      OverlappingNodes_proportion(i,1)=OverlappingNodes_Size(i,1)/community_size(i,1); % Overlapping Nodes in each Community
end
AvgOverlapping=median(community_size(:,1));
maxOverlapping=max(community_size(:,1));
minOverlapping=min(community_size(:,1));
Avgvall=find(community_size(:,1)==AvgOverlapping);
Maxvall=find(community_size(:,1)==maxOverlapping);
minvall=find(community_size(:,1)==minOverlapping);
proportion.Maxval=mean(OverlappingNodes_proportion(Maxvall(:,1),1));
proportion.Avgval=mean(OverlappingNodes_proportion(Avgvall(:,1),1));
proportion.Minval=mean(OverlappingNodes_proportion(minvall(:,1),1));
proportion.Total=mean(OverlappingNodes_proportion(:,1));
comm_count_node=zeros(size(NodesExistincom,1),1);
for i=1:size(NodesExistincom,1)
    for j=1:size(Community,1)
        if ismember(NodesExistincom(i,1),Community(j,1:size(Community,2)-sum(Community(j,:)==0)))==1
            comm_count_node(i,1)=comm_count_node(i,1)+1;
        end
    end
end
for i=1:size(Community,1)
    comm_proportion_node(i,1)=comm_count_node(i,1)/size(Community,1);
end
