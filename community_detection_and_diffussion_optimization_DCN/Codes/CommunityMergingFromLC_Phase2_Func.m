function [Final_Community_merged,Final_Community_merged_size]=CommunityMergingFromLC_Phase2_Func(LC,q)
w=1;
%%=====extracting the nodes of the CommunityMerged sets=========
for m=1:size(LC,1)
    for k=1:size(LC,2)
        if LC(m,k)~=0 
            Community_merged(q,w)=LC(m,k);      % nodes of the final community
            w=w+1;
        end
    end
end
w=1;
% %  communityMerger size (nodes in communitymerged)
Community_merged_size(q,1)=size(Community_merged,2)-sum(Community_merged(q,:)==0);
a=unique(Community_merged(q,1:Community_merged_size(q,1)));
for n=1:size(Community_merged,2)
    if n<=size(unique(Community_merged(q,1:Community_merged_size(q,1))),2)
        aa=a(n);
    else
        aa=0;
    end
    Final_Community_merged(q,n)= aa; % Final ComminityMerged nodes of the second phase
end
Final_Community_merged_size=size(Final_Community_merged,2)-sum(Final_Community_merged(q,:)==0); % the size of the Final Community Merged