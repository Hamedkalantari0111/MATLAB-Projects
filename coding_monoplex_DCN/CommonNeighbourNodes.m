function [NC]=CommonNeighbourNodes(w,N,NC)
[dd ee]=ismember (N(1,1:size(N,2)-sum(N(1,:,w)==0),w),N(2,1:size(N,2)-sum(N(2,:,w)==0),w));
[AA BB]=sort(dd(1,:),'desc');
CC(1:size(AA,2)-sum(AA(1,:)==0),1)=N(1,BB(1:size(AA,2)-sum(AA(1,:)==0)),w);
if size(CC,1)>0
    NC(1,1:size(CC,1)-sum(CC(:,1)==0))=CC(1:size(CC,1)-sum(CC(:,1)==0),1);
else
    NC=0;
end
