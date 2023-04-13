function [N]=NeighbourNodes(w,Adja_Mat,LC)
AA=0;
BB=0;
CC=0;
for i=1:size(LC,2)
    [AA BB]=sort(Adja_Mat(LC(1,i),:),'desc');
    CC(i,1:size(AA,2)-sum(AA(1,:)==0))=BB(1,1:size(AA,2)-sum(AA(1,:)==0));
    N(i,1:size(CC,2)-sum(CC(i,:)==0),w)=CC(i,1:size(CC,2)-sum(CC(i,:)==0));
    AA=0;
    BB=0;
end
