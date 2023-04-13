function [edge_Sorted,CC]= NodeProbing_Func(a,Adja_Mat)
x=1;
tab=0;
newlist=a;
for i=1:size(a,1)
    dd=0;
    ee=0;
    v=0;
    if size(newlist)>0
        [AA BB]=sort(Adja_Mat(newlist(1,1),:),'desc');
        CC(1:size(AA,2)-sum(AA(1,:)==0),i)=BB(1,1:size(AA,2)-sum(AA(1,:)==0));
        [dd ee]=ismember (newlist(:,1),CC(:,i));
    v = nonzeros(ee);
    for k=1:size(v,1)
        edge_Sorted(x,1:2)=[newlist(1,1),CC(v(k,1),i)];
        x=x+1;
    end
    tab=[tab;newlist(1,1);edge_Sorted(:,2)];
    tabue=nonzeros(tab);
    flag=~ismember(a(:,1),tabue(:,1));
    index=find(flag);
    flag=0;
    newlist=0;
    newlist=a(index,:);
    index=0;
    end
end






