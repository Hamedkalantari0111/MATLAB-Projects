function[Subscr_edge,Subscript_edge]= SubscriptEdgesFunc(local_Community,Comm_loc_edge_index,Loc_edge_index)
Sub=0;
Subscr_edge=zeros(size(local_Community,1),size(local_Community,1));
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        for k=1:Comm_loc_edge_index(i,1)
            for l=1:Comm_loc_edge_index(j,1) 
                if i~=j && Loc_edge_index(i,k)==Loc_edge_index(j,l)
                    Sub=Sub+1;
                end
            end
            Subscr_edge(i,j)=Subscr_edge(i,j)+Sub;
            Sub=0;
        end
        Subscript_edge(i,j)=Subscr_edge(i,j);% tedade yalhaye eshteraki beine har 2 local community
    end
end