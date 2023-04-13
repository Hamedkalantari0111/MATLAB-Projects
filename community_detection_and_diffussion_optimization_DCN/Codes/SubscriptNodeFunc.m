function[Subscr_node,Subscript_node]= SubscriptNodeFunc(local_Community,Commun_size)
Sub=0;
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        for k=1:Commun_size(i,1)
            for l=1:Commun_size(j,1)
                if i~=j && local_Community(i,k)==local_Community(j,l)
                    Sub=Sub+1;
                end
            end
            Subscr_node(i,j,k)=Sub;
            Sub=0;
        end
        Subscript_node(i,j)=sum(Subscr_node(i,j,:)); % nodes eshteraki beine har 2 local community
    end
end