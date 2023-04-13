function[Local_edge_index,Comm_loc_edge_index]= LocalEdgeFunc(local_Community,Commun_size,Adja_Mat,Input,InClm,LayerInputNo,F)
No=1;
Local_edge_index=0;
for i=1:size(local_Community,1) % for each community
    for j=1:Commun_size(i,1)% for each node in local community
        for k=1:Commun_size(i,1) % be ezaye har node digar an satr
            if Adja_Mat(local_Community(i,j),local_Community(i,k))==1 % "Adja_Mat" is considered!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                for l=1+F:LayerInputNo+F 
                    if local_Community(i,j)==Input(l,1) && local_Community(i,k)==Input(l,2)
                        Local_edge_index(i,No)=Input(l,InClm+2); % edge index in each local community
                        No=No+1;
                        break
                    end                    
                end
            end
        end
    end
    No=1;
end
for i=1:size(Local_edge_index,1)
    Comm_loc_edge_index(i,1)=size(Local_edge_index,2)-sum(Local_edge_index(i,:)==0); % the number of the edges in each local community
end