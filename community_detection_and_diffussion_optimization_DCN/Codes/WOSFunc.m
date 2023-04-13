function [WOS,WOS_p1,WOS_p2]=WOSFunc(local_Community,Commun_size,Comm_loc_edge_index,Subscript_node,Alfa,Subscript_edge)
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        if Commun_size(i,1)<=Commun_size(j,1)
            M_p1_down=Commun_size(i,1);
        else
            M_p1_down=Commun_size(j,1);
        end
        if  Comm_loc_edge_index(i,1)<= Comm_loc_edge_index(j,1)
            M_p2_down= Comm_loc_edge_index(i,1);
        else
            M_p2_down= Comm_loc_edge_index(j,1);
        end
        WOS_p1(i,j)=Alfa*((Subscript_node(i,j))/M_p1_down);
        WOS_p2(i,j)=(1-Alfa)*((Subscript_edge(i,j))/M_p2_down);
        WOS(i,j)= WOS_p1(i,j)+ WOS_p2(i,j); % Weighted Overlapping Score
    end
end