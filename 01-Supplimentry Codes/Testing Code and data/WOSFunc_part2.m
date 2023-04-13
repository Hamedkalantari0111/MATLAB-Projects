function [WOS_NCA,NCA,NCA_p1,NCA_p2]=WOSFunc_part2(local_Community,Commun_size,Comm_loc_edge_index,Subscript_node,Alfa,Subscript_edge,Gama,fi,Nodes_NotinCom_i_incom_j,Adja_Mat,NodesInCommunity,CCC)

  %================== part 1 of WOS_NCA "just WOS"========================
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        if Commun_size(i,1)<=Commun_size(j,1)
            M_p1_down(i,j)=Commun_size(i,1);
        else
            M_p1_down(i,j)=Commun_size(j,1);
        end
        if  Comm_loc_edge_index(i,1)<= Comm_loc_edge_index(j,1)
            M_p2_down(i,j)= Comm_loc_edge_index(i,1);
        else
            M_p2_down(i,j)= Comm_loc_edge_index(j,1);
        end
        WOS_p1(i,j)=Alfa*((Subscript_node(i,j))/M_p1_down(i,j));
        WOS_p2(i,j)=(1-Alfa)*((Subscript_edge(i,j))/M_p2_down(i,j));
        WOS(i,j)= WOS_p1(i,j)+ WOS_p2(i,j); % Weighted Overlapping Score
    end
end
%===========================part 2 of WOS_NCA means just "NCA"=============
for i=1:size(local_Community,1)
    for j=1:size(local_Community,1)
        if i~=j
            NCA_P1_up=Nodes_NotinCom_i_incom_j(i,j)+Nodes_NotinCom_i_incom_j(j,i);
            if Commun_size(i,1)<=Commun_size(j,1)
                NCA_p1_down1(i,j)=Commun_size(i,1);
            else
                NCA_p1_down1(i,j)=Commun_size(j,1);
            end
            NCA_p1_down=NCA_p1_down1(i,j)+Nodes_NotinCom_i_incom_j(i,j);
            NCA_p1(i,j)=NCA_P1_up/NCA_p1_down;
        end
    end
end
  AAAA=0;
  for i=1:size(NodesInCommunity,1)
      for j=1:size(NodesInCommunity,1)
          if i~=j
              for k=1:size(NodesInCommunity,2)-sum(NodesInCommunity(i,:)==0)
                  for o=1:size(NodesInCommunity,2)-sum(NodesInCommunity(j,:)==0)
                      if  CCC(i,k,j)==0 || CCC(j,o,i)==0
                          AAAA=AAAA+Adja_Mat(NodesInCommunity(i,k),NodesInCommunity(j,o));
                      end
                  end
              end
          end
          NCA_p2_up(i,j)= AAAA;
          NCA_p2_down(i,j)=M_p2_down(i,j)+NCA_p2_up(i,j);
          NCA_p2(i,j)= NCA_p2_up(i,j)/NCA_p2_down(i,j);
          AAAA=0;
      end
  end
  
  for i=1:size(NodesInCommunity,1)
      for j=1:size(NodesInCommunity,1)
          NCA(i,j)=(NCA_p1(i,j)+ NCA_p2(i,j));
          WOS_NCA(i,j)=Gama*WOS(i,j)+(1-Gama).^2* NCA(i,j);
      end
  end
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  