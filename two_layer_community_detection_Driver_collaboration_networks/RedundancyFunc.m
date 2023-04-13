function [ro_Total,ro,pc,pc_2hat,proportion]=RedundancyFunc(FinalCommunity_phase4,FinalCommunity_size_phase4,Adja_Mat,Adja_Mat_L2,LayerInputNo1,LayerInputNo2,Input,InClm)
F=0;
[Loc_edge_index_Layer1,Comm_loc_edge_index_Layer1]= LocalEdgeFunc(FinalCommunity_phase4,FinalCommunity_size_phase4,Adja_Mat,Input,InClm,LayerInputNo1,F); 
F=LayerInputNo1;
[Loc_edge_index_Layer2,Comm_loc_edge_index_Layer2]= LocalEdgeFunc(FinalCommunity_phase4,FinalCommunity_size_phase4,Adja_Mat_L2,Input,InClm,LayerInputNo2,F); 
count=0;
pc_2hat=0;
for i=1:size(Loc_edge_index_Layer1,1)
    for k=1:Comm_loc_edge_index_Layer1(i,1)
        for m=1:Comm_loc_edge_index_Layer2(i,1)
            if Input(Loc_edge_index_Layer1(i,k),6)==Input(LayerInputNo1+Loc_edge_index_Layer2(i,m),6) && Input(Loc_edge_index_Layer1(i,k),6) ~=0
                count=count+1;
                pc_2hat(i,1)=count; % yalhaye com c ke dar har 2 layer mojod mibashand
            else
                pc_2hat(i,1)=count; % yalhaye com c ke dar har 2 layer mojod mibashand
            end
        end
    end
    count=0;
end
for i=1:size(Loc_edge_index_Layer1,1)
    pc(i,1)=Comm_loc_edge_index_Layer1(i,1)+Comm_loc_edge_index_Layer2(i,1); % yalhaye com c ke hade aghal dar 1 layer mojod mibashand
end
Non_count=0;
for i=1:size(pc,1)
    if  pc_2hat(i,1) ~=0
        ro(i,1)= (2*pc_2hat(i,1))/(2*pc(i,1));
    else
        ro(i,1)=0;
        Non_count=Non_count+1;
    end
end
ro_Total=(1/(size(pc,1)-Non_count))*sum(ro(:,1))


Avgro=median(FinalCommunity_size_phase4(:,1));
maxro=max(FinalCommunity_size_phase4(:,1));
minro=min(FinalCommunity_size_phase4(:,1));
Avgvall=find(FinalCommunity_size_phase4(:,1)==Avgro);
Maxvall=find(FinalCommunity_size_phase4(:,1)==maxro);
minvall=find(FinalCommunity_size_phase4(:,1)==minro);
proportion.Maxval=mean(ro(Maxvall(:,1),1));
proportion.Avgval=mean(ro(Avgvall(:,1),1));
proportion.Minval=mean(ro(minvall(:,1),1));
proportion.Total=mean(ro(:,1));