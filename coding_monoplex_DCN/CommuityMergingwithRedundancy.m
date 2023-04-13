function [Ro_Final_Community]=CommuityMergingwithRedundancy(FinalCommunity_L1,Size_FinalCommunity_L1,FinalCommunity_L2,Size_FinalCommunity_L2,Adja_Mat,Adja_Mat_L2,LayerInputNo,Input,InClm)

A=0;
w=1;
FinalCom_Temp=0;
for i=1:size(FinalCommunity_L1,1)
    for j=1:size(FinalCommunity_L2,1)
        FinalCom_Temp(w,1:(Size_FinalCommunity_L1(i,1)+Size_FinalCommunity_L2(j,1)))=[FinalCommunity_L1(i,1:Size_FinalCommunity_L1(i,1)),FinalCommunity_L2(j,1:Size_FinalCommunity_L2(j,1))];
        w=w+1;
    end
end
for k=1:size(FinalCom_Temp,1)
    FinalCom_Temp_size(k,1)=size(FinalCom_Temp(k,:),2)-sum(FinalCom_Temp(k,:)==0);
end
for m=1:size(FinalCom_Temp,1)
    [ro_Total,ro,pc,pc_2hat]=RedundancyFunc(FinalCom_Temp,FinalCom_Temp_size,Adja_Mat,Adja_Mat_L2,LayerInputNo(1,1),LayerInputNo(1,2),Input,InClm);
    Ro_Final_Community(m,1)=ro_Total;
end
B=0;
