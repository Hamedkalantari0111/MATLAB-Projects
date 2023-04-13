function[TwoLayersCommunity_Local,TwoLayersCommunity_Local_Size]= TwoLayers_Intersect_Func(Final_Community,Final_Community_size,Final_Community_L2,Final_Community_size_L2)
for i=1:size(Final_Community_size,1)
    for j=1:size(Final_Community_size_L2,1)
        S=intersect(Final_Community(i,1:Final_Community_size(i,1)),Final_Community_L2(j,1:Final_Community_size_L2(j,1)));
        Len_S=size(S,2);
        TwoLayersCommunity_Local(i,1:Len_S,j)=S(1,:);
        TwoLayersCommunity_Local_Size(i,j)=Len_S;
    end
    S=0;
end