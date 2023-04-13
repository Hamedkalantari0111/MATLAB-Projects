function[nodes_not_in_com]= NodesnotinCommsFunc_phase3(Final_Community,Final_Community_size,Nodes)
w=1;
for i=1:size(Final_Community,1)
    for j=1:Final_Community_size(i,1)
        FCN(w,1)=Final_Community(i,j); %Final_community_Nodes  not unique
        w=w+1;
    end
end
Final_community_Nodes=unique(FCN,'sorted'); %Final_community_Nodes unique

[a bb]=ismember(1:Nodes,Final_community_Nodes);
[c dd]=sort(a,'descend');
ee=sum(c(1,:)==0);
ff=size(c,2)-sum(c(1,:)==0);
if ee==0
    nodes_not_in_com(1,1)=0; % be jaye k man 1 gozashtam!!!!
else
    for k=1:ee
        nodes_not_in_com(k,1)=dd(1,k+ff); % Nodes not in the communities
    end
end