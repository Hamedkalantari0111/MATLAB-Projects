function [Final_Community,Final_Community_size]=FinalCommunity_Phase2(LC,WOS,Beta) % LC=local Communtiy

k_in=1;
l_in=1;
T=1;
tabu_Index=zeros(1,1);   %tabu=communityhaye ghabli ke bayad hazf shavad dar interation badi
lenght=1;
q=1;
for i=1:size(WOS,1)
    if ismember(i,tabu_Index)==0
        Community(q,l_in)=i;          %=========final Communities in phase 2 =======
        l_in=l_in+1;
    for j=1:size(WOS,2)
        if ismember(j,tabu_Index)==0
            if  j>=i && WOS(i,j)>=Beta  %============= inja bayad WOS iterative bashad va static nabashad. chon comnunity jadid update mishvad dar inja
                Community(q,l_in)=j;
                l_in=l_in+1;
            else
                Community(q,l_in)=0;
            end
        else
            Community(q,l_in)=0;
        end
    end
    final_Community_set_size(q,1)=size(Community,2)-sum(Community(q,:)==0);
    for m=1:final_Community_set_size(q,1)
        if Community(q,m)~=0
        tabu_Index(1,lenght)= Community(q,m);       %tabu=communityhaye ghabli ke bayad hazf shavad dar interation badi
        lenght=lenght+1;
        end
    end
    l_in=1;
    q=q+1;
    end
end
w=1;
     %=====final Community=extracting the nodes of the final community sets=========
for i=1:size(Community,1)
    for j=1:final_Community_set_size(i,1)
        for k=1:size(LC,2)
            if LC(Community(i,j),k)~=0 
                Final_Comm_set(i,w)=LC(Community(i,j),k);      % nodes of the final community
                w=w+1;
            end
        end
    end
    w=1;
end
% final community size (nodes in community)
for i=1:size(Final_Comm_set,1)
    Final_Comm_size(i,1)=size(Final_Comm_set,2)-sum(Final_Comm_set(i,:)==0);
end
for i=1:size(Final_Comm_set,1)
    a=unique(Final_Comm_set(i,1:Final_Comm_size(i,1)));
    for j=1:size(Final_Comm_set,2)
        if j<=size(unique(Final_Comm_set(i,1:Final_Comm_size(i,1))),2)
            aa=a(j);
        else
            aa=0;
        end
            Final_Community(i,j)= aa; % Final Comminity nodes of the second phase 
    end
end
for i=1:size(Final_Community,1)
    Final_Community_size(i,1)=size(Final_Community,2)-sum(Final_Community(i,:)==0); % the size of the Final Community
end