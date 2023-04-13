function [final_Community_set_final,finall_Community_set_size,Community]=FinalCommunity_Phase2_WOS_NCA_Itrtv(LC,WOS,Beta,AdjaMat,Input,InClm,LayerInputNo,F,Nodes,Alfa,Gama) % LC=local Communtiy

k_in=1;
l_in=1;
T=1;
tabu_Index=zeros(1,1);   %tabu=communityhaye ghabli ke bayad hazf shavad dar interation badi
lenght=1;
q=1;
WOS_orig=WOS;
for i=1:size(LC,1) % LC??????????
    LC_New=0;
    if ismember(i,tabu_Index)==0
        Community(q,l_in)=i;          %=========final Communities in phase 2 =======
        Communityy(q,1)=i;
        l_in=l_in+1;
        [a b]=sort(WOS_orig(i,:),'Desc');
        for y=1:size(a,2)
            if a(1,y)>Beta;
                WOS_Rank(i,y)=b(1,y); % bare avval baraye har comm hesab minokin baraye marhale bad
            else
              if y==1
                 WOS_Rank(i,y)=0;
              else
                break
              end
            end
        end
        LC_New=LC(i,:);
    if WOS_Rank(i,1)~=0
    for j=1:size(WOS_Rank(i,1:size(WOS_Rank,2)-sum(WOS_Rank(i,:)==0)),2)
        if ismember(WOS_Rank(i,j),tabu_Index)==0
            if j>1
                LC_New=0;
                LC_New(1,1:size(final_Community_set_final,2)-sum(final_Community_set_final(q,:)==0))=final_Community_set_final(q,1:size(final_Community_set_final,2)-sum(final_Community_set_final(q,:)==0));
                LC_New(2,1:size(LC(WOS_Rank(i,j),:),2))=LC(WOS_Rank(i,j),:);
                [WOS_New,NCA,NCA_p1,NCA_p2]=WOS_NCA_Function(LC_New,AdjaMat,Input,InClm,LayerInputNo,F,Nodes,Alfa,Gama); % my new WOS function!!!!!!!!!
                WOS_orig(i,WOS_Rank(i,j))=WOS_New(1,2);
            end
            if  WOS_Rank(i,j)>=i && WOS_orig(i,WOS_Rank(i,j))>=Beta;  %============= inja bayad WOS iterative bashad va static nabashad. chon comnunity jadid update mishvad dar inja
                Community(q,l_in)=WOS_Rank(i,j);
                Communityy(q,2)=WOS_Rank(i,j);
                LC_New(2,1:size(LC(WOS_Rank(i,j),:),2))=LC(WOS_Rank(i,j),:);
                [Community_merged,Community_merged_size]=CommunityMergingFromLC_Phase2_Func(LC_New,q); % updating the community q (community(q))%and calculating "WOS_NCA"
                final_Community_set(1,1:size(Community_merged,2)-sum(Community_merged(q,:)==0))=Community_merged(q,1:size(Community_merged,2)-sum(Community_merged(q,:)==0));
                final_Community_set_final(q,1:size(final_Community_set,2)-sum(final_Community_set(1,:)==0))=final_Community_set;
                LC_New(1,1:size(final_Community_set_final,2)-sum(final_Community_set_final(q,:)==0))=final_Community_set_final(q,1:size(final_Community_set_final,2)-sum(final_Community_set_final(q,:)==0));
                final_Community_set=0;
                l_in=l_in+1;
            else
                Community(q,l_in)=0;
                final_Community_set=LC_New(1,1:size(LC_New,2)-sum(LC_New(1,:)==0));
                final_Community_set_final(q,1:size(final_Community_set,2)-sum(final_Community_set(1,:)==0))=final_Community_set;
            end
        else
            Community(q,l_in)=0;
            final_Community_set=LC_New(1,1:size(LC_New,2)-sum(LC_New(1,:)==0));
            final_Community_set_final(q,1:size(final_Community_set,2)-sum(final_Community_set(1,:)==0))=final_Community_set;
        end
    end
      final_Community_set_size(q,1)=size(Community,2)-sum(Community(q,:)==0);
      finall_Community_set_size(q,1)=size(final_Community_set_final,2)-sum(final_Community_set_final(q,:)==0);
    else
      final_Community_set_final(q,1:size(LC,2)-sum(LC(i,:)==0))=LC(i,1:size(LC,2)-sum(LC(i,:)==0));
      final_Community_set_size(q,1)=size(Community,2)-sum(Community(q,:)==0);
      finall_Community_set_size(q,1)=size(final_Community_set_final,2)-sum(final_Community_set_final(q,:)==0);
    end
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