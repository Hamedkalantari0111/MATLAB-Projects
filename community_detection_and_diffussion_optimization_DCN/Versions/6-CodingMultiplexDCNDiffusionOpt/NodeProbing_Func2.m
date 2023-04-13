function [edge_Sorted,CC]= NodeProbing_Func2(a,Adja_Mat)

x=1;
edge_Sorted=zeros(1,1);
for i=1:size(a,1) %  D: Could be D3 or Db (for both two layers!)
    for j=1: size(a,1)
        Y=[a(i,1),a(j,1)];
        for u=1: size(edge_Sorted,1)
            [EE,OO]=ismember (Y(1,:),edge_Sorted(u,:)); % for omit the doublicate edges!!!! in edge sort!
            if sum(EE(1,:))==2
                EEE=2;
            break
            end
        end
        if Adja_Mat(a(i,1),a(j,1))==1 && sum(EE(1,:))~=2 
            edge_Sorted(x,1:2)=[a(i,1),a(j,1)]; % edges sorted with respect to the aboved "Centrality measure"!!!!
            x=x+1;
        end
    end
end
% % edge_Sorted=zeros(1,1:2);
% % %=======================================
% % for i=1:size(a,1)
% %     [AA BB]=sort(Adja_Mat(a(i,1),:),'desc');
% %     CC(1:size(AA,2)-sum(AA(1,:)==0),i)=BB(1,1:size(AA,2)-sum(AA(1,:)==0));
% %     [dd ee]=ismember (a(:,1),CC(:,i));
% %     v = nonzeros(ee);
% %     for k=1:size(v,1)
% %         edge_Sortedd(x,2:3)=[a(i,1),CC(v(k,1),i)];
% %         x=x+1;
% %     end
% % end
% % edge_Sortedd(:,1)=1:size(edge_Sortedd,1);
% % edge_Sortedd(:,4)=0;
% % edge_Sortedd(:,5)=0;
% % for m=1:size(edge_Sortedd,1)
% %     a=0;
% %     b=0;
% %     [a b]=find(edge_Sortedd(:,3)==edge_Sortedd(m,2));
% %     for n=1:size(a,1)
% %         if edge_Sortedd(a(n,1),4)==0 && edge_Sortedd(a(n,1),5)==0
% %             edge_Sortedd(a(n,1),4)=0;
% %             edge_Sortedd(m,4)=1;
% %             edge_Sortedd(m,5)=1;
% %         end   
% %         if (ismember(edge_Sortedd(a(n,1),2),edge_Sortedd(m,3))==1) && edge_Sortedd(a(n,1),4)==0
% %             edge_Sortedd(a(n,1),5)=3;
% %         end
% %     end
% % end
% % [aa bb]=find(edge_Sortedd(1:size(edge_Sortedd,1),4)==1)
% % EdgeSort=edge_Sortedd(aa(:,1),:);
% % [ff gg]=find(edge_Sortedd(1:size(edge_Sortedd,1),4)==0)
% % for l=1:size(ff,1)
% %     if ismember(edge_Sortedd(ff(l,1),2),edge_Sortedd(aa(:,1),3))==0
% %         EdgeSort(size(EdgeSort,1)+1,:)=edge_Sortedd(ff(l,1),:);
% %     end
% % end
% tab=0;
% newlist=a;
% for i=1:size(a,1)
%     dd=0;
%     ee=0;
%     v=0;
%     if size(newlist)>0
%         [AA BB]=sort(Adja_Mat(newlist(1,1),:),'desc');
%         CC(1:size(AA,2)-sum(AA(1,:)==0),i)=BB(1,1:size(AA,2)-sum(AA(1,:)==0));
%         [dd ee]=ismember (newlist(:,1),CC(:,i));
%     v = nonzeros(ee);
%     for k=1:size(v,1)
%         edge_Sorted(x,1:2)=[newlist(1,1),CC(v(k,1),i)];
%         x=x+1;
%     end
%     tab=[tab;newlist(1,1);edge_Sorted(:,2)];
%     tabue=nonzeros(tab);
%     flag=~ismember(a(:,1),tabue(:,1));
%     index=find(flag);
%     flag=0;
%     newlist=0;
%     newlist=a(index,:);
%     index=0;
%     end
% end






