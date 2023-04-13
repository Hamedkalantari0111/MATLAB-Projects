clear all
clc
input=[1,2,4,7,12,14,25,27,30,33,35,40];
Cluster=3;
M=2;
K=1;
N=1000;
E=0.01; % Epsilon diff Value
Cores = zeros(N+1,Cluster);
Distance=zeros(N,Cluster,size(input,2)); %Distance of each nodes from each core
D=zeros(N+1,Cluster,size(input,2)); %D measure
Utility=zeros(N+1,Cluster,size(input,2)); %membership Function of each nodes to each core
Diff=zeros(N+1,Cluster); %Diff od utility each nodes to each core
noo=zeros(N+1,Cluster); %meghdare "noo" baraye har core
% eejade marakeze avalie
for k=1:Cluster
    Cores(1,k)= input(k*floor(size(input,2)/Cluster));
end

%mohasebe Utility va Distance and D=(norme oghlidosi fasele har node az core)
for o=1:size(input,2)
    for k=1:Cluster
    Distance(1,k,o)= abs(input(1,o)-Cores(1,k));
    D(1,k,o)= sqrt((input(1,o)-Cores(1,k))^2);
    end
    for n=1:Cluster
    Utility(1,n,o)=1-(Distance(1,n,o)/sum(Distance(1,:,o)));
    end
end

% %Algorithm l=Itration, i=Core, j=Nodes, 
Cores_sub=0;
Cores_down=0;
noo_sub=0;
noo_down=0;
for l=1:N
  for i=1:Cluster
      for j=1:size(input,2)
          noo_sub    = noo_sub    + (Utility(l,i,j)^M)*D(l,i,j);
          noo_down   = noo_down   + (Utility(l,i,j)^M);
          Cores_sub  = Cores_sub  + ((Utility(l,i,j)^M)*input(1,j));
          Cores_down = Cores_down + (Utility(l,i,j)^M);
      end
      Cores(l+1,i)=Cores_sub/Cores_down;
      noo(l,i)=K*(noo_sub/noo_down);
      Cores_sub=0;
      Cores_down=0;
      noo_sub=0;
      noo_down=0;
  end
  for j=1:size(input,2)
      for i=1:Cluster
          if input(1,j)==Cores(l+1,i);
              D(l+1,i,j)=0;
              Utility(l+1,i,j)=1;
          else
              D(l+1,i,j)= sqrt((input(1,j)-Cores(l+1,i))^2);
              Utility(l+1,i,j)=1/(1+((D(l+1,i,j))/(noo(l,i)))^(1/(M-1)));
              Diff(l+1,i)= Diff(l+1,i)+sqrt((Utility(l+1,i,j)-Utility(l,i,j)).^2);
          end
      end
  end
  if l>=3 && (sum(Diff(l+1,:))/(size(Diff,2)))<= E 
      break
  end
end
% multi-center type2 pcm 
% V= cor emembers
noo_sub=0;
noo_down=0;
 for i=1:Cluster
    for j=1:size(input,2)
       noo_sub= noo_sub+(Utility(l+1,i,j)^M)*D(l+1,i,j);
       noo_down= noo_down+(Utility(l+1,i,j)^M);
    end
       noo(l+1,i)=K*(noo_sub/noo_down);
       q=1;
       for  j=1:size(input,2)
        if D(l+1,i,j)<=noo(l+1,i)
            Core_member(i,q,1)=input(1,j);
            q=q+1;
        end
       end
       C(i,1)=q-1;
 end
Itration=1000;
final_Diff=zeros(Itration,Cluster);
VV=zeros(Itration,Cluster,size(C,1),size(C,1));
Dist=zeros(Itration,Cluster,size(C,1),size(C,1));
utility_hat=zeros(Itration,Cluster,size(C,1),size(C,1));
Utility_hathat=zeros(Itration,Cluster,size(C,1));
landa=zeros(Itration,Cluster,size(C,1),size(C,1));
% Core_member=zeros(Itration,Cluster,size(C,1));
beta_up=0;
beta_down=0;
beta=zeros(Cluster,1);
Disttdown=zeros(Itration,Cluster,size(C,1),size(C,1));
Disttup=zeros(Itration,Cluster,size(C,1),size(C,1));
utility_hathat_normalized=zeros(Itration,Cluster,size(C,1));
norm_down=zeros(Itration,Cluster);
part1_up=zeros(Itration,Cluster,size(C,1));
part2_up=zeros(Itration,Cluster,size(C,1));
part1_down=zeros(Itration,Cluster,size(C,1));
part3_down=zeros(Itration,Cluster,size(C,1));
% Core_member=zeros(Itration,Cluster,size(C,1));
 Nooo=1; %Nooo is the fraction of beta
 EE=0.002;
 EU=0.002;
 BIG=100;
 sum=0;
 CC=0;
 FF=zeros(Itration,1);
 sum_up=0;
 for t=1:Itration
     for i=1:Cluster
         for j=1:C(i,1)
             for jprim=1:C(i,1)
                 if j>jprim || j<jprim
                     VV(t,i,j,jprim)=sqrt((Core_member(i,j,t)-Core_member(i,jprim,t)).^2);
                     xx=Nooo*noo(l+1,i)/VV(t,i,j,jprim);
                     if xx==0
                         landa(t,i,j,jprim)=Big;
                         landa(t,i,j,jprim)=xx;
                     end
                 end
             end
         end
     end
     % mohasebe U^ %lower utility_hat_normalized
     for i=1:Cluster
         for j=1:C(i,1)
             for k=1:size(input,2)
                 Dist(t,i,j,k)= sqrt((input(1,k)-Core_member(i,j,t)).^2); % Dist=Distance new
                 utility_hat(t,i,j,k)=1/(1+(Dist(t,i,j,k)/noo(l+1,i)).^(1/(M-1)));
             end
         end
     end

     % mohasebe U^^ %upper =utility_hathat_normalized
     Distt_down=0;
     Distt_up=0;
     sub=0;
     for i=1:Cluster
         for j=1:C(i,1)
             for k=1:C(i,1)
                 Distt_up= Distt_up+ sqrt((Core_member(i,j,t)-Core_member(i,k,t)).^2); 
                 for p=1:C(i,1)
                     Distt_down=Distt_down+ sqrt((Core_member(i,p,t)-Core_member(i,k,t)).^2); 
                 end
                 Disttdown(t,i,j,k)=Distt_down;
                 Distt_down=0;
             end
             Disttup(t,i,j)= Distt_up;
             for k=1:C(i,1)
                 sub=sub+Disttup(t,i,j)/Disttdown(t,i,j,k);
             end
             Utility_hathat(t,i,j)=(sub)^(-1);
             sub=0;
             Distt_up=0;
             CC=CC+ Utility_hathat(t,i,j);
         end
         norm_down(t,i)=CC;
         CC=0;
         for j=1:C(i,1)
             utility_hathat_normalized(t,i,j)= Utility_hathat(t,i,j)/norm_down(t,i);
         end
     end
     %Calculate Core Membera of the cluster
     AA=0;
     BB=0;
     CC=0;
     DD=0;
     for i=1:Cluster
         for j=1:C(i,1)
             for k=1:size(input,2)
                 AA= AA + ((utility_hat(t,i,j,k)).^M)*input(1,k);
                 BB= BB + ((utility_hat(t,i,j,k)).^M);
             end
             part1_up(t,i,j)=AA;
             AA=0;
             part1_down(t,i,j)=BB;
             BB=0;
             for l=1:C(i,1)
                 CC= CC +(((utility_hathat_normalized(t,i,j).^M) + (landa(t,i,j,l).^M))*Core_member(i,l,t));
                 DD= DD + (landa(t,i,j,l).^M);
             end
             part2_up(t,i,j)=CC;
             CC=0;
             part3_down(t,i,j)=DD;
             DD=0;
             Core_member(i,j,t+1)=(part1_up(t,i,j)+part2_up(t,i,j))/(part1_down(t,i,j)+(C(i,1)*(utility_hathat_normalized(t,i,j).^M))+part3_down(t,i,j))
         end
             for j=1:C(i,1)
                 for k=1:size(input,2)
                     if t==1
                          final_Diff(t,i)=2*EU;
                     else
                         final_Diff(t,i)= final_Diff(t,i)+sqrt(((utility_hat(t,i,j,k)-utility_hat(t-1,i,j,k)).^2)+((utility_hathat_normalized(t,i,j)-utility_hathat_normalized(t-1,i,j)).^2));
                     end
                 end
             end
         FF(t) = FF(t) + final_Diff(t,i);
     end
     if t>=3 && FF(t)/(size(final_Diff,2))<= EU
         break
     end
 end
 
     % mohasebe beta new
 for i=1:Cluster
     for k=1:size(input,2)
         for j=1: C(i,1)
             if input(1,k)==Core_member(i,j,t);
                 Distance_new(i,j,k)=0;
             else
                 Distance_new(i,j,k)= sqrt((input(1,k)-Core_member(i,j,t)).^2); %------------- t ya t+1 ???????????????????? baraye core member
             end
             beta_up   = beta_up   +((utility_hat(t,i,j,k).^M)* Distance_new(i,j,k));
             beta_down = beta_down + (utility_hat(t,i,j,k));
         end
     end
     beta(i)=beta_up/beta_down;
 end
 for i=1:Cluster
     for j=1:C(i,1)
         sum_up= sum_up + Core_member(i,j,t);
     end
     Core_member_new(i)=sum_up/C(i,1);
     sum_up=0;
 end
 
% V= cor emembers
% noo_sub=0;
% noo_down=0;
% for i=1:Cluster
%     q=1;
%     for j=1:size(input,2)
%         if D(l+1,i,j)<=noo(l+1,i)
%             Core_member(i,q,1)=input(1,j);
%             q=q+1;
%         end
%     end
%     C(i,1)=q-1;
% end
% 
% for j=1:size(input,2)
%     for i=1:Cluster
%         Distance2(1,i,j)= sqrt((input(1,j)-Cores(1,k))^2);
%         if input(1,j)==Cores(s,i);
%             D(s,i,j)=0;
%             Utility(l+1,i,j)=1;
%         else
%             D(l+1,i,j)= sqrt((input(1,j)-Cores(l+1,i))^2);
%         end
%     end
% end
