clear all
clc
input=[1,2,4,7,12,14,25,27,30,33,35,40];
Cluster=3;
M=3;
K=1;
N=1000;
E=0.001; % Epsilon diff Value
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
        if D(l+1,i,j)<noo(l+1,i)
            V(i,q)=input(1,j);
            q=q+1;
        end
       end
       C(i,1)=q-1;
 end
 VV=0;
 Nooo=0.79; %Nooo is the fraction of beta
 Itration=10000;

 for i=1:Cluster
     for j=1:C(i,1)
         for jprim=1:C(i,1)
             if j>jprim || j<jprim
             VV(i,j,jprim)=sqrt((V(i,j)-V(i,jprim)).^2);
             landa(i,j,jprim)=Nooo*noo(l+1,i)/VV(i,j,jprim);
             end
         end
     end

%      VV=0;
 end
 
 % mohasebe U^ %lower utility_hat_normalized
 
 for i=1:Cluster
     for j=1:C(i,1)
         for k=1:size(input,2)
             Dist(i,j,k)= sqrt((input(1,k)-V(i,j)).^2); % Dist=Distance new
             utility_hat(i,j,k)=1/(1+(Dist(i,j,k)/noo(l+1,i)).^(1/(M-1)));
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
             Distt_up= Distt_up+ sqrt((V(i,j)-V(i,k)).^2); 
             for p=1:C(i,1)
                 Distt_down=Distt_down+ sqrt((V(i,p)-V(i,k)).^2); 
             end
             Disttdown(i,j,k)=Distt_down;
             Distt_down=0;
         end
         Disttup(i,j)= Distt_up;
         for k=1:C(i,1)
             sub=sub+Disttup(i,j)/Disttdown(i,j,k);
         end
         Utility_hathat(i,j)=sub;
         sub=0;
         Distt_up=0;
     end
     norm_down(i)=sum(Utility_hathat(i,:));
     for j=1:C(i,1)
         utility_hathat_normalized(i,j)= Utility_hathat(i,j)/norm_down(i);
     end
  end
  
%Calculate Core Membera of the cluster Wij
AA=0;
BB=0;
CC=0;
DD=0;
for i=1:Cluster
    for j=1:C(i,1)
        for k=1:size(input,2)
            AA= AA + ((utility_hat(i,j,k)).^M)*input(1,k);
            BB= BB + ((utility_hat(i,j,k)).^M);
        end
        part1_up(i,j)=AA;
        AA=0;
        part1_down(i,j)=BB;
        BB=0;
        for l=1:C(i,1)
            CC= CC +(((utility_hathat_normalized(i,j).^M) + (landa(i,j,l).^M))*V(i,l));
            DD= DD + (landa(i,j,l).^M);
        end
        part2_up(i,j)=CC;
        CC=0;
        part3_down(i,j)=DD;
        DD=0;
        Core_member(i,j)=(part1_up(i,j)+part2_up(i,j))/(part1_down(i,j)+(C(i,1)*(utility_hathat_normalized(i,j).^M))+part3_down(i,j))
    end
end
%Algorithm l=Itration, i=Core, j=Nodes, 



















