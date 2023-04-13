clear all
clc
input=[1,2,4,7,12,14,25,27,30,33,35,40];
Cluster=3;
M=2;
K=1;
N=1000;
E=0.01; % Epsilon Value
Cores = zeros(N+1,Cluster);
Distance=zeros(N,Cluster,size(input,2)); %Distance of each nodes from each core
D=zeros(N+1,Cluster,size(input,2)); %D measure
Utility=zeros(N+1,Cluster,size(input,2)); %membership Function of each nodes to each core
Diff=zeros(N+1,Cluster); %Diff od utility each nodes to each core
noo=zeros(N,Cluster); %meghdare "noo" baraye har core
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