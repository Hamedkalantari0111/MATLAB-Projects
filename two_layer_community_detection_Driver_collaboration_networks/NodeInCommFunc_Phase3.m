function[N_in,N_in_size]= NodeInCommFunc_Phase3(Final_Community,Final_Community_size,Adja_Mat,Nodes)
k=0;
N_in=0;
for i=1:size(Final_Community,1) % i= Comunity counter
    for j=1:Final_Community_size(i,1) % j= node counter in each communtiy
        for l=1:Nodes
            if Adja_Mat(Final_Community(i,j),l)==1
                k=k+1;
                N_in(i,k,j)=l; % N=the neighbour of node i (node mojod dar har com)i=com, j=node in com, k= numbers of the neighbours
            end
        end
        k=0;
    end
    k=0;
end
  %==========size of the "neighbour nodes in community matrix"=======
for i=1:size(N_in,1)
    for j=1:size(N_in,3)
        N_in_size(i,j)=size(N_in,2)-sum(N_in(i,:,j)==0); 
    end
end