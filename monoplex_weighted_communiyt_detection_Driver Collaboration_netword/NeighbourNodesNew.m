function [N]=NeighbourNodesNew(w,Adja_Mat,Adja_Mat_L2,LC,Nodes)
k=0;
for i=1:size(LC,2)
    for j=1:Nodes
        if Adja_Mat(LC(1,i),j)==1
            k=k+1;
            N_L1(i,k,w)=j; % N=Neighbour nodes of the above LC(i,j) in layer #1
        end
    end
    k=0;
end
k=0;
for i=1:size(LC,2)
    for j=1:Nodes
        if Adja_Mat_L2(LC(1,i),j)==1
            k=k+1;
            N_L2(i,k,w)=j; % N=Neighbour nodes of the above LC(i,j) in layer #2
        end
    end
    k=0;
end

for i=1:size(LC,2)
    S=union(N_L1(i,1:size(N_L1,2)-sum(N_L1(i,:,w)==0),w),N_L2(i,1:size(N_L2,2)-sum(N_L2(i,:,w)==0),w));
    N(i,1:size(S,2),w)=S;
end

    