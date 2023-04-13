function[NeighNodesLocalCom]= NeighNodeLocalCommFunc(local_Community,Nodes,Adja_Mat)
k=0;
 for i=1:size(local_Community,1)
     for l=1:size(local_Community,2)
         for j=1:Nodes 
             if local_Community(i,l)~=0
                 if Adja_Mat(local_Community(i,l),j)==1
                     k=k+1;
                     AA(i,k)=j;
                 end
             end
         end
     end
     k=0;
     BB=unique(AA(i,:));
     CC=size(BB,2);
     NeighNodesLocalCom(i,1:CC)=BB; % NeighNodesLocalCom= all of the neighbour node of the community nodes
 end