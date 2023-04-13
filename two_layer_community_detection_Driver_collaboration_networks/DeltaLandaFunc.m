function[Delta,landa]= DeltaLandaFunc(LC,N,Nodes,w)
Delta=zeros(Nodes,Nodes); %if both two nodes i,j where in the Community =1 , O.W=0;
landa=zeros(Nodes,Nodes); %if only one of the nodes was in the Community=1 , O.W=0;
CC=zeros(size(LC,2),size(N,2));
DD=zeros(size(LC,2),size(N,2));
for i=1:size(LC,2)
    [C,D]=ismember(N(i,:,w),LC(1,:));
    CC(i,:)=C; % index of Neighbour node in "N"
    DD(i,:)=D; %Neighbour Node
end
for i=1:size(CC,1)
    for j=1:size(CC,2)
        if CC(i,j)==1 
            Delta(LC(1,i),LC(1,DD(i,j)))=1; %if both two nodes i,j (of the edge) were in the Community =1 , O.W=0;
            landa(LC(1,i),LC(1,DD(i,j)))=0; %if only one of the nodes (of the edge) was in the Community=1 , O.W=0;
        elseif N(i,j,w)~=0
            Delta(LC(1,i),N(i,j,w))=0;
            Delta(N(i,j,w),LC(1,i))=0;
            landa(LC(1,i),N(i,j,w))=1;
            landa(N(i,j,w),LC(1,i))=1;
        end
    end
end