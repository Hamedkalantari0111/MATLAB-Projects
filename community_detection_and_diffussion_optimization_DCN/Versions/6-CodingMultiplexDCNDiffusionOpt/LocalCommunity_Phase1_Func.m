function [local_Community,local_Community_size]=LocalCommunity_Phase1_Func(edge_Sorted,Adja_Mat,Nodes,Resolution)
WW=1;
XX=1;
ZZ=1;
LC_Matrix=0;
MM=0;
for w=1:size(edge_Sorted,1)
    NC=0;
    LC=zeros(1,2);
    LC(1,:)=edge_Sorted(w,:); % Local Comunity
    [N]=NeighbourNodes(w,Adja_Mat,LC); %% N=Neighbour nodes of the above LC(i,j)
    %==================Find NC (Neighbour Community)=====================
    [NC]=CommonNeighbourNodes(w,N,NC); %% N=Neighbour nodes of the above LC(i,j)
    if size(NC,2)>=Resolution
        %==============find M measure===============
        [Delta,landa]=DeltaLandaFunc(LC,N,Nodes,w); %calculate Delta and Landa
        %============== the M measure============
        [MM]=MFunction(Nodes,Adja_Mat,Delta,landa); % Calculating M measure
        M(w,1)=MM; % M measure
        LC_temp=zeros(1,1);
        for m=1:size(NC,2)
            LC_temp=LC;
            LC_temp(1,size(LC_temp,2)+1)=NC(1,m);
            [N_temp]=NeighbourNodes(w,Adja_Mat,LC_temp); %% N_temp=Neighbour nodes of the above LC(i,j)
            %===============find M measure================
            [Delta_temp,landa_temp]=DeltaLandaFunc(LC_temp,N_temp,Nodes,w); %calculate Delta and Landa
            %==========the M measure===========
           [MM_temp]=MFunction(Nodes,Adja_Mat,Delta_temp,landa_temp); % Calculating M measure
            M_temp(m,w)=MM_temp; % M measure
            %==========assigning node to LC===========
            if  M_temp(m,w)>M(w,1)
                for y=1:size(LC_temp,2)
                    local_C(WW,y)=LC_temp(1,y);
                end
                LC=LC_temp;
                M(w,1)=M_temp(m,w);
                local_M(WW,1)=M_temp(m,w);
                WW=WW+1;
            end
        end
        local_Community(XX,1:size(local_C,2))=local_C(WW-1,:); % local_Community: Local Communities from phase 1 extracted
        XX=XX+1;
        ZZ=XX;
    end
end
for i=1:size(local_Community,1)
    local_Community_size(i,1)=size(local_Community,2)-sum(local_Community(i,:)==0); % the size of the Final Community
end