function[M_Measure]= Mmeasure_phase3Func(Final_Community,Nodes,Adja_Mat,delta,landa)
for k=1:size(Final_Community,1)
    for i=1:Nodes
        for j=1:Nodes
            if i~=j
                MMM_in(i,j,k)=Adja_Mat(i,j)*delta(i,j,k);
                MMM_out(i,j,k)=Adja_Mat(i,j)*landa(i,j,k);
            end
        end
        MM1_in(i,k)=sum(MMM_in(i,:,k));
        MM1_out(i,k)=sum(MMM_out(i,:,k));
    end
    MM_in(k,1)=0.5*sum(MM1_in(:,k));
    MM_out(k,1)=0.5*sum(MM1_out(:,k));
    M_Measure(k,1)=MM_in(k,1)/MM_out(k,1); % M measur of the fina communities
end