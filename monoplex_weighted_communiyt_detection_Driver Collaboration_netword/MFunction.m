function[M]= MFunction(Nodes,Adja_Mat,Delta,landa)
MM_in(1:Nodes,1:Nodes)=Adja_Mat(1:Nodes,1:Nodes).*Delta(1:Nodes,1:Nodes);
MM_out(1:Nodes,1:Nodes)=Adja_Mat(1:Nodes,1:Nodes).*landa(1:Nodes,1:Nodes); 
M1_in(1:Nodes)=sum(MM_in(1:Nodes,:));
M1_out(1:Nodes)=sum(MM_out(1:Nodes,:));
M_in=0.5*sum(M1_in(:));
M_out=0.5*sum(M1_out(:));
M=M_in/M_out;