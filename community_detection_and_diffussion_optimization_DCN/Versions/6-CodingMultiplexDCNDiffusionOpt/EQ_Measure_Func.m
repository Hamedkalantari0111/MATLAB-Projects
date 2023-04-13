function [Q]=EQ_Measure_Func(Input,Nodes,FinalCommunity,Community_size,Adja_Mat)
Input_Size=size(Input,1);
s = Input(:,1); 
t = Input(:,2); 
G = graph(s,t);
D = degree(G); % Degree of the nodes!
L=0;
for i=1:Nodes 
    for j=1: size(FinalCommunity,1)
        if ismember(i,FinalCommunity(j,:))==1
            L=L+1;
            O(i,1)=L; % the number of communities which node i belongs to.
        end
    end
    L=0;
end
QQ=0;
%(1/EQ_factor)*(EQ_part1-(EQ_part2_up/(2*Input_Size)))
for i=1:size(FinalCommunity,1)
    for j=1:Community_size(i,1) % node u
        for k=1:Community_size(i,1) % node v
            if FinalCommunity(i,j)~=FinalCommunity(i,k)
                EQ_factor=(O(FinalCommunity(i,j),1) * O(FinalCommunity(i,k),1));
                EQ_part1=Adja_Mat(FinalCommunity(i,j),FinalCommunity(i,k));
                EQ_part2_up=D(FinalCommunity(i,j),1)*D(FinalCommunity(i,k),1);
                Result=(1/EQ_factor)*(EQ_part1-(EQ_part2_up/(2*Input_Size))); %(1/1)=(1/EQ_factor) %%(1/(2*Input_Size))*(1/EQ_factor)*(EQ_part1-(EQ_part2_up/(2*Input_Size))); %(1/1)=(1/EQ_factor)
                QQ=QQ+Result;
                end
        end
    end
end
Q=0;
Q=(1/(2*Input_Size))*QQ;