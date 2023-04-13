clear all
clc
tic;
 Input=dlmread('zakhari.txt');
 Inputphase1=dlmread('zakhariphase1.txt');
% Inputphase1=dlmread('abrarphase1.txt');
 maxi=max(max(Input));
 [InRow, InClm]=size(Input);
 Adja_Mat = zeros(maxi,maxi); 
 for i=1:InRow
     Adja_Mat(Input(i,1),Input(i,2))= 1;
     Adja_Mat(Input(i,2),Input(i,1))= 1; 
 end
% maxi=163;
% Adja_Mat =dlmread('abrar.txt');
%=======================phase1=====================%
sortnodes_as = sortrows(Inputphase1,2);
for i=1:maxi
    for j=1:2
    sortnodes_ds(maxi+1-i,j)=sortnodes_as(i,j);
    end
end
t=1;
j=1;
k=1;
SE=zeros((maxi/2),2);
for i=1:maxi
        if ismember(sortnodes_ds(i,1),SE)~=1;
           if i==maxi
              SN(k,1)= sortnodes_ds(i,1); 
           end
            for t=1:(maxi-i)
                 if Adja_Mat(sortnodes_ds(i,1),sortnodes_ds(i+t,1))==1;
                  SE(j,1)= sortnodes_ds(i,1);
                  SE(j,2)= sortnodes_ds(i+t,1);
                  j=j+1;
                  break
                 end
                 if i+t>=maxi;
                     SN(k,1)= sortnodes_ds(i,1);
                     k=k+1;
                 end
            end
        end
end
%=======================phase2=====================%
maxj=size(SE,1);
LC=SE;
o=1;
for i=1:maxi
         for j=1:maxj
            if LC(j,1)~=0;
                if Adja_Mat(LC(j,1),i)==1
                    if Adja_Mat(LC(j,2),i)==1
                        LC(j,2+o)=i;
                        o=o+1;
                    end
                end
            end
        end
end
k=1;
for i=1:size(LC,1);
    for j=1:size(LC,2);
        if LC(i,j)~=0;
            LCphase2(i,k)=LC(i,j);
            k=k+1;
        end
    end
    k=1;
end
%=======================phase3=====================%
alpha=0.8;
beta=0.4;
for k=1:size(LCphase2,1);
    e=1;
    for i=1:size(LCphase2,2);
        for j=1:(size(LCphase2,2)-i);
            if LCphase2(k,i+j)~=0
                if Adja_Mat(LCphase2(k,i),LCphase2(k,i+j))==1;
                    edges(k,e)=LCphase2(k,i);
                    edges(k,e+1)=LCphase2(k,i+j);
                    e=e+2;
                end
            end
        end
    end
end
sub_edges=zeros(size(LCphase2,1),size(LCphase2,1));
for k=1:size(edges,1);
    for i=1:(size(edges,2)-1);
        for p=(k+1):size(edges,1);
            oo=intersect(LCphase2(k,:),LCphase2(p,:));
            mm=size(oo,2);
            if ismember(0,oo)==1;
                mm=mm-1;
            end
            if mm>=2;
                for j=1:((size(edges,2)/2)-1);
                    if edges(p,((2*j-1)))>0 && edges(k,i)>0;
                        if edges(k,i)==edges(p,((2*j)-1))&& edges(k,i+1)==edges(p,((2*j)-1)+1);
                            sub_edges(k,p)= sub_edges(k,p)+1;
                        end
                    end
                end
            end         
         end
    end
end

L=zeros(size(LCphase2,1),1);
WOS=zeros(size(LCphase2,1),size(LCphase2,1));
uu=zeros(size(LCphase2,1),1);
for i=1:size(LCphase2,1);
    for u=1:size(LCphase2,2);
       if LCphase2(i,u)>0;
           uu(i)= uu(i)+1;                  
       end
    end
end
kk=zeros(size(LCphase2,1),1);
for i=1:size(LCphase2,1);
    for u=1:size(edges,2);
       if edges(i,u)>0;
           kk(i)= kk(i)+1;                  
       end
    end
end
yy=zeros(size(LCphase2,1),size(LCphase2,1));
for k=1:size(LCphase2,1);
    for p=(k+1):size(LCphase2,1);
        for i=1:size(LCphase2,2);
            for j=1:size(LCphase2,2);
                if LCphase2(k,i)>0 && LCphase2(p,j)>0;                    
                    if Adja_Mat(LCphase2(k,i),LCphase2(p,j))>0;
                        yy(k,p)=yy(k,p)+1;
                        if ismember(LCphase2(k,i),LCphase2(p,:))==1 && ismember(LCphase2(p,j),LCphase2(k,1:(i-1)))==1;
                          yy(k,p)=yy(k,p)-1;  
                        end
                    end
                end
            end
        end
    end
end
           
k=1;
m=1;
for i=1:size(LCphase2,1);
    for j=1:(size(LCphase2,1)-i);
       if isempty(intersect(LCphase2(i,:),LCphase2(i+j,:)))==0;
               eshterak_nodes(i,i+j)=size(intersect(LCphase2(i,:),LCphase2(i+j,:)),2);
               if ismember(0,intersect(LCphase2(i,:),LCphase2(i+j,:)))==1;
                  eshterak_nodes(i,i+j)=eshterak_nodes(i,i+j)-1; 
               end
               part1(i,i+j)=(eshterak_nodes(i,i+j)/min(uu(i),uu(i+j)));
               part2(i,i+j)=(yy(i,i+j)/(min((kk(i)/2),(kk(i+j)/2))+(yy(i,i+j)-sub_edges(i,i+j))));
               %part2(i,i+j)=(sub_edges(i,i+j)/min((kk(i)/2),(kk(i+j)/2)));
               WOS(i,i+j)=(alpha*part1(i,i+j))+((1-alpha)*part2(i,i+j));
             if WOS(i,i+j)>=beta;
                 for u=1:size(LCphase2,2);
                     if LCphase2(i,u)>0;
                        Cphase3(k,m)=LCphase2(i,u);
                        m=m+1;
                     end
                 end
                  for u=1:size(LCphase2,2);
                     if LCphase2(i+j,u)>0;
                        Cphase3(k,m)=LCphase2(i+j,u);
                        m=m+1;
                     end
                 end
                 k=k+1;
                 L(i,i+j)=1;
                 L(i+j,i)=1;
                 m=1;
             end             
        end
    end
end
SUM=0;
for i=1:size(LCphase2,1);
    for j=1:(size(LCphase2,1)-i);
        SUM=SUM+WOS(i,i+j);
    end
end
mean_WOS=SUM/(((size(LCphase2,1)*size(LCphase2,1))-(size(LCphase2,1)))/2);
for i=1:size(LCphase2,1);
    if ismember(1,L(i,:))==0;
        for u=1:size(LCphase2,2);
            if LCphase2(i,u)>0;
                Cphase3(k,m)=LCphase2(i,u);
                m=m+1;
            end
        end
        k=k+1;
        m=1;
    end
end
Cphase3_ver1=zeros(size(Cphase3,1),size(Cphase3,2));
ii=1;
for i=1:size(Cphase3,1);
    ii=1;
    for u=1:size(Cphase3,2);
        if ismember(Cphase3(i,u),Cphase3_ver1(i,:))==0;
            Cphase3_ver1(i,ii)=Cphase3(i,u);
            ii=ii+1;
        end
    end
end
jj=1;
cnodes=zeros(size(Cphase3_ver1,1),1);
for i=1:size(Cphase3_ver1,1);
    for u=1:size(Cphase3_ver1,2);
            if Cphase3_ver1(i,u)>0;
               cnodes(i)=cnodes(i)+1;
            end
    end
end
tarkib=zeros(size(Cphase3_ver1,1),size(Cphase3_ver1,1));
eshterak_Cphase3=zeros(size(Cphase3_ver1,1),size(Cphase3_ver1,1));
for i=1:size(Cphase3_ver1,1);
    for j=1:(size(Cphase3_ver1,1)-i);
        eshterak_Cphase3(i,i+j)=size(intersect(Cphase3_ver1(i,:),Cphase3_ver1(i+j,:)),2);
        if ismember(0,intersect(Cphase3_ver1(i,:),Cphase3_ver1(i+j,:)))==1;
           eshterak_Cphase3(i,i+j)=eshterak_Cphase3(i,i+j)-1;
        end
    end     
end
Cphase3_ver2(1,1:size(Cphase3_ver1,2))=zeros(1,1:size(Cphase3_ver1,2));
for i=1:(size(Cphase3_ver1,1)-1);
     for j=1;
         if ismember(1,tarkib(:,i))~=1;
            if eshterak_Cphase3(i,i+j)>=(0.6*min(cnodes(i),cnodes(i+j)));
               a=size(union(Cphase3_ver1(i,:),Cphase3_ver1(i+j,:)),2);
               Cphase3_ver2(jj,1:a)=union(Cphase3_ver1(i,:),Cphase3_ver1(i+j,:));            
               tarkib(i,i+j)=1;
                             
            end
         end
     end
    for j=2:(size(Cphase3_ver1,1)-i);
        if ismember(1,tarkib(:,i))~=1;
            if eshterak_Cphase3(i,i+j)>=(0.6*min(cnodes(i),cnodes(i+j)));
               if size(Cphase3_ver2,1)<jj;
                    Cphase3_ver2(jj,1:size(Cphase3_ver1(i+j,:),2))=zeros(1,size(Cphase3_ver1(i+j,:),2));
                end
                a=size(union(Cphase3_ver2(jj,:),Cphase3_ver1(i+j,:)),2);
                Cphase3_ver2(jj,1:a)=union(Cphase3_ver2(jj,:),Cphase3_ver1(i+j,:));            
                tarkib(i,i+j)=1;                             
            end
        end
    end 
    if ismember(1,tarkib(i,:))==1;
        jj=jj+1; 
    else if ismember(1,tarkib(i,:))~=1 && ismember(1,tarkib(:,i))~=1;                       
        b=size(Cphase3_ver1(i,:),2);
        Cphase3_ver2(jj,1:b)=Cphase3_ver1(i,:);
        jj=jj+1;
        end
        if i==(size(Cphase3_ver1,1)-1);
            if ismember(1,tarkib(i+1,:))~=1 && ismember(1,tarkib(:,i+1))~=1;                       
                b=size(Cphase3_ver1(i+1,:),2);
                Cphase3_ver2(jj,1:b)=Cphase3_ver1(i+1,:);
            end
        end
    end 
end
for i=1:size(Cphase3_ver2,1);    
    if Cphase3_ver2(i,1)==0;
        for u=2:size(Cphase3_ver2,2);
        Cphase3_ver2(i,u-1)=Cphase3_ver2(i,u);
        end
    end
end
%=======================phase4=====================%
M_in=zeros(size(Cphase3_ver2,1),1);
M_out=zeros(size(Cphase3_ver2,1),1);
eq=0.4;
for i=1:size(Cphase3_ver2,1);
    for u=1:size(Cphase3_ver2,2);
        if Cphase3_ver2(i,u)~=0;
            for j=1:size(Inputphase1,1);
                if Adja_Mat(Cphase3_ver2(i,u),Inputphase1(j,1))==1 && ismember(Inputphase1(j,1),Cphase3_ver2(i,:))==0;
                    M_out(i)= M_out(i)+1;
                end
            end
        end

    end
end
for i=1:size(Cphase3_ver2,1);
    for u=1:size(Cphase3_ver2,2)
        if Cphase3_ver2(i,u)~=0;
            for j=1:size(Inputphase1,1);
                if Adja_Mat(Cphase3_ver2(i,u),Inputphase1(j,1))==1 && ismember(Inputphase1(j,1),Cphase3_ver2(i,:))==1;
                    M_in(i)= M_in(i)+1;           
                end
            end
        end

    end
end
for i=1:size(Cphase3_ver2,1);
modularity_ver1(i)=M_in(i)/M_out(i);
end

for d=1:size(SN,1);
    M_in=zeros(size(Cphase3_ver2,1),1);
    M_out=zeros(size(Cphase3_ver2,1),1);
    tedad_nodes=zeros(size(Cphase3_ver2,1),1);
        for i=1:size(Cphase3_ver2,1);
            for u=1:size(Cphase3_ver2,2);
                    if Cphase3_ver2(i,u)>0;
                       tedad_nodes(i)=tedad_nodes(i)+1;
                    end
            end
        end
    if ismember(SN(d),Cphase3_ver2)==0;
        for i=1:size(Cphase3_ver2,1);
            Cphase3_ver2(i,(tedad_nodes(i)+1))=SN(d);
            for u=1:(tedad_nodes(i)+1);
                if Cphase3_ver2(i,u)~=0;
                    for j=1:size(Inputphase1,1);
                        if Adja_Mat(Cphase3_ver2(i,u),Inputphase1(j,1))==1 && ismember(Inputphase1(j,1),Cphase3_ver2(i,:))==0;
                            M_out(i)= M_out(i)+1;
                        end
                    end
                end

            end
        end  
       for i=1:size(Cphase3_ver2,1);           
            for u=1:(tedad_nodes(i)+1);
                if Cphase3_ver2(i,u)~=0;
                    for j=1:size(Inputphase1,1);
                        if Adja_Mat(Cphase3_ver2(i,u),Inputphase1(j,1))==1 && ismember(Inputphase1(j,1),Cphase3_ver2(i,:))==1;
                            M_in(i)= M_in(i)+1;           
                         end
                    end
                end

            end
       end    
        for i=1:size(Cphase3_ver2,1);
            modularity_ver2(i)=M_in(i)/M_out(i);
            diff(i)=modularity_ver2(i)-modularity_ver1(i);            
        end
        [z,x]=max(diff);
        for i=1:size(Cphase3_ver2,1);            
             if i~=x;
                Cphase3_ver2(i,(tedad_nodes(i)+1))=0;
             end            
        end 
    end
end
for i=1:size(Cphase3_ver2,1);
    for u=2:(size(Cphase3_ver2,2)-1);
        if Cphase3_ver2(i,u)==Cphase3_ver2(i,u-1);        
        Cphase3_ver2(i,u)=Cphase3_ver2(i,u+1);
            if u+1==size(Cphase3_ver2,2);
               Cphase3_ver2(i,u+1)=0;
            end
        end
    end
end
%=======================EQ=====================%
m_edges=(sum(sum(Adja_Mat))/2);
degree=zeros(size(Adja_Mat,2),1);
O_count=zeros(size(Adja_Mat,2),1);
for i=1:size(Adja_Mat,2);
    degree(i)=sum(Adja_Mat(i,:));
end
for i=1:size(Adja_Mat,2);
    degree(i)=sum(Adja_Mat(i,:));
end
for i=1:size(Adja_Mat,2);
    for j=1:size(Cphase3_ver2,1);
        if ismember(Inputphase1(i,1),Cphase3_ver2(j,:))==1;
           O_count(i)=O_count(i)+1;        
        end
    end
end
EQ_ver1=0;
for c=1:size(Cphase3_ver2,1);
    for u=1:(tedad_nodes(c)-1);        
        for v=(u+1):tedad_nodes(c);
            if Cphase3_ver2(c,u)~=0 && Cphase3_ver2(c,v)~=0;
            EQ_ver1=EQ_ver1+((1/(O_count(u,1)*O_count(v,1)))*(Adja_Mat(u,v)-((degree(u,1)*degree(v,1))/(2*m_edges))));            
            end
        end
    end
end
EQ=(0.5/m_edges)*(EQ_ver1)+eq;
toc
