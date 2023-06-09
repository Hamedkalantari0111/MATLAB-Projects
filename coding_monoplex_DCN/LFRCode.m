function [LFRIs,Adjmat,charmat,M_parameter,d_average,C_Star]=LFRCode(kmin,kmax,gama,miu,N,smin,smax,betta)
% generating LFR benchmark datasets and their real group identifiers
% according to "Benchmark graphs for testing community detection algorithms"
% AmirMohsen KarimiMajd - 28/5/92
ns=smax-smin+1;
nsdist=[(smin:smax)' zeros(ns,1)];
nsdist(:,2)=nsdist(:,1).^(-betta);
nsdist(:,2)=nsdist(:,2)./sum(nsdist(:,2));
for i=2:ns
    nsdist(i,2)=nsdist(i-1,2)+nsdist(i,2);
end
nk=kmax-kmin+1;
nkdist=[(kmin:kmax)' zeros(nk,1)];
nkdist(:,2)=nkdist(:,1).^(-gama);
nkdist(:,2)=nkdist(:,2)./sum(nkdist(:,2));
for i=2:nk
    nkdist(i,2)=nkdist(i-1,2)+nkdist(i,2);
end
csvector=zeros(2,1);
i=1;
while sum(csvector)<N
    r=rand;
    candid=nsdist(nsdist(:,2)>=r,:);
    csvector(i)=candid(1,1);
    i=i+1;
end
ssc=sum(csvector);
lc=length(csvector);
for i=1:lc-1
    d=ceil((ssc-N)*(csvector(i)-smin)/ssc);
    csvector(i)=csvector(i)-d;    
end
ssc=sum(csvector);
csvector(lc)=csvector(lc)+N-ssc;
adjmat=zeros(N);
charmat=[(1:N)' zeros(N,3)];
j=1;
ll=csvector(j);
for i=1:N
    if i<=csvector(j)
        charmat(i,2)=j;
    else
        j=j+1;
        charmat(i,2)=j;
        csvector(j)=csvector(j-1)+csvector(j);
    end
    r=rand;
    candid=nkdist(nkdist(:,2)>=r,:);
    charmat(i,3)=ceil((1-miu)*candid(1,1));
    charmat(i,4)=candid(1,1)-charmat(i,3);
end
for i=1:lc
    s=charmat(charmat(:,2)==i,:);
    ss=s(:,3:4);
    ss=flipud(sortrows(ss,1));
    charmat(s(1,1):s(end,1),:)=[s(:,1:2) ss];
end
for i=1:N
    sa=charmat(charmat(:,2)==charmat(i,2),1);
    s=sa;
    ls=length(s);
    s(s==i)=0;
    t=charmat(charmat(:,2)~=charmat(i,2),1);
    lt=length(t);
    g=1;
    while charmat(i,3)>0
        j=ls-g+1;
        if s(j)>0
            if charmat(s(j),3)>0
                charmat(s(j),3)=charmat(s(j),3)-1;
                charmat(i,3)=charmat(i,3)-1;
                adjmat(i,s(j))=1;
                adjmat(s(j),i)=1;
                g=g+1;
            else
                g=g+1;
            end
        else
            g=g+1;
        end
        if g==ls
            break;
        end
    end
    
    while charmat(i,4)>0
        j=lt-g+1;
        if t(j)>0
            if charmat(t(j),4)>0
                charmat(t(j),4)=charmat(t(j),4)-1;
                charmat(i,4)=charmat(i,4)-1;
                adjmat(i,t(j))=1;
                adjmat(t(j),i)=1;
                g=g+1;
            else
                g=g+1;
            end
        else
            g=g+1;
        end
        if g==lt
            break;
        end
    end
    a=kmin-sum(adjmat(i,:));
    if a>0
        b=(1-adjmat(i,sa)).*sa';
        b=b(b>0);
        for j=1:a
            adjmat(i,b(j))=1;
            adjmat(b(j),i)=1;
        end
    end
end
LFRIs=[charmat(:,1:2) adjmat];
Adjmat=adjmat;
charmat;
%csvwrite('LFRIs.csv', [charmat(:,1:2) adjmat])
M_parameter=sum(sum(adjmat))/2;
%disp(sum(sum(adjmat))/2) % is "M"
d_average=sum(sum(adjmat))/350;
%disp(sum(sum(adjmat))/350) % is "d average"
C_Star=lc;
%disp(lc) % is "C*"