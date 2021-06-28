%ALL=xlsread('C:\Users\carcei01\desktop\WORK PROJECTS\MATERNAL CARE\Co-housing data\V85\V85_12102017_all units timestamps');

ALL; %the array of all channel timestamps
bt; %behavioral timestamps
MI=zeros(size(ALL,2),1); % empty array for storing the mod indexes for all channels
%P=zeros(size(ALL,2),1); %empty array for storing p-values for all cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% place max limits on behavioral episodes
b=bt;
% for o=1:length(b)
%     if (b(o,2)-b(o,1))>=40;
%         b(o,2)=b(o,1)+40;
%     else
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adj=0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for calculating MI of the REAL behavioral episode
c=zeros(size(b,1),1);
c1=zeros(size(b,1),1);
c(:)=(b(:,1)-adj);
c1(:)=c(:)-(b(:,2)-b(:,1));
%c1(:)=c(:)-0.2;

for j=1:size(ALL,2);

a=ALL(:,j); %unit time stamps

a(isnan(a))=[];
d=zeros(size(b)); %matrix for average frequencies 

    
% calculating firing rate during the behavior episode
for i=1:size(b,1);
    x=find(a>=(b(i,1)) & a<=b(i,2));
    
   %x=find(a>(b(i,1)-1) & a<(b(i,1)+10));
  if isempty(x)
      d(i,1)=0/(b(i,2)-b(i,1));
  else
      d(i,1)=length(x)/(b(i,2)-b(i,1));
      %d(i,1)=length(x)/10;
  end 
end

% calculating the firing rate during the baseline
for i=1:size(b,1);
  y=find(a>=c1(i,1) & a<=c(i,1));
    if isempty(y)
      d(i,2)=0/50;
    else
    d(i,2)=length(y)/(c(i,1)-c1(i,1));
   
    end
   
end


% calculate mean modulation index for each cell
mi=((d(:,1)-d(:,2)))*100./(d(:,2)+d(:,1));
mi(find(isnan(mi)))=[];
mmi=mean(mi,1);
mmi(find(isnan(mmi)))=[];

% add data for each cell to the existing matrix
MI(j,:)=mmi;
realMI=MI;

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for randomization, to calculate significance
B=zeros(size(b));
bh=b;
r = randi([-500 500],1,1000);
RMI= zeros(size(ALL,2),size(r,2));

for g=1:size(r,2)
    bh=b+r(g);
    
   B(:,1)=bh(:,1);
   B(:,2)=bh(:,2);
  


% determine the interval for the baseline
C=zeros(size(B,1),1);
C1=zeros(size(B,1),1);
C(:)=B(:,1)-adj;
C1(:)=C(:)-(B(:,2)-B(:,1));
%C1(:)=C(:)-5;

for j=1:size(ALL,2);

a=ALL(:,j); %unit time stamps

a(isnan(a))=[];
d=zeros(size(B)); %matrix for average frequencies 

    
% calculating firing rate during the behavior episode
for i=1:size(b,1);
    x=find(a>=(B(i,1)) & a<=B(i,2));
    
   %x=find(a>(b(i,1)-1) & a<(b(i,1)+10));
  if isempty(x)
      d(i,1)=0/(B(i,2)-B(i,1));
  else
      d(i,1)=length(x)/(B(i,2)-B(i,1));
      %d(i,1)=length(x)/10;
  end 
end

% calculating the firing rate during the baseline
for i=1:size(b,1);
  y=find(a>=C1(i,1) & a<=C(i,1));
    if isempty(y)
      d(i,2)=0/50;
    else
    d(i,2)=length(y)/(C(i,1)-C1(i,1));
   
    end
   
end

%calculate significance for each cell
%[p,h]=ranksum(d(:,1),d(:,2));

% calculate mean modulation index for each cell
mi=((d(:,1)-d(:,2)))*100./(d(:,2)+d(:,1));
mi(find(isnan(mi)))=0;
mmi=mean(mi,1);
mmi(find(isnan(mmi)))=0;


% add data for each cell to the existing matrix
MI(j,:)=mmi;
%P(j,:)=p;

end
RMI(:,g)=MI;
%MI(find(isnan(MI)))=[];

end

RMI;  
S=zeros(size(RMI,1),1);
for f=1:size(RMI,1)
  S(f,1)=prctile(RMI(f,:),95);
  S(f,2)=prctile(RMI(f,:),5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate p-vals from permutation tests
Above=zeros(size(RMI,1),1);
for e=1:size(ALL,2)
    if realMI(e,1)>=0
    Above(e,1)=size(find(RMI(e,:)>realMI(e,1)),2);
    else
    Above(e,1)=size(find(RMI(e,:)<realMI(e,1)),2);
    end
    
end


pval(:,1)=Above(:,1)./1000;

sortpvals=sort(pval,'ascend');
FDR=zeros(length(sortpvals),3);
FDR(:,1)=sortpvals;
FDR(:,2)=(1:1:size(FDR,1));
FDR(:,3)=(FDR(:,2)/size(FDR,1))*0.05

% mat2clip([MI,P])
mat2clip(FDR)
clear all
