function [Stat,BestValue,XTarget,pdm]=DLRao(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X)
% fhd,nPop,nVar,VarMin,VarMax,MaxIt,X
%%Input parameters
%%fhd----------------objective function
%%nPop---------------population size 
%%nVar---------------the number of variables
%%VarMin-------------the lower boundaries of variables
%%VarMin-------------the upper boundaries of variables
%%MaxIt--------------the maximum number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Output parameters
%%BestCost-----------convergence curve
%%BestValue----------the optimal fitness value
%%XTarget------------the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nPop
%     X(i,:)=VarMin+rand(1,nVar).*(VarMax-VarMin); %Initial population
    f(i) = fhd(X(i,:));
end
pdm=[];
gen=1;% Initial the number of iterations
[Best_Cost,ind]=min(f);% 
Div=dispersion(X);
Stat(1,:)=[Best_Cost median(f) iqr(f) Div];
XTarget=X(ind,:);M_X=rand(1,nVar);
for gen=1:MaxIt % Main loop
    [t,tindex]=sort(f);
    Best=X(tindex(1),:); 
    worst=X(tindex(end),:);
    xnew=zeros(nPop,nVar);
    for ii=1:nPop;
        uF(tindex(ii))=ii/nPop;
    end
    M=mean(xnew);
    for i=1:nPop
        xnew(i,:)=CL(i,uF,Best,worst,M,X,nPop,nVar);
    end
    for i=1:nPop
        xnew(i,:) = BC(xnew(i,:),VarMin,VarMax,nVar);%max(min(xnew(i,:),VarMax),VarMin);% boundary limit
        fnew(i,:) = fhd(xnew(i,:));
    end
    for i=1:nPop
        if(fnew(i)<f(i))
            X(i,:) = xnew(i,:);
            f(i) = fnew(i);
            pdm=[pdm;X(i,:) f(i)];
        end
    end
    [Best_Cost,ind]=min(f);
    Div=dispersion(X);
    Stat(gen,:)=[Best_Cost median(f) iqr(f) Div];
    XTarget=X(ind,:);
end
BestValue=min(f);
end

function newsol=CL(i,uF,Best,worst,M,pop,nPop,nVar)
j=randperm(nPop,1);
k=randperm(nPop,1);
while k==j||j==i||k==i
    j=randperm(nPop,1);
    k=randperm(nPop,1);
end
if uF(i)>=1/2
    [F1,F2]=Ffactors(i,j,uF);
    newsol = pop(i,:)+ F1*(Best-pop(i,:))+F2*(pop(k,:)-pop(j,:));
else
    r1=rand(1,nVar);
    r2=rand(1,nVar);
    delta1=r1.*(Best-worst);
    if rand<0.5
        newsol = pop(i,:)+delta1+r2.*(pop(k,:)-pop(j,:));
    else
        newsol = pop(i,:)+delta1-r2.*(M-pop(k,:));
    end
end
end 

function [F1,F2]=Ffactors(i,j,uF)

F1=uF(i)+0.1*trnd(1,1,1);% ?从柯西分布
F2=uF(j)+0.1*trnd(1,1,1);
while F2<=0
      F2=uF(j)+0.1*trnd(1,1,1);
end
while F1<=0
    F1=uF(i)+0.1*trnd(1,1,1);
end
if F1>=1 || F2>=1
    F1=1;
    F2=1;
end
end