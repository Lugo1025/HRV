function [m, p, DistHist]=vqsplit(X,L)

% Entradas:
% X: matriz de vectores característicos.
% L: tamaño del codebook (numero de particiones)
% 
% Salidas:
% M: codebook con los cenroides de salida
% P: peso de cada cluster
% DH: distorcion de cada iteracion
%
% Algoritmo:
%1. Hallar la media
% 2.Dividir cada centroide en 2
% 3. Asigna cada data al centroide
% 4. Encontrar los centroides
% 5. Calcular las distancias
% 6. Si la distancia no ha cambiado mucho
%       el numero de centroids es mas pequeño que L2 y regresamos a paso 2
%       sino pasamos al paso 7
%       else Goto 7
%    
% 7. El numero de centroides es mas grande
%    Descartar el centroide con mayor distorsion
% 8. Fin


e=.01; 
eRed=0.75; 
DT=.005; 
DTRed=0.75; 
MinPop=0.10; 
             


d=size(X,1); 
N=size(X,2); 
isFirstRound=1;

if numel(L)==1
    M=mean(X,2); 
    CB=[M*(1+e) M*(1-e)]; 
else
    CB=L; 
    L=size(CB,2);
    e=e*(eRed^fix(log2(L)));
    DT=DT*(DTRed^fix(log2(L)));
end

LC=size(CB,2);

Iter=0;
Split=0;
IsThereABestCB=0;
maxIterInEachSize=20; 
EachSizeIterCounter=0;
while 1
    %Calculo de distancias
    [minIndx, dst]=VQIndex(X,CB); 

    ClusterD=zeros(1,LC);
    Population=zeros(1,LC);
    LowPop=[];
    % Encontrar centroides (Media de cada cluster)
    for i=1:LC
        Ind=find(minIndx==i);
        if length(Ind)<MinPop*N/LC 
            LowPop=[LowPop i];
        else
            CB(:,i)=mean(X(:,Ind),2);
            Population(i)=length(Ind);
            ClusterD(i)=sum(dst(Ind));
        end        
    end
    if ~isempty(LowPop)
        [temp MaxInd]=maxn(Population,length(LowPop));
        CB(:,LowPop)=CB(:,MaxInd)*(1+e); 
        CB(:,MaxInd)=CB(:,MaxInd)*(1-e);
        
        %re-train
        [minIndx, dst]=VQIndex(X,CB);

        ClusterD=zeros(1,LC);
        Population=zeros(1,LC);
        
        for i=1:LC
            Ind=find(minIndx==i);
            if ~isempty(Ind)
                CB(:,i)=mean(X(:,Ind),2);
                Population(i)=length(Ind);
                ClusterD(i)=sum(dst(Ind));
            else 
                CB(:,i)=X(:,fix(rand*N)+1);
                disp('A random vector was assigned as a codeword.')
                isFirstRound=1;
            end                
        end
    end
    Iter=Iter+1;
    if isFirstRound 
        TotalDist=sum(ClusterD(~isnan(ClusterD)));
        DistHist(Iter)=TotalDist;
        PrevTotalDist=TotalDist;        
        isFirstRound=0;
    else
        TotalDist=sum(ClusterD(~isnan(ClusterD)));  
        DistHist(Iter)=TotalDist;
        PercentageImprovement=((PrevTotalDist-TotalDist)/PrevTotalDist);
        if PercentageImprovement>=DT 
            PrevTotalDist=TotalDist; 
            isFirstRound=0;
        else
            EachSizeIterCounter=0;
            if LC>=L 
                if L==LC 
                    disp(TotalDist)
                    break
                else 
                    [temp, Ind]=min(Population); 
                    NCB=zeros(d,LC-1);
                    NCB=CB(:,setxor(1:LC,Ind(1)));
                    CB=NCB;
                    LC=LC-1;
                    isFirstRound=1;
                end
            else
                CB=[CB*(1+e) CB*(1-e)];
                e=eRed*e; 
                DT=DT*DTRed; 
                LC=size(CB,2);
                isFirstRound=1;
                Split=Split+1;
                IsThereABestCB=0; 
                disp(LC)
            end
        end
    end    
    if ~IsThereABestCB
        BestCB=CB;
        BestD=TotalDist;
        IsThereABestCB=1;
    else 
        if TotalDist<BestD
            BestCB=CB;
            BestD=TotalDist;
        end
    end
    EachSizeIterCounter=EachSizeIterCounter+1;
    if EachSizeIterCounter>maxIterInEachSize 
        EachSizeIterCounter=0;
        CB=BestCB; 
        IsThereABestCB=0;
        if LC>=L 
            if L==LC 
                disp(TotalDist)
                break
            else 
                [temp, Ind]=min(Population);
                NCB=zeros(d,LC-1);
                NCB=CB(:,setxor(1:LC,Ind(1)));
                CB=NCB;
                LC=LC-1;
                isFirstRound=1;
            end
        else 
            CB=[CB*(1+e) CB*(1-e)];
            e=eRed*e; 
            DT=DT*DTRed; 
            LC=size(CB,2);
            isFirstRound=1;
            Split=Split+1;
            IsThereABestCB=0;
            disp(LC)
        end
    end        
    disp(TotalDist)
    p=Population/N;
    save CBTemp CB p DistHist
end
m=CB;

p=Population/N;

disp(['Iterations = ' num2str(Iter)])
disp(['Split = ' num2str(Split)])

function [v, i]=maxn(x,n)


if nargin<2
    [v, i]=max(x); 
else
    n=min(length(x),n);
    [v, i]=sort(x);
    v=v(end:-1:end-n+1);
    i=i(end:-1:end-n+1);    
end
        
function [I, dst]=VQIndex(X,CB) 


L=size(CB,2);
N=size(X,2);
LNThreshold=64*10000;

if L*N<LNThreshold
    D=zeros(L,N);
    for i=1:L
        D(i,:)=sum((repmat(CB(:,i),1,N)-X).^2,1);
    end
    [dst I]=min(D);
else
    I=zeros(1,N);
    dst=I;
    for i=1:N
        D=sum((repmat(X(:,i),1,L)-CB).^2,1);
        [dst(i) I(i)]=min(D);
    end
end
    
function [I, dist]=VQLSFSpectralIndex(X,CB,W)


if nargin<3
    L=256;
    W=ones(L,1);
else
    if isscalar(W)
        L=W;
        W=ones(L,1);
    elseif isvector(W)
        W=W(:);
        L=length(W);
    else
        error('Argumento invalido!')
    end
end

NX=size(X,2);
NCB=size(CB,2);

AX=lsf2lpc(X);
ACB=lsf2lpc(CB);


D=zeros(NCB,1);

w=linspace(0,pi,L+1);
w=w(1:end-1);
N=size(AX,2)-1;
WFZ=zeros(N+1,L);
IMAGUNIT=sqrt(-1);
for k=0:N
    WFZ(k+1,:)=exp(IMAGUNIT*k*w);
end

SCB=zeros(L,NCB);
for i=1:NCB
    SCB(:,i)=(1./abs(ACB(i,:)*WFZ));
end

I=zeros(1,NX);
dist=zeros(1,NX);
for j=1:NX
    SX=(1./abs(AX(j,:)*WFZ))';    
    for i=1:NCB
        D(i)=sqrt(sum(((SX-SCB(:,i)).^2).*W));
    end
    [dist(j), I(j)]=min(D);
end









