function [IndexOfFragmentPoint,Similarity] = SimilarityOfFile( filename )
tempBuffer=imread(filename);

Similarity(1)=0;
height=size(tempBuffer,1);
wid=size(tempBuffer,2);
for i=1:1:height-1
    Similarity(i)= fdiffMU_1D(tempBuffer(i,:,:),tempBuffer(i+1,:,:))/wid;
end

for i=1:1:height-2
    Similarity1(i)= abs(Similarity(i+1)-Similarity(i));
end



ma=max(Similarity1);
NumOfFragment=0;


IndexOfMax2(1)=1;
IndexOfMax1(1)=1;

LenOfMax2=0;
j=1;
Itarate=0;
while NumOfFragment<10 && Itarate<=3

    IndexOfMax2=IndexOfMax1;
   
    
    Itarate=Itarate+1;
    LenOfMax2=j;
    
    ma=ma/2;
    IndexOfMax1=find( Similarity1 > ma);
    
    NumOfFragment=1;
    len=size(IndexOfMax1,2);
    
    j=1
    IndexOfMax1(1)=IndexOfMax1(1);
    for i=1:1:len-1
        if IndexOfMax1(i+1)-IndexOfMax1(i)>10
            j=j+1;
            NumOfFragment=NumOfFragment+1;
            IndexOfMax1(j)=IndexOfMax1(i+1);
            
        end
    end
    
    
end

X=1:1:height-2 ;
Y=1:1:height-2 ;
%plot(X,Similarity(Y));
hold on;
plot(X,Similarity1(Y),'r');

Y=1:1:LenOfMax2;
plot(IndexOfMax2(Y),Similarity1(IndexOfMax2(Y)),'b*');


IndexOfFragmentPoint=IndexOfMax2(1:LenOfMax2);
for k=1:LenOfMax2
    if mod(IndexOfFragmentPoint(k),2)==0
        IndexOfFragmentPoint(k)=IndexOfFragmentPoint(k)-1;
    end
end


end

