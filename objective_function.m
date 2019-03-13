function [cost] = objective_function( UB,col,nvars,A,x,i)
    costu=0;
    [row col]=size(A);
    cent=zeros(UB-1,col);
    suqcent=zeros(UB-1,col);
    calcent=zeros(UB-1,col);
    sumAcent=zeros(UB-1,1);
    costn=zeros(UB-1,col);
    suqa=zeros(row,UB-1);
    sumsuqA=zeros(row,1);
    d=zeros(1,1);
    t=zeros(1,1);
    e=zeros(1,1);
    b=zeros(1,1);
    r=zeros(1,1);
    suqA=A.^2;
    for sumA=1:row                                      %A.^2
        sumsuqA(sumA,1)=sum(suqA(sumA,:));
    end
    suqA=sqrt(sumsuqA);
    for ii=1:UB-1                                       %for centroids
       w=0;
        for iii=1:nvars                                 %for elements inside clusters
            
            if x(i,iii)==ii
                w=w+1;
                calcent(ii,:)=A(iii,:)+calcent(ii,:);   %for sumation centroids for each clusters  
                
            end
        end
        cent(ii,:)=calcent(ii,:)/w;                      %for calculate initial centroids for each clusters
    end
    suqcent=cent.^2;
    for sumcent=1:UB-1                                   %B.^2
        sumAcent(sumcent,1)=sum(suqcent(sumcent,:));
    end
    sumAcent=sqrt(sumAcent);
    for c=1:UB-1                                         %to clalculate cost
        ni=0;
        
        for cc=1:nvars                                   %to calculate each clusters
            
            if x(i,cc)==c
                
                 ni=ni+1;
                costn(c,:)=cent(c,:).*A(cc,:);          %A.*B
                suqa(cc,c)=sum(costn(c,:));             %sum A.*B
                e=suqA(cc,:).*sumAcent(c,:);            %sum A.^2.*B.^2
                b=suqa(cc,c)/e;                         %A.*b/A.^2.*B.^2
                d=b+d;                                  %sumA.*b/A.^2.*B.^2
            end
            
        end
        if ni==0
            ni=1;
            r=d/ni;                                         % (A.*b/A.^2.*B.^2)/ni
        elseif isnan(d)
            d=0;
            r=d/ni;
        else
            r=d/ni;
        end
        t=r+t;                                          %sum (A.*b/A.^2.*B.^2)/ni
        d=0;        
    end                           
     costu(1,1)=t/(UB-1);                                  %cost all clusters for first solution
     cost=1-costu;
end

