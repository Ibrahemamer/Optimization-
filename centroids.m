function [cent] = centroids( UB,col,nvars,A,x)
    [row col]=size(A);
    cent=zeros(UB-1,col);
    
    calcent=zeros(UB-1,col);
    
    for ii=1:UB-1               %for centroids
       w=0;
        for iii=1:nvars          %for elements inside clusters
            
            if x(1,iii)==ii
                w=w+1;
                calcent(ii,:)=A(iii,:)+calcent(ii,:);   %for calculate centroids for each clusters  
                
            end
        end
        cent(ii,:)=calcent(ii,:)/w;                      %for calculate centroids for each clusters
    end
