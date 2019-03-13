function  DE(image,truelabels)

% % -------------------------------------------------------------------------
% % original parametets 
% %pop_size=30; % set population size
cr = 0.2;   %the best out of 0.002, 0.9 , 0.02
F = 0.8;  % the best out of 0.1, 0.9 , 0.5, 

% -------------------------------------------------------------------------
runs=20;
BFarray = zeros (runs,1);
error= zeros (runs,1);
% FMacro=zeros(runs,1);
% FMicro=zeros(runs,1);
% MaxP=zeros(runs,1);
% MaxR=zeros(runs,1);
MaxF=zeros(runs,1);
[row col]=size(image);
UB=3;
LB=1;
a_hi = UB; %boundary values of variables (gene)
a_lo = LB;
%paramters for Differential Evolution

% -------------------------------------------------------------------------
 pop_size=row; % set population size

for ii=1:20 %iteration counter initialized
%GENERATE INITIAL POPULATION
a = fix((a_hi-a_lo).*rand(pop_size,10)+a_lo); 
par = a;
[ro co]=size(par);
for iter= 1:co
costassignment(iter,:)=par(:,iter);
end;

for iter= 1:co
fitness(iter,1)=objective_function( UB,col,row,image, costassignment(iter,:),1);
%fitness(iter,1)=Fuzzyobjective_function( UB,col,row,image, costassignment(iter,:),1);
 %fitness(iter,1)=Euclidian( UB,col,row,image, costassignment(iter,:),1);
end;
 
%for iu=1:30 
 for j=1:co
       % DE STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Pick 3 parents
       r(1) = floor(rand()* co) + 1;
       while (r(1)==j)
	     r(1) = floor(rand()* co) + 1;
       end
       r(2) = floor(rand()* co) + 1;
       while ((r(2)==r(1))||(r(2)==j))
         r(2) = floor(rand()* co) + 1;
       end
       r(3) = floor(rand()* co) + 1;
       while ((r(3)==r(2))||(r(3)==r(1))||(r(3)==j))
         r(3) = floor(rand()* co) + 1;
       end
       %mutant vector
       v=fix(par(:,r(1))+F*(par(:,r(2))-par(:,r(3))));
        for iter= 1:row
            if (v(iter,1) < a_lo ||v(iter,1)>a_hi)
                v(iter,1)=fix(a_lo +(a_hi-a_lo).*rand);
            end;
        end;
       
     
%recombination to generate offspring
      u=par(:,j);
       rv=rand(1,row);
       for i=1:row
           if rv(1,i)<cr        
            u(i,1)=v(i,:);
           end
       end
  
 %recompute fitness value of offsprings 
 costassignment(1,:)=u;
   new_fitness=objective_function( UB,col,row,image, costassignment(1,:),1);
         %new_fitness=Fuzzyobjective_function( UB,col,row,image, costassignment(1,:),1);
     % new_fitness= Euclidian( UB,col,row,image, costassignment(1,:),1);
%        new_fitness =- tr_op(global_mean,B,C,im_size,image,u(1,1),u(1,2),u(1,3),u(1,4));
     if(new_fitness <= fitness(j))
         par(:,j)=u;
          fitness(j)=new_fitness;
     end
 end
%end

 [Fmin,index]=min(fitness);
      Xmin= par(:,index );
       BFitness = Fmin;
       BestGen=Xmin; 
% fprintf('Best ADDC: %4.2f %% \n',BFitness);   
BFarray(ii) = BFitness;
%       hold on
%       pause(0.01)
%             
%             figure(1)
%             title('Memetic ADDC Cluster')
%             xlabel('Iterations')
%             ylabel('function value of ADDC')
%             plot(iu,BFitness,'--r.','LineWidth',2,...
%                 'MarkerEdgeColor','k')
%  
% for iter= 1:co
%     up_par(iter,:)= par(:,iter);
% end;  
%  
%  [Fmin,index]=min(fitness);
%       Xmin=up_par(index,:);
%        BFitness = Fmin;
%        BestGen=Xmin; 
% fprintf('Best ADDC: %4.2f %% \n',BFitness);   
% BFarray(ui) = BFitness;

[C, label] = ind2cluster(BestGen);
Rerror = valid_errorate(label, truelabels);
error(ii)=Rerror;
%truedata(1,:)= transpose(truelabels);       
EVAL = Evaluate( truelabels',label',UB );
MaxF(ii)= EVAL;
% [EVAL1, EVAL2, precision, recall] = Ev( truelabels',label',UB );     
% FMacro(ii)=EVAL1;
% FMicro(ii)=EVAL2;
% MaxP(ii)=precision;
% MaxR(ii)=recall;
end

%Step 2: Generate a table of data
T = table(BFarray, MaxF,error, 'VariableNames',{'ADDC','FMeasure','Errorratio'});

%T = table(BFarray, FMacro, FMicro,MaxP,MaxR,error,'VariableNames',{'ADDC','FMacro','FMicro','Pmacro','Rmacro','error'});

% Step 3: Write in a file
writetable(T,'C:\Users\32599642\Desktop\exceldata.xls');

end 


   
