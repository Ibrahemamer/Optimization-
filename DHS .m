%%DHS
function [BestGen,BestFitness,gx]=DHS(A,UB,LB,HMS, truelabels)

% This code has been written with Matlab 7.0
% You can modify the simple constraint handlening method using more efficient
% methods. A good review of these methods can be found in:
% Carlos A. Coello Coello,"Theoretical and numerical constraint-handling techniques used with evolutionary algorithms: a survey of the state of the art"

global NVAR NG NH MaxItr HMCR PARmin PARmax bwmin bwmax;
global HM NCHV fitness PVB BW ;
global BestIndex WorstIndex BestFit WorstFit BestGen currentIteration;

% [row col]=size(A);
% truelabels = A(:,1); 
%    A = A(:,2:col);
   [row col]=size(A);

NVAR=row;         %number of variables
NG=6;           %number of ineguality constraints
NH=0;           %number of eguality constraints
MaxItr=1000;    % maximum number of iterations
HMCR=0.6;       % harmony consideration rate  0< HMCR <1
PARmin=0.45;      % minumum pitch adjusting rate
PARmax=0.9;      % maximum pitch adjusting rate
bwmin=0.001;    % minumum bandwidth
bwmax=1.0;      % maxiumum bandwidth
PVB=zeros(NVAR,2);
PVB=[UB LB];   % range of variables
F=0.2;
% /**** Initiate Matrix ****/
HM=zeros(HMS,NVAR);
NCHV=zeros(1,NVAR);
BestGen=zeros(1,NVAR);
fitness=zeros(1,HMS);
BW=zeros(1,NVAR);
gx=zeros(1,NG);

error =zeros(20,1);
Max_F=zeros(20,1);
ADDC =zeros(20,1);
% INITIALIZATION and randeual and fetness function of cost

cost=zeros(1,1);
HM=fix(LB+((UB-LB).*rand(HMS,NVAR)));

%for i=1:HMS                        %for pop size
 %   cost(i,1)  = objective_function( UB,col,NVAR,A,xx,i);
%end

    function cost=Fitness(xx,i)
       cost  = objective_function( UB,col,NVAR,A,xx,i); 
    end

% warning off MATLAB:m_warning_end_without_block

MainHarmony;


% /*********************************************/

    function initialize
        % randomly initialize the HM
        for i=1:HMS
            
            
            fitness(i) = Fitness(HM,i);
        end
    end

%/*******************************************/

    function MainHarmony
        % global NVAR NG NH MaxItr HMS HMCR PARmin PARmax bwmin bwmax;
        % global HM NCHV fitness PVB BW gx currentIteration;
      for iii=1:20  
        initialize;
        currentIteration  = 0;
        
        while(StopCondition(currentIteration))
            
            PAR=(PARmax-PARmin)/(MaxItr)*currentIteration+PARmin;
            coef=log(bwmin/bwmax)/MaxItr;
            for pp =1:NVAR
                BW(pp)=bwmax*exp(coef*currentIteration);
            end
            % improvise a new harmony vector
            for i =1:NVAR
                ran = rand(1);
                if( ran < HMCR ) % memory consideration
                index1 = randint(1,HMS);
                index2 = randint(1,HMS);
                index3 = randint(1,HMS);
                     while (index1==index2||index1==index3||index2==index3)
                     index1 = randint(1,HMS);
                     index2 = randint(1,HMS);
                     index3 = randint(1,HMS);
                     break
                     end 
                      NCHV(i)=fix(HM(index1,i)+F*(HM(index2,i)-HM(index3,i)));
                       for iter= 1:row
                       if (NCHV(i) < LB ||NCHV(i)>UB)
                       NCHV(i)=fix(LB +(UB-LB).*rand);
                       end;
                       end;
                else
                    NCHV(i) = randval( PVB(1,1), PVB(1,2) ); % random selection
                end
            end
            
            newFitness = objective_function( UB,col,NVAR,A,NCHV,1);
            UpdateHM( newFitness,NCHV );
            
            currentIteration=currentIteration+1;
        end
          %hold off
        %BestFitness = min(fitness);
          BestFitness = min(fitness);
          ADDC(iii)=BestFitness;
         [C, label] = ind2cluster(BestGen);
       e= valid_errorate(label, truelabels);
       error(iii)=e;
        EVAL = Evaluate( label,truelabels',UB ); 
        Max_F(iii)=EVAL;
    end
        T = table(Max_F, ADDC, error,'VariableNames',{'Fmeasure','ADDC','Error'});
       %T = table( BFarray, 'VariableNames',{'ADDC'});
       % Step 3: Write in a file
       writetable(T,'C:\Users\32599642\Desktop\DHS.xls');
    end 
    
    
% /*****************************************/

    function UpdateHM( NewFit,X )
        % global NVAR MaxItr HMS ;
        % global HM NCHV BestGen fitness ;
        % global BestIndex WorstIndex BestFit WorstFit currentIteration;
        
        if(currentIteration==0)
            BestFit=fitness(1);
            for i = 1:HMS
                if( fitness(i) < BestFit )
                    BestFit = fitness(i);
                    BestIndex =i;
                end
            end
            
            WorstFit=fitness(1);
            WorstIndex =1;
            for i = 1:HMS
                if( fitness(i) > WorstFit )
                    WorstFit = fitness(i);
                    WorstIndex =i;
                end
            end
        end
        if (NewFit< WorstFit)
            
            if( NewFit < BestFit )
                HM(WorstIndex,:)=X;
                BestGen=X;
                fitness(WorstIndex)=NewFit;
                BestIndex=WorstIndex;
            else
                HM(WorstIndex,:)=X;
                fitness(WorstIndex)=NewFit;
            end
            
            
            WorstFit=fitness(1);
            WorstIndex =1;
            for i = 1:HMS
                if( fitness(i) > WorstFit )
                    WorstFit = fitness(i);
                    WorstIndex =i;
                end
            end
            
        end
        
    end % main if
[Fmin,index]=min(fitness);Xmin=HM(index,:);
BFitness = Fmin;
BestGen=Xmin;
%               
%    labelid=BestGen;
%     [C, label] = ind2cluster(labelid);
%     valid_errorate(label, truelabels);
% disp(['Best Fitnessn of HS: ', num2str(Fmin)])
% disp(['HS: ', num2str(BestGen)])
end %function

% /*****************************************/
function val1=randval(Maxv,Minv)
    val1=fix(rand(1)*(Maxv-Minv)+Minv);
end

function val2=randint(Maxv,Minv)
    val2=round(rand(1)*(Maxv-Minv)+Minv);
end
% /*******************************************/

function val=StopCondition(Itr)
    global MaxItr;
    val=1;
    if(Itr>=MaxItr)
        val=0;
    end
end

% /*******************************************/




