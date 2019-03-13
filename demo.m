%THIS IS A DEMO FOR IMAGE ENHANCEMENT USING MEMETIC ALGORITHM EMPLOYING A 
%HYBRID METAHEURISTIC TECHNIQUE USING DIFFERENTIAL EVOLUTION AND SIMULATED
%ANNEALING

% clc;
% clear all;
global BestGen truelabels UB HMS LB BestFitness; 
load Reuters_reduced.csv;
A=Reuters_reduced;
[row col]=size(A);
truelabels = A(:,1);
A = A(:,2:col);
[row col]=size(A);
 %A = (A-min(A(:))) ./ (max(A(:)-min(A(:))));
UB=11;          %6 event crime 
LB=1;
HMS=20;
%MEMETIC(image,truelabels);\
%MDHS;
   fprintf('Press 1 to do the HS t\n');
   fprintf('Press 2 to do the DHS\n');
     fprintf('Press 3 to do the HS_KM\n');
       fprintf('Press 4 to do the MDHS_Shrink\n');
   fprintf('Press 5 to do the KM\n');
   fprintf('Press 6 to do the  Adaptive DHS\n');
   fprintf('Press 7 to do the  Adaptive MDHS\n');
   fprintf('Press 8 to CMDHS_Interleaved\n');
      fprintf('Press 9 for a Fuzzy clustering \n');
      fprintf('Press 10 for a Fuzzy clustering \n');
      fprintf('Press 11 for a chaotic CMDHS \n');   
      fprintf('Press 12 for a  GBHS \n');  
 tic()
 x = input('Enter a number: ');
%            
switch  x
% 
    case 1
        HS (A,UB,LB,HMS, truelabels);
    case 2
       DHS(A,UB,LB,HMS, truelabels);
    case 3
      MHS_Interleaved(A,UB,LB,HMS, truelabels)     
    case 4
      MDHS_Interleaved(A,UB,LB,HMS, truelabels)
    case 5
      KM(A,UB,LB);
    case 6 
       Adaptive_DHS(A,UB,LB,HMS, truelabels)
    case 7
      Adaptive_MDHS_Interleaved(A,UB,LB,HMS, truelabels)
    case 8 
       CMDHS_Interleaved(A,UB,LB,HMS, truelabels)
    case 9 
       FuzzyDHS(A,UB,LB,HMS, truelabels)
    case 10
        ACMDHS(A,UB,LB,HMS, truelabels)
    case 11 
        CMDHS_Interleaved_cls(A,UB,LB,HMS, truelabels)
    case 12 
       GBHS(A,UB,LB,HMS, truelabels)
%      
%    
    otherwise
      warning('Unexpected local search input')
end 
%  truelabels;
% [C, label] = ind2cluster(BestGen);
% label;
% valid_errorate(label, truelabels);
% EVAL = Evaluate( label,truelabels',UB ); 

%        fprintf('the F-measure after 1000 iterations equals to: %f\n',EVAL);
%        fprintf('the F-measure after 1000 iterations equals to: %f\n',BestFitness);

toc();
%MaxF(ii)= EVAL;
%[EVAL1, EVAL2, precision, recall] = Ev( truelabels',label',UB );    
%         
        toc();
         