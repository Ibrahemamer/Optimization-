 function [FMarco,Fmicro ,pp, r] = Evaluate( truedata,labelt,UB )
% This fucntion evaluates the performance of a classification model by 
% calculating the common performance measures: Accuracy, Sensitivity, 
% Specificity, Precision, Recall, F-Measure, G-mean.
% Input: ACTUAL = Column matrix with actual class labels of the training
%                 examples
%        PREDICTED = Column matrix with predicted class labels by the
%                    classification model
% Output: EVAL = Row matrix with all the performance measures
pp=0;
r=0;
f_measure=0;
precision=0;
alpha=length(truedata);
[truedata] = sort(truedata);
[labelt] = sort(labelt);
for i=1:UB-1 
idx = (truedata()==i);

p = length(truedata (:,idx));
nn = labelt(:,idx)==i;
if nn==0;
    n=0;
else
n=sum(nn~=0);
end
N = p+n;
%[truelabels, S] = sort(truelabels);
tp = sum(truedata(idx)==labelt(idx));
tn = sum(truedata(~idx)==labelt(~idx));
fp = n-tn;
fn = p-tp;

tp_rate = tp/p;
tn_rate = tn/n;

%accuracy = ((tp+tn)/(tp+tn+fp+fn))*100;
%c=c+accuracy;
%sensitivity = tp_rate *100;
%s=s+sensitivity;
%specificity = tn_rate*100;
%sc=sc+specificity;

precision =  (tp/(n));
 if isnan(precision)
    precision=0;
end
pp=pp+precision;

recall = tp/(p);
 if isnan(recall)
    recall=0;
end
r=r+recall;

f1=(2*precision*recall)/(precision + recall);
 if isnan(f1)
    f1=0;
end 
f_measure =(f_measure+ f1); 
end 
FMarco=(f_measure/(UB-1));
pp=pp/(UB-1);
r= r/(UB-1);
Fmicro= (2*pp*r)/(pp + r);
fprintf('P_Macro: %4.4f %% \n', pp);
fprintf('R_Macro: %4.4f %% \n', r);
fprintf('FMacro: %4.4f %% \n', FMarco);
fprintf('F_Micro: %4.4f %% \n', Fmicro);

 end 
