#' RVConsensusOPLS
#' Consensus OPLS-DA with RV coefficients weighting and 
#' DQ2 computation for discriminant analysis
#' 
#' @param 
#' @param 
#' @param 
#' @param 
#'
#' @return 

#' @examples
#' 

RVConsensusOPLS <- function(data = collection,
                            Y,
                            1, 10, 100,
                            "nfold",
                            "da", 0){
  # Evaluate time ellapse
  execution_time <- system.time(
    
  )
  return(c(execution_time, model))
}

# ---------- MATLAB ---------- #
ntable=size(collection,2);
nrow=size(collection(1).d,1) ; 
W_mat=zeros(nrow,nrow);
preProcK='mc';
preProcY='mc';
cvFrac=0.75;

if strcmp(modelType,'reg')
Yc=koplsScale(Y,'mc','no');
else
  Yc.X=Y;
end

for ta=1:ntable
temp=koplsKernel(collection(ta).d,[],'p',1);
%xnorm(ta)=norm(collection(ta).d,'fro')^2;
xnorm(ta)=norm(temp,'fro');
AMat{ta}=temp/xnorm(ta);
RV(ta)=(RV_modified(AMat{ta},Yc.X)+1)/2;
W_mat=W_mat+RV(ta)*AMat{ta};
end

modelCV=ConsensusOPLSCV(W_mat,Y,A,maxOrtholvs,nrcv,cvType,preProcK,preProcY,cvFrac,modelType,verbose);

Ylarg=size(Y,2);

if strcmp(modelType,'da') % Search for the optimal model based on DQ2
for i=0:maxOrtholvs
for j=1:Ylarg
[dqq(i+1,j),PRESSD(i+1,j)] = DQ2(modelCV.cv.AllYhat(:,Ylarg*i+j),Y(:,j)); % Compute DQ2 index
end
end
dq2=mean(dqq,2);
index=A; %Minimum model size
while index<maxOrtholvs+A && dq2(index+1)-dq2(index)>0.01  % 1 percent DQ2 increase criterion
index=index+1;
end
modelCV.cv.DQ2Yhat=dq2;
modelCV.cv.OrthoLVsOptimalNum=index-A;

else % Search for the optimal model based on Q2Yhat
index=A; %Minimum model size
while index<maxOrtholvs+A && modelCV.cv.Q2Yhat(index+1)-modelCV.cv.Q2Yhat(index)>0.01  % 1 percent Q2 increase criterion
index=index+1;
end
modelCV.cv.OrthoLVsOptimalNum=index-A;
end

if modelCV.cv.OrthoLVsOptimalNum==0
OrthoLVsNum=1;
else
  OrthoLVsNum=modelCV.cv.OrthoLVsOptimalNum;
end

% Recompute the optimal model
modelCV.koplsModel=koplsModel(W_mat,Y,A,OrthoLVsNum,preProcK,preProcY);

% Adjust Yhat to the selected model size
modelCV.cv.Yhat=modelCV.cv.AllYhat(:,(Ylarg*A)+(OrthoLVsNum*A):(Ylarg*A)+(OrthoLVsNum*A)+Ylarg-1);

% Compute the blocks contributions for the selected model
for j=1:ntable
for k=1:A
lambda(j,k)=modelCV.koplsModel.T(:,k)'*AMat{j}*modelCV.koplsModel.T(:,k);
    end
	for l=1:OrthoLVsNum
        lambda(j,l+A)=modelCV.koplsModel.To(:,l)'*AMat{j}*modelCV.koplsModel.To(:,l);
end
end

modelCV.koplsModel.lambda_raw=lambda;

for nb=1:size(lambda,2)
lambda(:,nb)=lambda(:,nb)/sum(lambda(:,nb));
end
modelCV.koplsModel.lambda=lambda;

% Compute the loadings for the selected model size

for ta=1:ntable
for m=1:A
loadings{ta,m}=collection(ta).d'*modelCV.koplsModel.T(:,m)/(modelCV.koplsModel.T(:,m)'*modelCV.koplsModel.T(:,m));
end
for n=1:OrthoLVsNum
loadings{ta,n+m}=collection(ta).d'*modelCV.koplsModel.To(:,n)/(modelCV.koplsModel.To(:,n)'*modelCV.koplsModel.To(:,n));
end
end

modelCV.RV=RV;
modelCV.koplsModel.loadings=loadings;
%disp('Done!')
