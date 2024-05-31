function [modelCV]=RVConsensusOPLS(collection,Y,A,maxOrtholvs,nrcv,cvType,modelType,verbose,Sample_names, Blocs_Names)
%
%Consensus OPLS-DA with RV oefficients weighting and DQ2 computation for
%discriminant analysis
%JB 2012-2013
%
% modelCV=RVConsensusOPLS(collection,Y,1,10,100,'nfold','da',0);
%
% Reference:
% J. Boccard, D.N. Rutledge
% A consensus OPLS-DA strategy for multiblock Omics data fusion
% Analytica Chimica Acta, 769, 30-39, 2013

tStart = clock; 
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

%add for results comparables
Yc.X = koplsDummy(Yc.X);

for ta=1:ntable
	temp=koplsKernel(collection(ta).d,[],'p',1);
    %xnorm(ta)=norm(collection(ta).d,'fro')^2;
    xnorm(ta)=norm(temp,'fro');
	AMat{ta}=temp/xnorm(ta);
    RV(ta)=(RV_modified(AMat{ta},Yc.X)+1)/2;
	W_mat=W_mat+RV(ta)*AMat{ta};
	
	if(ta==1)
        temp_table = array2table(temp, 'RowNames', Sample_names, 'VariableNames', Sample_names);
		writetable(temp_table, "20230901934_ProMetIs_RVConsensusOPLS_00_temp.txt", 'Delimiter', '\t', 'WriteRowNames', true);	
	end
end
AMat_table = cell(1, 7);
for i = 1:numel(AMat)
    AMat_table_tmp = array2table(AMat{i}, 'RowNames', Sample_names, 'VariableNames', Sample_names);
    AMat_table{i} = AMat_table_tmp;
    filename = ['20230901934_ProMetIs_RVConsensusOPLS_02_Amat_', num2str(i), '.txt'];
    writetable(AMat_table_tmp, filename, 'Delimiter', '\t', 'WriteRowNames', true);
end
RV_table = array2table(RV, 'VariableNames', Blocs_Names);
writetable(RV_table, "20230901934_ProMetIs_RVConsensusOPLS_01_RV.txt", 'Delimiter', '\t');
%W_mat_table = array2table(W_mat, 'RowNames', Sample_names, 'VariableNames', Sample_names);
%writetable(W_mat_table, "20230901934_ProMetIs_RVConsensusOPLS_03_Wmat.txt", 'Delimiter', '\t', 'WriteRowNames', true);


% K = W_mat; A = A; oax = maxOrtholvs; nrcv = nrcv; cvType = cvType;
% preProcK = preProcK; preProcY = preProcY; cvFrac = cvFrac; modelType = modelType;
% verbose = verbose;
Y = Yc.X;
modelCV=ConsensusOPLSCV(W_mat,Y,A,maxOrtholvs,nrcv,cvType,preProcK,preProcY,cvFrac,modelType,verbose);

Ylarg=size(Y,2);

if strcmp(modelType,'da') % Search for the optimal model based on DQ2
    for i=0:maxOrtholvs
        for j=1:Ylarg
         [dqq(i+1,j),PRESSD(i+1,j)] = DQ2(modelCV.cv.AllYhat(:,Ylarg*i+j),Y(modelCV.cv.cvTestIndex,j)); % Compute DQ2 index
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

%writematrix(dqq, "RVConsensusOPLS_4_dqq.txt")

if modelCV.cv.OrthoLVsOptimalNum==0
    OrthoLVsNum=1;
else
    OrthoLVsNum=modelCV.cv.OrthoLVsOptimalNum;
end

% Recompute the optimal model
K = W_mat;Y = Y;A = A;nox = OrthoLVsNum;preProcK;preProcY;
modelCV.koplsModel=koplsModel(W_mat,Y,A,OrthoLVsNum,preProcK,preProcY);
%writematrix(W_mat, "RVConsensusOPLS_5_optModel_W_mat.txt")

writematrix(modelCV.koplsModel.Cp, "20230901934_ProMetIs_ConsensusOPLSCV_01_Cp.txt")
writematrix(modelCV.koplsModel.Sp, "20230901934_ProMetIs_ConsensusOPLSCV_02_Sp.txt")
writematrix(modelCV.koplsModel.Up, "20230901934_ProMetIs_ConsensusOPLSCV_03_Up.txt")

AllYhat_table = array2table(modelCV.cv.AllYhat, 'RowNames', Sample_names, 'VariableNames', ["1_p" "2_p" "1_po1" "2_po1"]);
writetable(AllYhat_table, '20230901934_ProMetIs_ConsensusOPLSCV_04_AllYhat.txt', 'Delimiter', '\t', 'WriteRowNames', true);

writematrix(modelCV.cv.Q2Yhat, "20230901934_ProMetIs_ConsensusOPLSCV_05_Q2Yhat.txt")
writematrix(modelCV.cv.DQ2Yhat, "20230901934_ProMetIs_ConsensusOPLSCV_06_DQ2Yhat.txt")

Tp_table = cell(1, numel(modelCV.koplsModel.Tp));
for i = 1:numel(modelCV.koplsModel.Tp)
    Tp_table_tmp = array2table(modelCV.koplsModel.Tp{i}, 'RowNames', Sample_names);
    Tp_table{i} = Tp_table_tmp;
    filename = ['20230901934_ProMetIs_ConsensusOPLSCV_07_Tp_', num2str(i), '.txt'];
    writetable(Tp_table_tmp, filename, 'Delimiter', '\t', 'WriteRowNames', true);
end

scores = [modelCV.koplsModel.T modelCV.koplsModel.To];
scores_table = array2table(scores, 'RowNames', Sample_names, 'VariableNames', ["p_1" "o_1"]);
writetable(scores_table, '20230901934_ProMetIs_ConsensusOPLSCV_08_scores.txt', 'Delimiter', '\t', 'WriteRowNames', true);

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

lambda_raw_table = array2table(modelCV.koplsModel.lambda_raw, 'RowNames', Blocs_Names, 'VariableNames', ["p_1" "o_1"]);
writetable(lambda_raw_table, '20230901934_ProMetIs_RVConsensusOPLS_04_lambdaRaw.txt', 'Delimiter', '\t', 'WriteRowNames', true);

for nb=1:size(lambda,2)
	lambda(:,nb)=lambda(:,nb)/sum(lambda(:,nb));
end
modelCV.koplsModel.lambda=lambda;

lambda_table = array2table(modelCV.koplsModel.lambda, 'RowNames', Blocs_Names, 'VariableNames', ["p_1" "o_1"]);
writetable(lambda_table, '20230901934_ProMetIs_RVConsensusOPLS_05_lambda.txt', 'Delimiter', '\t', 'WriteRowNames', true);

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

loadings_table = cell(1, numel(Blocs_Names));
for i = 1:numel(Blocs_Names)
    loadings_table_tmp = array2table([modelCV.koplsModel.loadings{i,1} modelCV.koplsModel.loadings{i,2}], 'VariableNames', ["p_1" "o_1"]);
    loadings_table{i} = loadings_table_tmp;
    filename = ['20230901934_ProMetIs_RVConsensusOPLS_06_loadings_', num2str(i), '.txt'];
    writetable(loadings_table{i}, filename, 'Delimiter', '\t');
end

tStop = clock;
tElapsed = etime(tStop, tStart)
disp('Done!')