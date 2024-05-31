function [PermRes] = RVConsensusOPLSPerm(collection,Y,nbruns,PredLVs,maxOrtholvs, Sample_names, Blocs_Names, plot)
% Y random permutations for model validity estimation
%
%Leave-One-Out CV
%
% X: independant data kernel
% Y: dependant data block
% nbruns: number of random permutations
% maxOrtholvs: maximum number of orthogonal LVs to compute
%
% The first model is the original one
%
% JB 2012-2018
%
% Reference:
% J. Boccard, D.N. Rutledge
% A consensus OPLS-DA strategy for multiblock Omics data fusion
% Analytica Chimica Acta, 769, 30-39, 2013


if nargin==5
    plot=1;
end

NbObs=size(collection(1).d,1);
% [long larg]=size(X);

modelCV=RVConsensusOPLS(collection,Y,PredLVs,maxOrtholvs,NbObs,'nfold','da',0, Sample_names, Blocs_Names);
PermRes.lvnum(1)=modelCV.cv.OrthoLVsOptimalNum+PredLVs;
PermRes.R2val(1)=modelCV.koplsModel.R2Yhat(end);
PermRes.DQ2val(1)=modelCV.cv.DQ2Yhat(PermRes.lvnum(1));
PermRes.Q2val(1)=modelCV.cv.Q2Yhat(PermRes.lvnum(1));
PermRes.PredAc(1)=modelCV.da.tot_sens;
PermRes.Y{1}=Y;
PermRes.RV(1)=RV_modified(Y,Y);

h = waitbar(0,'Random Permutations');
for i=1:nbruns
    %clc
    waitbar(i/nbruns,h, ['Y random permutations: ' num2str(nbruns)]);
    %disp(['Y random permutation run # ' num2str(i)]);
    Ys=Y(randperm(length(Y)),:);
        
    modelCV=RVConsensusOPLS(collection,Ys,PredLVs,maxOrtholvs,NbObs,'nfold','da',0, Sample_names, Blocs_Names);
    PermRes.lvnum(i+1)=modelCV.cv.OrthoLVsOptimalNum+PredLVs; 
    PermRes.R2val(i+1)=modelCV.koplsModel.R2Yhat(end);
    PermRes.DQ2val(i+1)=modelCV.cv.DQ2Yhat(PermRes.lvnum(i+1));
    PermRes.Q2val(i+1)=modelCV.cv.Q2Yhat(PermRes.lvnum(i+1));
    PermRes.PredAc(i+1)=modelCV.da.tot_sens;
    PermRes.Y{i+1}=Ys;
    PermRes.RV(i+1)=RV_modified(Y,Ys);
end

if plot
    figure; scatter(PermRes.RV, PermRes.R2val,'g')
    %PR2 = polyfix(PermRes.RV(2:end),PermRes.R2val(2:end),PermRes.RV(1),PermRes.R2val(1),1);
    %R2line=refline(PR2(1), PR2(2));
    %set(R2line, 'Color', 'g')
    hold on
    scatter(PermRes.RV, PermRes.Q2val,'b')
    %PQ2 = polyfix(PermRes.RV(2:end),PermRes.Q2val(2:end),PermRes.RV(1),PermRes.Q2val(1),1);
    %Q2line=refline(PQ2(1), PQ2(2));
    %set(Q2line, 'Color', 'b')
end

close(h);