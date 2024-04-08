%%Import data set
%Open data files
data = load('20230901934_step2_validation_data_ProMetIS_collection.mat');
data = data.data_ProMetIS;
%Open column names
VarNames = load('20230901934_step2_validation_data_ProMetIS_VarNames.mat');
VarNames = VarNames.data_ProMetIS_VarNames;
%Open response variable Y (factor)
Y = load('20230901934_step2_validation_Y_factor.mat');
Y = Y.Y;
%Open sample names
Sample_names = load('20230901934_step2_validation_SampleNames.mat');
Sample_names = Sample_names.SampleNames;
%Extract bloc names
Blocs_Names = fieldnames(data);

%%Data preprocessing
%Build matlab collection
collection_ProMetIs(1) = matrix2saisir(data.c18acquity_neg); 
collection_ProMetIs(2) = matrix2saisir(data.c18acquity_pos); 
collection_ProMetIs(3) = matrix2saisir(data.c18hypersil_pos); 
collection_ProMetIs(4) = matrix2saisir(data.hilic_neg); 
%Change variables names in the collection
collection_ProMetIs(1).i = VarNames.c18acquity_neg;
collection_ProMetIs(2).i = VarNames.c18acquity_pos;
collection_ProMetIs(3).i = VarNames.c18hypersil_pos;
collection_ProMetIs(4).i = VarNames.hilic_neg;
% Save collection in .mat file
save('20230901934_step2_validation_collection_matlab.mat', 'collection_ProMetIs')

%%ConsensusOPLS model
%Model parameters
LVsPred=1; % Number of predictive component(s)
LVsOrtho=1; % Maximum number of orthogonal components = only 1 in matlab
CVfolds=size(collection_ProMetIs(1).d,1); % Number of cross-validation folds

% Check inputs
collection = collection_ProMetIs;
Y = Y; % same format as in R
A = LVsPred; 
maxOrtholvs = LVsOrtho;
nrcv = CVfolds; 
cvType = 'nfold';
modelType = 'da'; 
verbose = 0;
Sample_names = Sample_names;
Blocs_Names = Blocs_Names;

%Run model
RVConsensusOPLSModelCV=RVConsensusOPLS_modif(collection,Y,A,maxOrtholvs,nrcv, cvType,modelType,verbose, Sample_names, Blocs_Names);

disp('RVConsensusOPLS Results ');
disp(['R2: ' num2str(RVConsensusOPLSModelCV.koplsModel.R2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
disp(['Q2: ' num2str(RVConsensusOPLSModelCV.cv.Q2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
disp(['DQ2: ' num2str(RVConsensusOPLSModelCV.cv.DQ2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
disp(' ');
disp('Confusion Matrix:');
disp(RVConsensusOPLSModelCV.da.confusionMatrix);

VIP=MBVIP(collection, Y, RVConsensusOPLSModelCV);

%% Plot the results
% Consensus Score plot
figure; plot(RVConsensusOPLSModelCV.koplsModel.T(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.To(:,1),'k.');
title('ConsensusOPLS Score plot');
xlabel('Predictive component')
ylabel('Orthogonal component')
text(RVConsensusOPLSModelCV.koplsModel.T(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.To(:,1), ...
    Sample_names(:,1), ...
    'HorizontalAlignment', 'left','VerticalAlignment', ...
    'bottom','FontSize', 8)
%saveas(gcf,'20230901934_step2_validation_scoresPlot.png')
 
% Block contributions to the predictive component
figure; bar(RVConsensusOPLSModelCV.koplsModel.lambda(:,1));
set(gca,'xticklabel',Blocs_Names)
ylabel('Weight');
title('Block contributions to the predictive component');
%saveas(gcf,'20230901934_step2_validation_BlocContribToPredCompo.png')
 
% Block contributions to the first orthogonal component
figure; bar(RVConsensusOPLSModelCV.koplsModel.lambda(:,2));
set(gca,'xticklabel',Blocs_Names)
ylabel('Weight');
title('Block contributions to the first orthogonal component');
%saveas(gcf,'20230901934_step2_validation_BlocContribTo1stOrthoCompo.png')

% Block contributions predictive vs. orthogonal
figure; scatter(RVConsensusOPLSModelCV.koplsModel.lambda(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.lambda(:,2));
xlabel('Predictive component')
ylabel('Orthogonal component')
text(RVConsensusOPLSModelCV.koplsModel.lambda(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.lambda(:,2), ...
    Blocs_Names, 'HorizontalAlignment','left','VerticalAlignment', ...
    'bottom','FontSize',8)
title('Block contributions');
%saveas(gcf,'20230901934_step2_validation_BlocContribPredVSOrthoCompo.png')

% Loading plots (one for each table)
figure;
title('ConsensusOPLS Loading plot');
xlabel('Predictive component');
ylabel('Orthogonal component');
hold all;
for i=1:length(collection)
   scatter(RVConsensusOPLSModelCV.koplsModel.loadings{i,1}, ...
       RVConsensusOPLSModelCV.koplsModel.loadings{i,2},'.', ...
       'DisplayName',[Blocs_Names{i}]);
end
legend(gca,'show')
%saveas(gcf,'20230901934_step2_validation_LoadingPlot.png')


%% Permutations
[PermRes]=RVConsensusOPLSPerm(collection,Y,1000,LVsPred,LVsOrtho);
 
%% Build results matrices
loadings_VIP_bloc1 = table(VarNames.c18acquity_neg, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{1,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{1,2}, ...
    VIP{1,1}, ...
    'VariableNames', {'c18acquity_neg' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc1, '20230901934_step2_validation_loadings_VIP_c18acquity_neg.txt');
 
loadings_VIP_bloc2 = table(VarNames.c18acquity_pos, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{2,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{2,2}, ...
    VIP{1,2}, ...
    'VariableNames', {'c18acquity_pos' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc2, '20230901934_step2_validation_loadings_VIP_c18acquity_pos.txt');
 
loadings_VIP_bloc3 = table(VarNames.c18hypersil_pos, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{3,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{3,2}, ...
    VIP{1,3}, ...
    'VariableNames', {'c18hypersil_pos' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc3, '20230901934_step2_validation_loadings_VIP_c18hypersil_pos.txt');

loadings_VIP_bloc4 = table(VarNames.hilic_neg, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{4,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{4,2}, ...
    VIP{1,4}, ...
    'VariableNames', {'hilic_neg' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc4, '20230901934_step2_validation_loadings_VIP_hilic_neg.txt');

