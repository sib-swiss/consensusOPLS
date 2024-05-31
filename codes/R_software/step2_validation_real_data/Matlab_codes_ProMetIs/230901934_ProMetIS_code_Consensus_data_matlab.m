%%Import data set
%Ouverture du fichier de données
data = load('20230901934_step2_validation_data_ProMetIs_collection.mat');
data = data.data_ProMetIS;
%Ouverture du fichier de noms des colonnes des données
VarNames = load('20230901934_step2_validation_data_ProMetIs_VarNames.mat');
VarNames = VarNames.data_ProMetIS_VarNames;
%Ouverture du fichier du facteur Y
Y = load('20230901934_step2_validation_data_ProMetIs_Y_factor.mat');
Y = Y.Y;
%Ouverture du fichier des noms d'echantillons
Sample_names = load('20230901934_step2_validation_data_ProMetIs_SampleNames.mat');
Sample_names = Sample_names.SampleNames;
%Extraction des noms de blocs
Blocs_Names = fieldnames(data);

%%Data preprocessing
%Construction de la collection
collection_ProMetIs(1) = matrix2saisir(data.c18acquity_neg);
collection_ProMetIs(2) = matrix2saisir(data.c18acquity_pos);
collection_ProMetIs(3) = matrix2saisir(data.c18hypersil_pos);
collection_ProMetIs(4) = matrix2saisir(data.hilic_neg);
%Modification des noms de variables de la collection
collection_ProMetIs(1).i = VarNames.c18acquity_neg;
collection_ProMetIs(2).i = VarNames.c18acquity_pos;
collection_ProMetIs(3).i = VarNames.c18hypersil_pos;
collection_ProMetIs(4).i = VarNames.hilic_neg;
% Sauvegarder la structure en format .mat
save('20230901934_ProMetIs_collection_matlab.mat', 'collection_ProMetIs')

%%ConsensusOPLS model
%Definition des paramètres du modèle
LVsPred=1; % Number of predictive component(s)
LVsOrtho=1; % Maximum number of orthogonal components
CVfolds=size(collection_ProMetIs(1).d,1); % Number of cross-validation folds

collection = collection_ProMetIs;
Y = Y; 
A = LVsPred; 
maxOrtholvs = LVsOrtho;
nrcv = CVfolds; 
cvType = 'nfold';
modelType = 'da'; 
verbose = 0;

RVConsensusOPLSModelCV=RVConsensusOPLS(collection,Y,A,maxOrtholvs,nrcv,cvType,modelType,verbose,Sample_names, Blocs_Names);

%tElapsed =1.0354

disp('RVConsensusOPLS Results ');
disp(['R2: ' num2str(RVConsensusOPLSModelCV.koplsModel.R2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
%0.51319
disp(['Q2: ' num2str(RVConsensusOPLSModelCV.cv.Q2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
%-0.033051
disp(['DQ2: ' num2str(RVConsensusOPLSModelCV.cv.DQ2Yhat(RVConsensusOPLSModelCV.cv.OrthoLVsOptimalNum+1))]);
%-0.072917
disp(' ');
disp('Confusion Matrix:');
%    0.6154    0.3846
%    0.4231    0.5769
disp(RVConsensusOPLSModelCV.da.confusionMatrix);

VIP=MBVIP(collection, Y, RVConsensusOPLSModelCV);

%% Plot the results
% Consensus Score plot
figure; plot(RVConsensusOPLSModelCV.koplsModel.T(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.To(:,1),'k.');
title('ConsensusOPLS Score plot');
%axis([-0.5 0.5 -1 1]);
xlabel('Predictive component')
ylabel('Orthogonal component')
text(RVConsensusOPLSModelCV.koplsModel.T(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.To(:,1), ...
    Sample_names(:,1), ...
    'HorizontalAlignment', 'left','VerticalAlignment', ...
    'bottom','FontSize', 8)
saveas(gcf,'20230901934_step2_validation_data_ProMetIS_scoresPlot.png')
 
% Block contributions to the predictive component
figure; bar(RVConsensusOPLSModelCV.koplsModel.lambda(:,1));
set(gca,'xticklabel',Blocs_Names)
ylabel('Weight');
title('Block contributions to the predictive component');
saveas(gcf,'20230901934_step2_validation_data_ProMetIS_BlocContribToPredCompo.png')
 
% Block contributions to the first orthogonal component
figure; bar(RVConsensusOPLSModelCV.koplsModel.lambda(:,2));
set(gca,'xticklabel',Blocs_Names)
ylabel('Weight');
title('Block contributions to the first orthogonal component');
saveas(gcf,'20230901934_step2_validation_data_ProMetIS_BlocContribTo1stOrthoCompo.png')

% Block contributions predictive vs. orthogonal
figure; scatter(RVConsensusOPLSModelCV.koplsModel.lambda(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.lambda(:,2));
%axis([0 0.5 0 0.5]);
xlabel('Predictive component')
ylabel('Orthogonal component')
text(RVConsensusOPLSModelCV.koplsModel.lambda(:,1), ...
    RVConsensusOPLSModelCV.koplsModel.lambda(:,2), ...
    Blocs_Names, 'HorizontalAlignment','left','VerticalAlignment', ...
    'bottom','FontSize',8)
title('Block contributions');
saveas(gcf,'20230901934_step2_validation_data_ProMetIS_BlocContribPredVSOrthoCompo.png')

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
saveas(gcf,'20230901934_step2_validation_data_ProMetIS_LoadingPlot.png')


%% Permutations
[PermRes]=RVConsensusOPLSPerm(collection,Y,1000,LVsPred,LVsOrtho);
 
%% Build results matrices
scores = table(RVConsensusOPLSModelCV.koplsModel.T, RVConsensusOPLSModelCV.koplsModel.To, 'VariableNames', {'Pred' 'Ortho'});
writetable(scores, '20230901934_step2_validation_data_ProMetIS_scores.txt');
 
loadings_VIP_bloc1 = table(VarNames.Ceramids, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{1,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{1,2}, ...
    VIP{1,1}, ...
    'VariableNames', {'Ceramids' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc1, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_Ceramids.txt');
 
loadings_VIP_bloc2 = table(VarNames.Eicosanoids, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{2,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{2,2}, ...
    VIP{1,2}, ...
    'VariableNames', {'Eicosanoids' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc2, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_Eicosanoids.txt');
 
loadings_VIP_bloc3 = table(VarNames.NegLcMs, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{3,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{3,2}, ...
    VIP{1,3}, ...
    'VariableNames', {'NegLcMs' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc3, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_NegLcMs.txt');

loadings_VIP_bloc4 = table(VarNames.NegSfcMs, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{4,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{4,2}, ...
    VIP{1,4}, ...
    'VariableNames', {'NegSfcMs' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc4, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_NegSfcMs.txt');

loadings_VIP_bloc5 = table(VarNames.Nmr, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{5,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{5,2}, ...
    VIP{1,5}, ...
    'VariableNames', {'Nmr' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc5, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_Nmr.txt');

loadings_VIP_bloc6 = table(VarNames.PosLcMs, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{6,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{6,2}, ...
    VIP{1,6}, ...
    'VariableNames', {'PosLcMs' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc6, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_PosLcMs.txt');

loadings_VIP_bloc7 = table(VarNames.PosSfcMs, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{7,1}, ...
    RVConsensusOPLSModelCV.koplsModel.loadings{7,2}, ...
    VIP{1,7}, ...
    'VariableNames', {'PosSfcMs' 'Loadings_p' 'Loadings_o' 'VIP'});
writetable(loadings_VIP_bloc7, '20230901934_step2_validation_data_ProMetIS_loadings_VIP_PosSfcMs.txt');