% %% Load Zeisel data, paste feature importance data into cell and reshape
% clear; clc;
% matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles';
% figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
% % Creating Zeisel_FeatureImportance.mat
% % % paste in Excel data from Harry 
% % zeisel_fi_all = reshape(zeisel_fi_all,200,20);
% % zeisel_fi_all = zeisel_fi_all(:,11:end);
% % zeisel_fi_all = cell2mat(zeisel_fi_all);
% % % paste in Excel data from Harry
% % zeisel_fi_long = reshape(zeisel_fi_long,200,20);
% % zeisel_fi_long = zeisel_fi_long(:,11:end);
% % zeisel_fi_long = cell2mat(zeisel_fi_long);
% % % paste in Excel data from Harry
% % zeisel_fi_short = reshape(zeisel_fi_short,200,20);
% % zeisel_fi_short = zeisel_fi_short(:,11:end);
% % zeisel_fi_short = cell2mat(zeisel_fi_short);
% % save('Zeisel_FeatureImportance.mƒat', ...
% %     'zeisel_fi_all','zeisel_fi_short','zeisel_fi_long');
% 
% % Average cell types over broader classes in Zeisel ontology
% % zeisel_classes = unique(zeisel_ontology_cns(:,2));
% % zeisel_fi_all_class = zeros(length(zeisel_classes),size(zeisel_fi_all,2));
% % zeisel_fi_long_class = zeisel_fi_all_class;
% % zeisel_fi_short_class = zeisel_fi_all_class;
% % zeisel_fi_out_class = zeisel_fi_all_class;
% % zeisel_fi_in_class = zeisel_fi_all_class;
% % 
% % for i = 1:length(zeisel_classes)
% %     classinds = ismember(zeisel_ontology_cns(:,2),zeisel_classes{i});
% %     zeisel_fi_all_class(i,:) = mean(zeisel_fi_all(classinds,:));
% %     zeisel_fi_short_class(i,:) = mean(zeisel_fi_short(classinds,:));
% %     zeisel_fi_long_class(i,:) = mean(zeisel_fi_long(classinds,:));
% %     zeisel_fi_out_class(i,:) = mean(zeisel_fi_out(classinds,:));
% %     zeisel_fi_in_class(i,:) = mean(zeisel_fi_in(classinds,:));
% % end
% % % Manually tweak zeisel_classes here
% zeisel_fi_all_class_mean = mean(zeisel_fi_all_class,2);
% [~,zeisel_fi_all_reorderinds] = sort(zeisel_fi_all_class_mean,'descend');
% zeisel_fi_all_class_reorder = zeisel_fi_all_class(zeisel_fi_all_reorderinds,:);
% zeisel_classes_all_reorder = zeisel_classes(zeisel_fi_all_reorderinds);
% 
% zeisel_fi_short_class_mean = mean(zeisel_fi_short_class,2);
% [~,zeisel_fi_short_reorderinds] = sort(zeisel_fi_short_class_mean,'descend');
% zeisel_fi_short_class_reorder = zeisel_fi_short_class(zeisel_fi_short_reorderinds,:);
% zeisel_classes_short_reorder = zeisel_classes(zeisel_fi_short_reorderinds);
% 
% zeisel_fi_long_class_mean = mean(zeisel_fi_long_class,2);
% [~,zeisel_fi_long_reorderinds] = sort(zeisel_fi_long_class_mean,'descend');
% zeisel_fi_long_class_reorder = zeisel_fi_long_class(zeisel_fi_long_reorderinds,:);
% zeisel_classes_long_reorder = zeisel_classes(zeisel_fi_long_reorderinds);
% 
% zeisel_fi_out_class_mean = mean(zeisel_fi_out_class,2);
% [~,zeisel_fi_out_reorderinds] = sort(zeisel_fi_out_class_mean,'descend');
% zeisel_fi_out_class_reorder = zeisel_fi_out_class(zeisel_fi_out_reorderinds,:);
% zeisel_classes_out_reorder = zeisel_classes(zeisel_fi_out_reorderinds);
% 
% zeisel_fi_in_class_mean = mean(zeisel_fi_in_class,2);
% [~,zeisel_fi_in_reorderinds] = sort(zeisel_fi_in_class_mean,'descend');
% zeisel_fi_in_class_reorder = zeisel_fi_in_class(zeisel_fi_in_reorderinds,:);
% zeisel_classes_in_reorder = zeisel_classes(zeisel_fi_in_reorderinds);
% 
% 
% %% Load Tasic data, paste feature importance data into cell and reshape

% % % paste in Excel data from Harry
% % tasic_fi_all = reshape(tasic_fi_all,25,20);
% % tasic_fi_all = tasic_fi_all(:,11:end);
% % tasic_fi_all = cell2mat(tasic_fi_all);
% % % paste in Excel data from Harry
% % tasic_fi_long = reshape(tasic_fi_long,25,20);
% % tasic_fi_long = tasic_fi_long(:,11:end);
% % tasic_fi_long = cell2mat(tasic_fi_long);
% % % paste in Excel data from Harry
% % tasic_fi_short = reshape(tasic_fi_short,25,20);
% % tasic_fi_short = tasic_fi_short(:,11:end);
% % tasic_fi_short = cell2mat(tasic_fi_short);
% % % paste in Excel data from Harry
% % tasic_fi_classification = reshape(tasic_fi_classification,25,20);
% % tasic_fi_classification = tasic_fi_classification(:,11:end);
% % tasic_fi_classification = cell2mat(tasic_fi_classification);
% % % 
% % save([matdir filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_all',...
% %     'tasic_fi_short','tasic_fi_long','tasic_fi_classification');
% load([matdir filesep 'Tasic_FeatureImportance.mat'])
% 
% % Average cell types over broader classes in Tasic ontology
% tasic_classes = unique(tasic_ontology_cns(:,2));
% tasic_fi_all_class = zeros(length(tasic_classes),size(tasic_fi_all,2));
% tasic_fi_long_class = tasic_fi_all_class;
% tasic_fi_short_class = tasic_fi_all_class;
% tasic_fi_classification_class = tasic_fi_all_class;
% tasic_fi_out_class = tasic_fi_all_class;
% tasic_fi_in_class = tasic_fi_all_class;
% for i = 1:length(tasic_classes)
%     classinds = ismember(tasic_ontology_cns(:,2),tasic_classes{i});
%     tasic_fi_all_class(i,:) = mean(tasic_fi_all(classinds,:));
%     tasic_fi_short_class(i,:) = mean(tasic_fi_short(classinds,:));
%     tasic_fi_long_class(i,:) = mean(tasic_fi_long(classinds,:));
%     tasic_fi_classification_class(i,:) = mean(tasic_fi_classification(classinds,:));
%     tasic_fi_out_class(i,:) = mean(tasic_fi_out(classinds,:));
%     tasic_fi_in_class(i,:) = mean(tasic_fi_in(classinds,:));
% end
% 
% tasic_fi_all_class_mean = mean(tasic_fi_all_class,2);
% [~,tasic_fi_classification_reorderinds] = sort(tasic_fi_all_class_mean,'descend');
% tasic_fi_classification_class_reorder = tasic_fi_all_class(tasic_fi_classification_reorderinds,:);
% tasic_classes_classification_reorder = tasic_classes(tasic_fi_classification_reorderinds);
% 
% tasic_fi_short_class_mean = mean(tasic_fi_short_class,2);
% [~,tasic_fi_short_reorderinds] = sort(tasic_fi_short_class_mean,'descend');
% tasic_fi_short_class_reorder = tasic_fi_short_class(tasic_fi_short_reorderinds,:);
% tasic_classes_short_reorder = tasic_classes(tasic_fi_short_reorderinds);
% 
% tasic_fi_long_class_mean = mean(tasic_fi_long_class,2);
% [~,tasic_fi_long_reorderinds] = sort(tasic_fi_long_class_mean,'descend');
% tasic_fi_long_class_reorder = tasic_fi_long_class(tasic_fi_long_reorderinds,:);
% tasic_classes_long_reorder = tasic_classes(tasic_fi_long_reorderinds);
% 
% tasic_fi_classification_class_mean = mean(tasic_fi_classification_class,2);
% [~,tasic_fi_classification_reorderinds] = sort(tasic_fi_classification_class_mean,'descend');
% tasic_fi_classification_class_reorder = tasic_fi_classification_class(tasic_fi_classification_reorderinds,:);
% tasic_classes_classification_reorder = tasic_classes(tasic_fi_classification_reorderinds);
% 
% tasic_fi_out_class_mean = mean(tasic_fi_out_class,2);
% [~,tasic_fi_out_reorderinds] = sort(tasic_fi_out_class_mean,'descend');
% tasic_fi_out_class_reorder = tasic_fi_out_class(tasic_fi_out_reorderinds,:);
% tasic_classes_out_reorder = tasic_classes(tasic_fi_out_reorderinds);
% 
% tasic_fi_in_class_mean = mean(tasic_fi_in_class,2);
% [~,tasic_fi_in_reorderinds] = sort(tasic_fi_in_class_mean,'descend');
% tasic_fi_in_class_reorder = tasic_fi_in_class(tasic_fi_in_reorderinds,:);
% tasic_classes_in_reorder = tasic_classes(tasic_fi_in_reorderinds);

%% Boxplot results
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
savenclose = 1;
FIBoxplot(matdir,'Both','All',savenclose,figdir); % All connections
FIBoxplot(matdir,'Both','Classification',savenclose,figdir); % All classification
FIBoxplot(matdir,'Tasic',{'Short','Long'},savenclose,figdir); % Tasic, short/long
FIBoxplot(matdir,'Tasic',{'Out','In'},savenclose,figdir); % Tasic, source/target
FIBoxplot(matdir,'Tasic',{'Out Long','In Long'},savenclose,figdir); % Tasic, source/target, long
FIBoxplot(matdir,'Tasic',{'Out Short','In Short'},savenclose,figdir); % Tasic, source/target, short
FIBoxplot(matdir,'Zeisel',{'Short','Long'},savenclose,figdir); % Zeisel, short/long
FIBoxplot(matdir,'Zeisel',{'Out','In'},savenclose,figdir); % Zeisel, source/target
FIBoxplot(matdir,'Zeisel',{'Out Long','In Long'},savenclose,figdir); % Zeisel, source/target, long
FIBoxplot(matdir,'Zeisel',{'Out Short','In Short'},savenclose,figdir); % Zeisel, source/target, short
FIBoxplot(matdir,'Zeisel',{'NeoToAll Out','NeoToAll In'},savenclose,figdir); % Zeisel, source/target, neo-to/from-all

%% Plot long vs. short FI scatterplot (all non-neuronals)
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
load([matdir filesep 'Tasic_Ontology.mat'],'tasic_ontology_cns');
load([matdir filesep 'Zeisel_Ontology.mat'],'zeisel_ontology_cns');
load([matdir filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_short','tasic_fi_long')
load([matdir filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_short','zeisel_fi_long')
savenclose = 1;

tasic_fi_short_mean = mean(tasic_fi_short,2);
tasic_fi_short_mean_norm = (tasic_fi_short_mean - min(tasic_fi_short_mean))...
    /(max(tasic_fi_short_mean) - min(tasic_fi_short_mean));
tasic_fi_long_mean = mean(tasic_fi_long,2);
tasic_fi_long_mean_norm = (tasic_fi_long_mean - min(tasic_fi_long_mean))...
    /(max(tasic_fi_long_mean) - min(tasic_fi_long_mean));

zeisel_fi_short_mean = mean(zeisel_fi_short,2);
zeisel_fi_short_mean_norm = (zeisel_fi_short_mean - min(zeisel_fi_short_mean))...
    /(max(zeisel_fi_short_mean) - min(zeisel_fi_short_mean));
zeisel_fi_long_mean = mean(zeisel_fi_long,2);
zeisel_fi_long_mean_norm = (zeisel_fi_long_mean - min(zeisel_fi_long_mean))...
    /(max(zeisel_fi_long_mean) - min(zeisel_fi_long_mean));

teglu_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Neo Glu')).';
oligo_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Oligo')).';
vasc_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Vasc')).';
immune_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Immune')).';
astro_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Astro')).';
other_ind_tasic = setdiff(1:length(tasic_ontology_cns(:,2)),[teglu_ind_tasic,...
    oligo_ind_tasic,vasc_ind_tasic,astro_ind_tasic,immune_ind_tasic]);
teglu_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Hip Neo Glu')).';
oligo_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Oligo')).';
vasc_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Vasc')).';
immune_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Immune')).';
astro_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Astro')).';
other_ind_zeisel = setdiff(1:length(zeisel_ontology_cns(:,2)),[teglu_ind_zeisel,...
    oligo_ind_zeisel,vasc_ind_zeisel,astro_ind_zeisel,immune_ind_zeisel]);

scatter_vals = cell(2,6);
supertypes_zeisel = unique(zeisel_ontology_cns(:,2));
cmap_raw = hsv(length(supertypes_zeisel));
c_neoglu = cmap_raw(ismember(supertypes_zeisel,'Hip Neo Glu'),:);
c_oligo = cmap_raw(ismember(supertypes_zeisel,'Oligo'),:);
c_vasc = cmap_raw(ismember(supertypes_zeisel,'Vasc'),:);
c_immune = cmap_raw(ismember(supertypes_zeisel,'Immune'),:);
c_astro = cmap_raw(ismember(supertypes_zeisel,'Astro'),:);
scatter_c = [c_neoglu; c_oligo; c_vasc; c_immune; [0 0 1]; [0.85 0.85 0.85]];
scatter_s = {'o','o'};

scatter_vals{1,1} = [tasic_fi_short_mean_norm(teglu_ind_tasic),...
    tasic_fi_long_mean_norm(teglu_ind_tasic)];
scatter_vals{1,2} = [tasic_fi_short_mean_norm(oligo_ind_tasic),...
    tasic_fi_long_mean_norm(oligo_ind_tasic)];
scatter_vals{1,3} = [tasic_fi_short_mean_norm(vasc_ind_tasic),...
    tasic_fi_long_mean_norm(vasc_ind_tasic)];
scatter_vals{1,4} = [tasic_fi_short_mean_norm(immune_ind_tasic),...
    tasic_fi_long_mean_norm(immune_ind_tasic)];
scatter_vals{1,5} = [tasic_fi_short_mean_norm(astro_ind_tasic),...
    tasic_fi_long_mean_norm(astro_ind_tasic)];
scatter_vals{1,6} = [tasic_fi_short_mean_norm(other_ind_tasic),...
    tasic_fi_long_mean_norm(other_ind_tasic)];
scatter_vals{2,1} = [zeisel_fi_short_mean_norm(teglu_ind_zeisel),...
    zeisel_fi_long_mean_norm(teglu_ind_zeisel)];
scatter_vals{2,2} = [zeisel_fi_short_mean_norm(oligo_ind_zeisel),...
    zeisel_fi_long_mean_norm(oligo_ind_zeisel)];
scatter_vals{2,3} = [zeisel_fi_short_mean_norm(vasc_ind_zeisel),...
    zeisel_fi_long_mean_norm(vasc_ind_zeisel)];
scatter_vals{2,4} = [zeisel_fi_short_mean_norm(immune_ind_zeisel),...
    zeisel_fi_long_mean_norm(immune_ind_zeisel)];
scatter_vals{2,5} = [zeisel_fi_short_mean_norm(astro_ind_zeisel),...
    zeisel_fi_long_mean_norm(astro_ind_zeisel)];
scatter_vals{2,6} = [zeisel_fi_short_mean_norm(other_ind_zeisel),...
    zeisel_fi_long_mean_norm(other_ind_zeisel)];

figure('Position',[0 0 325 450]); hold on;
s1 = scatter(scatter_vals{1,6}(:,1),scatter_vals{1,6}(:,2),75,scatter_s{1},'filled',...
    'MarkerEdgeColor',scatter_c(6,:),'MarkerFaceColor',scatter_c(6,:));
scatter(scatter_vals{2,6}(:,1),scatter_vals{2,6}(:,2),75,scatter_s{2},'filled',...
    'MarkerEdgeColor',scatter_c(6,:),'MarkerFaceColor',scatter_c(6,:));
s2 = scatter(scatter_vals{1,1}(:,1),scatter_vals{1,1}(:,2),75,scatter_s{1},'filled',...
    'MarkerEdgeColor',scatter_c(1,:),'MarkerFaceColor',scatter_c(1,:));
s3 = scatter(scatter_vals{1,2}(1),scatter_vals{1,2}(2),75,scatter_s{1},'filled',...
    'MarkerEdgeColor',scatter_c(2,:),'MarkerFaceColor',scatter_c(2,:));
s4 = scatter(scatter_vals{1,3}(1),scatter_vals{1,3}(2),75,scatter_s{1},'filled',...
    'MarkerEdgeColor',scatter_c(3,:),'MarkerFaceColor',scatter_c(3,:));
s5 = scatter(scatter_vals{1,4}(1),scatter_vals{1,4}(2),75,scatter_s{1},'filled',...
    'MarkerEdgeColor',scatter_c(4,:),'MarkerFaceColor',scatter_c(4,:));
s6 = scatter(scatter_vals{1,5}(1),scatter_vals{1,5}(2),75,scatter_s{1},'filled',...
    'MarkerEdgeColor',scatter_c(5,:),'MarkerFaceColor',scatter_c(5,:));
scatter(scatter_vals{2,1}(:,1),scatter_vals{2,1}(:,2),75,scatter_s{2},'filled',...
    'MarkerEdgeColor',scatter_c(1,:),'MarkerFaceColor',scatter_c(1,:));
scatter(scatter_vals{2,2}(:,1),scatter_vals{2,2}(:,2),75,scatter_s{2},'filled',...
    'MarkerEdgeColor',scatter_c(2,:),'MarkerFaceColor',scatter_c(2,:));
scatter(scatter_vals{2,3}(:,1),scatter_vals{2,3}(:,2),75,scatter_s{2},'filled',...
    'MarkerEdgeColor',scatter_c(3,:),'MarkerFaceColor',scatter_c(3,:));
scatter(scatter_vals{2,4}(:,1),scatter_vals{2,4}(:,2),75,scatter_s{2},'filled',...
    'MarkerEdgeColor',scatter_c(4,:),'MarkerFaceColor',scatter_c(4,:));
scatter(scatter_vals{2,5}(:,1),scatter_vals{2,5}(:,2),75,scatter_s{2},'filled',...
    'MarkerEdgeColor',scatter_c(5,:),'MarkerFaceColor',scatter_c(5,:));
ylim([0 1]); yticks(1); xlim([0 1]); xticks([0 1]);
ylabel('Long-range F.I.'); xlabel('Short-range F.I.'); 
set(gca,'box','on','FontName','DejaVunSans','FontSize',16);
l = legend([s2, s3, s4, s5, s6, s1], {'Tel Glu','Oligo','Vasc','Immune','Astro','Other'},...
    'box','off','Position',[0.6 0.675 0.25 0.25],'FontSize',14);
% text(0.05,0.925,'Tel Glu','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(1,:));
% text(0.65,0.5,'Oligo','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(2,:));
% text(0.9,0.2,'Vasc','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(3,:));

% testshort = []; testlong = [];
% for i = 1:2
%     for j = 1:4
% %         if size(scatter_vals{i,j},1) == 1
% %             testshort = [testshort; scatter_vals{i,j}(1)];
% %             testlong = [testlong; scatter_vals{i,j}(2)];
% %         else
%             testshort = [testshort; scatter_vals{i,j}(:,1)];
%             testlong = [testlong; scatter_vals{i,j}(:,2)];
% %         end
%     end
% end
% [r,p] = corr(testshort,testlong)
if savenclose
    print([figdir filesep 'FI_AllCTScatter'],'-dtiffn','-r300'); close;
end

%% Plot eigenvector centrality for all types, Tasic
load([matdir filesep 'Centrality_Results.mat'],'centrality_results');
load([matdir filesep 'Tasic_outstruct_ng606.mat'],'outstruct');
types = tasic_classes;
degree_cell = cell(2,length(types));
degree_res = centrality_results(ismember(centrality_results(:,2),'degree'),:);
eigen_cell = cell(2,length(types));
eigen_res = centrality_results(ismember(centrality_results(:,2),'eigenvector_centrality'),:);
degree_C = cell2mat(degree_res(1:424,4)); eigen_C = cell2mat(eigen_res(1:424,4));

for i = 1:length(types)
    ctinds = ismember(tasic_ontology_cns(:,2), types{i});
    celldensity_typei = outstruct(1).corrB(:,ctinds);   
    celldensity_typei_reg = zeros(426,sum(ctinds));
    for j = 1:sum(ctinds)
        [~,celldensity_typei_regj] = Voxel_To_Region(celldensity_typei(:,j),matdir);
        celldensity_typei_reg(:,j) = celldensity_typei_regj;
    end
    celldensity_typei_reg([12,(12+213)],:) = [];
    [~,scores_typei,~] = pca(celldensity_typei_reg);
    typei_score1 = scores_typei(:,1);   
    typei_score1_norm = (typei_score1 - min(typei_score1))/(max(typei_score1) - min(typei_score1));
    degree_cell{1,i} = types{i};
    eigen_cell{1,i} = types{i};
    degree_cell{2,i} = [typei_score1_norm, degree_C];                        
    eigen_cell{2,i} = [typei_score1_norm, eigen_C];    
end

types_plot = tasic_classes;
cmap = hsv(length(unique(tasic_ontology_cns(:,2))));
% nonneuronals = repmat({'Oligo','Endo','Macro','Astro'},1,2);
% nonneuronals_plot = {'Oligo','Endo','Macro','Astro'};
figure('Position',[0 0 1200 600]);
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(types_plot)
    nexttile; hold on; box on;
    datamat = eigen_cell{2,i};
    datamat([171,383],:) = []; % remove outlier region-pair
    scatter(datamat(:,1),datamat(:,2),40,'d','filled','MarkerEdgeColor',...
        cmap(i,:),'MarkerFaceColor',cmap(i,:));
    h = lsline;
    h.LineWidth = 3; h.Color = [0.75,0.75,0.25];
    ymax = 0.2; xmax = 1;
    ylim([0,ymax]); xlim([0,xmax]);
    yticks([0,ymax]); xticks([0,xmax]);
    [R,p] = corr(datamat(:,1),datamat(:,2),'type','Pearson');
    text(0.05*xmax,0.875*ymax,{sprintf("Pearson's R = %0.2f",R),sprintf("P value = %0.1d",p)},...
        'FontName','DejaVuSans','FontSize',16)
    title(types_plot{i})
    set(gca,'FontSize',18,'FontName','DejaVuSans');
end
xlabel(t,'Normalized cell-type density','FontSize',18,'FontName','DejaVuSans'); 
ylabel(t,'Eigenvector Centrality','FontSize',18,'FontName','DejaVuSans');
print([figdir filesep 'CentralityTasicAll'],'-dtiffn','-r300'); close;

%% Plot centrality measures for non-neuronals, Tasic
load([matdir filesep 'Centrality_Results.mat'],'centrality_results')
degree_cell = cell(2,25);
degree_res = centrality_results(ismember(centrality_results(:,2),'degree'),:);
eigen_cell = cell(2,25);
eigen_res = centrality_results(ismember(centrality_results(:,2),'eigenvector_centrality'),:);
celltypes = unique(degree_res(:,1));
for i = 1:length(celltypes)
    degree_cell{1,i} = celltypes{i};
    eigen_cell{1,i} = celltypes{i};
    ctinds = ismember(degree_res(:,1),celltypes{i});
    degree_cell{2,i} = cell2mat(degree_res(ctinds,3:4));
    eigen_cell{2,i} = cell2mat(eigen_res(ctinds,3:4));
end

nonneuronals = repmat({'Oligo','Endo','Macro','Astro'},1,2);
nonneuronals_plot = {'Oligo','Vasc','Immune','Astro'};
cmap = hsv(length(unique(tasic_ontology_cns(:,2))));
cmap = cmap([5,8,2,1],:);
cmap = repmat(cmap,2,1);
% nonneuronals = repmat({'Oligo','Endo','Macro','Astro'},1,2);
% nonneuronals_plot = {'Oligo','Endo','Macro','Astro'};
figure('Position',[0 0 1000 500]);
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
for i = 1:8
    nexttile; hold on; box on;
    if ismember(i,1:4)
        datamat = degree_cell{2,ismember(celltypes,nonneuronals{i})};
        datamat([171,383],:) = []; % remove outlier region-pair
        scatter(datamat(:,1),datamat(:,2),40,'o','filled','MarkerEdgeColor',...
            cmap(i,:),'MarkerFaceColor',cmap(i,:));
        h = lsline;
        h.LineWidth = 3; h.Color = [0.75,0.75,0.25];
        ymax = 2500; xmax = 1;
        ylim([0,ymax]); xlim([0,xmax]);
        yticks(ymax); xticks([0,xmax]);
        ytickformat('%.1f');
        ax = gca; ax.YAxis.Exponent = 3;
        if i == 1
            ylabel('Degree');
            text(-0.08,0,'0','FontSize',18,'FontName','DejaVuSans');
        else
            set(gca,'YTickLabel',[])
        end
        [rho,p] = corr(datamat(:,1),datamat(:,2));
        text(0.05*xmax,0.875*ymax,{sprintf("Pearson's R = %0.2f",rho),sprintf("P value = %0.1d",p)},...
            'FontName','DejaVuSans','FontSize',16)
        set(gca,'FontSize',18,'FontName','DejaVuSans');
        title(nonneuronals_plot{i},'FontSize',20,'FontName','DejaVuSans');
    else
        datamat = eigen_cell{2,ismember(celltypes,nonneuronals{i})};
        datamat([171,383],:) = []; % remove outlier region-pair
        scatter(datamat(:,1),datamat(:,2),40,'d','filled','MarkerEdgeColor',...
            cmap(i,:),'MarkerFaceColor',cmap(i,:));
        h = lsline;
        h.LineWidth = 3; h.Color = [0.75,0.75,0.25];
        ymax = 0.2; xmax = 1;
        ylim([0,ymax]); xlim([0,xmax]);
        yticks([0,ymax]); xticks([0,xmax]);
        if i == 5
            ylabel('Eigenvector Centrality');
        else
            set(gca,'YTickLabel',[])
        end
        [rho,p] = corr(datamat(:,1),datamat(:,2));
        text(0.05*xmax,0.875*ymax,{sprintf("Pearson's R = %0.2f",rho),sprintf("P value = %0.1d",p)},...
            'FontName','DejaVuSans','FontSize',16)
        set(gca,'FontSize',18,'FontName','DejaVuSans');
    end
end
xlabel(t,'Normalized cell-type density','FontSize',18,'FontName','DejaVuSans');
% print([figdir filesep 'NonneuronalCentrality_Tasic'],'-dtiffn','-r300'); close;

%% Plot eigenvector centrality for all types, Zeisel
load([matdir filesep 'Centrality_Results.mat'],'centrality_results');
load([matdir filesep 'Zeisel_outstruct_nG1360.mat'],'outstruct');
types = zeisel_classes;
degree_cell = cell(2,length(types));
degree_res = centrality_results(ismember(centrality_results(:,2),'degree'),:);
eigen_cell = cell(2,length(types));
eigen_res = centrality_results(ismember(centrality_results(:,2),'eigenvector_centrality'),:);
degree_C = cell2mat(degree_res(1:424,4)); eigen_C = cell2mat(eigen_res(1:424,4));

for i = 1:length(types)
    ctinds = ismember(zeisel_ontology_cns(:,2), types{i});
    celldensity_typei = outstruct(1).corrB(:,ctinds);   
    celldensity_typei_reg = zeros(426,sum(ctinds));
    for j = 1:sum(ctinds)
        [~,celldensity_typei_regj] = Voxel_To_Region(celldensity_typei(:,j),matdir);
        celldensity_typei_reg(:,j) = celldensity_typei_regj;
    end
    celldensity_typei_reg([12,(12+213)],:) = [];
    [~,scores_typei,~] = pca(celldensity_typei_reg);
    typei_score1 = scores_typei(:,1);   
    typei_score1_norm = (typei_score1 - min(typei_score1))/(max(typei_score1) - min(typei_score1));
    degree_cell{1,i} = types{i};
    eigen_cell{1,i} = types{i};
    degree_cell{2,i} = [typei_score1_norm, degree_C];                        
    eigen_cell{2,i} = [typei_score1_norm, eigen_C];    
end

types_plot = zeisel_classes;
cmap = hsv(length(unique(zeisel_ontology_cns(:,2))));
% nonneuronals = repmat({'Oligo','Endo','Macro','Astro'},1,2);
% nonneuronals_plot = {'Oligo','Endo','Macro','Astro'};
figure('Position',[0 0 1200 1000]);
t = tiledlayout(4,5,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(types_plot)
    nexttile; hold on; box on;
    datamat = eigen_cell{2,i};
    datamat([171,383],:) = []; % remove outlier region-pair
    scatter(datamat(:,1),datamat(:,2),40,'d','filled','MarkerEdgeColor',...
        cmap(i,:),'MarkerFaceColor',cmap(i,:));
    h = lsline;
    h.LineWidth = 3; h.Color = [0.75,0.75,0.25];
    ymax = 0.2; xmax = 1;
    ylim([0,ymax]); xlim([0,xmax]);
    yticks([0,ymax]); xticks([0,xmax]);
    [R,p] = corr(datamat(:,1),datamat(:,2),'type','Pearson');
    text(0.05*xmax,0.875*ymax,{sprintf("Pearson's R = %0.2f",R),sprintf("P value = %0.1d",p)},...
        'FontName','DejaVuSans','FontSize',16)
    title(types_plot{i})
    set(gca,'FontSize',18,'FontName','DejaVuSans');
end
xlabel(t,'Normalized cell-type density','FontSize',18,'FontName','DejaVuSans'); 
ylabel(t,'Eigenvector Centrality','FontSize',18,'FontName','DejaVuSans');
print([figdir filesep 'CentralityZeiselAll'],'-dtiffn','-r300'); close;

%% Plot centrality measures for non-neuronals, Zeisel
load([matdir filesep 'Centrality_Results.mat'],'centrality_results');
load([matdir filesep 'Zeisel_outstruct_nG1360.mat'],'outstruct');
types = {'Astro','Immune','Oligo','Vasc'};
degree_cell = cell(2,length(types));
degree_res = centrality_results(ismember(centrality_results(:,2),'degree'),:);
eigen_cell = cell(2,length(types));
eigen_res = centrality_results(ismember(centrality_results(:,2),'eigenvector_centrality'),:);
degree_C = cell2mat(degree_res(1:424,4)); eigen_C = cell2mat(eigen_res(1:424,4));

for i = 1:length(types)
    ctinds = ismember(zeisel_ontology_cns(:,2), types{i});
    celldensity_typei = outstruct(1).corrB(:,ctinds);   
    celldensity_typei_reg = zeros(426,sum(ctinds));
    for j = 1:sum(ctinds)
        [~,celldensity_typei_regj] = Voxel_To_Region(celldensity_typei(:,j),matdir);
        celldensity_typei_reg(:,j) = celldensity_typei_regj;
    end
    celldensity_typei_reg([12,(12+213)],:) = [];
    [~,scores_typei,~] = pca(celldensity_typei_reg);
    typei_score1 = scores_typei(:,1);   
    typei_score1_norm = (typei_score1 - min(typei_score1))/(max(typei_score1) - min(typei_score1));
    degree_cell{1,i} = types{i};
    eigen_cell{1,i} = types{i};
    degree_cell{2,i} = [typei_score1_norm, degree_C];                        
    eigen_cell{2,i} = [typei_score1_norm, eigen_C];    
end

nonneuronals = repmat({'Oligo','Astro','Vasc','Immune'},1,2);
nonneuronals_plot = {'Oligo','Astro','Vasc','Immune'};
cmap = hsv(length(unique(zeisel_ontology_cns(:,2))));
cmap = cmap([14,1,17,9],:);
cmap = repmat(cmap,2,1);
% nonneuronals = repmat({'Oligo','Endo','Macro','Astro'},1,2);
% nonneuronals_plot = {'Oligo','Endo','Macro','Astro'};
figure('Position',[0 0 1000 500]);
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
for i = 1:8
    nexttile; hold on; box on;
    if ismember(i,1:4)
        datamat = degree_cell{2,ismember(types,nonneuronals{i})};
        datamat([171,383],:) = []; % remove outlier region-pair
        scatter(datamat(:,1),datamat(:,2),40,'o','filled','MarkerEdgeColor',...
            cmap(i,:),'MarkerFaceColor',cmap(i,:));
        h = lsline;
        h.LineWidth = 3; h.Color = [0.75,0.75,0.25];
        ymax = 2500; xmax = 1;
        ylim([0,ymax]); xlim([0,xmax]);
        yticks(ymax); xticks([0,xmax]);
        ytickformat('%.1f');
        ax = gca; ax.YAxis.Exponent = 3;
        if i == 1
            ylabel('Degree');
            text(-0.08,0,'0','FontSize',18,'FontName','DejaVuSans');
        else
            set(gca,'YTickLabel',[])
        end
        [R,p] = corr(datamat(:,1),datamat(:,2),'type','Pearson');
        text(0.05*xmax,0.875*ymax,{sprintf("Pearson's R = %0.2f",R),sprintf("P value = %0.1d",p)},...
            'FontName','DejaVuSans','FontSize',16)
        set(gca,'FontSize',18,'FontName','DejaVuSans');
        title(nonneuronals_plot{i},'FontSize',20,'FontName','DejaVuSans');
    else
        datamat = eigen_cell{2,ismember(types,nonneuronals{i})};
        datamat([171,383],:) = []; % remove outlier region-pair
        scatter(datamat(:,1),datamat(:,2),40,'d','filled','MarkerEdgeColor',...
            cmap(i,:),'MarkerFaceColor',cmap(i,:));
        h = lsline;
        h.LineWidth = 3; h.Color = [0.75,0.75,0.25];
        ymax = 0.2; xmax = 1;
        ylim([0,ymax]); xlim([0,xmax]);
        yticks([0,ymax]); xticks([0,xmax]);
        if i == 5
            ylabel('Eigenvector Centrality');
        else
            set(gca,'YTickLabel',[])
        end
        [R,p] = corr(datamat(:,1),datamat(:,2),'type','Pearson');
        text(0.05*xmax,0.875*ymax,{sprintf("Pearson's R = %0.2f",R),sprintf("P value = %0.1d",p)},...
            'FontName','DejaVuSans','FontSize',16)
        set(gca,'FontSize',18,'FontName','DejaVuSans');
    end
end
xlabel(t,'Normalized cell-type density','FontSize',18,'FontName','DejaVuSans');
print([figdir filesep 'NonneuronalCentrality'],'-dtiffn','-r300'); close;

%% Brainframe images of non-neuronals, Tasic
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
addpath('/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe');
load([matdir filesep 'CellDensity_Tasic_nG606.mat'],'outstruct')
load([matdir filesep 'BrainFrame_Dependencies.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey','input_struct');
outstruct_tasic = outstruct;
input_struct.voxthresh = 0.65;
input_struct.nbin = 10;
input_struct.voxUreg = 0;
input_struct.xfac = 0.5; 
input_struct.pointsize = 0.1;
input_struct.savenclose = 1;
input_struct.bgcolor = 'w';
input_struct.img_format = 'tiffn';
input_struct.img_renderer = 'painters';
input_struct.img_directory = figdir;
img_views_ = [[0 0 1]; [-1 0 0]; [0 1 0]];
input_struct.img_views = img_views_;

% colors1 = [[1 0.5 0]; [1 1 1]];
% colors2 = [[1 1 0]; [1 0.5 0.5]];
colors_base = hsv(length(unique(tasic_ontology_cns(:,2))));
colors_base = colors_base([5,8,2,1],:);
colors1 = (colors_base + 2*ones(4,3))/3;
% colors1 = (2*colors_base + ones(4,3))/2;
colors2 = colors_base;
types = [24,25,23,22];

for i = 1:length(types)
    newVoxMap = zeros(size(GENGDmod));
    curinput = outstruct_tasic(1).corrB(:,types(i));
    for k = 1:length(structList)
        [~,voxinds] = ismember(structIndex{k},nonzerovox);
        allvox = nonzerovox(voxinds);
        curvox = curinput(voxinds);
        newVoxMap(allvox) = curvox;
    end
    datinput = newVoxMap;
    
    %stuff to create data goes here
    datinput = imresize3(datinput,[133 81 115]);
    input_struct.data = datinput;
    input_struct.regsUbins = 0;
    input_struct.img_labels = [classkey{types(i)}];
%     input_struct.voxthresh = voxthreshes(i);
    input_struct.cmap = twocolor(colors1(i,:),colors2(i,:),input_struct.nbin);
    brainframe(input_struct);   
end

%% Brainframe images of non-neuronals, Zeisel
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
addpath('/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe');
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'BrainFrame_Dependencies.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey','input_struct');
load([matdir filesep 'Zeisel_outstruct_nG1360.mat'],'outstruct');
load([matdir filesep 'Zeisel_Ontology.mat'],'zeisel_ontology_cns');

vts = [0.65,0.75,0.65,0.65];
input_struct.nbin = 10;
input_struct.voxUreg = 0;
input_struct.xfac = 0.5; 
input_struct.pointsize = 0.1;
input_struct.savenclose = 1;
input_struct.bgcolor = 'w';
input_struct.img_format = 'tiffn';
input_struct.img_renderer = 'painters';
input_struct.img_directory = figdir;
img_views_ = [0 0 1];
input_struct.img_views = img_views_;

% colors1 = [[1 0.5 0]; [1 1 1]];
% colors2 = [[1 1 0]; [1 0.5 0.5]];
cmap = hsv(length(unique(zeisel_ontology_cns(:,2))));
colors_base = cmap([14,1,17,9],:);
colors1 = (colors_base + 2*ones(4,3))/3;
% colors1 = (2*colors_base + ones(4,3))/2;
colors2 = colors_base;
types = {'Oligo','Astro','Vasc','Immune'};

for i = 1:length(types)
    newVoxMap = zeros(size(GENGDmod));
    ctinds = ismember(zeisel_ontology_cns(:,2), types{i});
    celldensity_typei = outstruct(1).corrB(:,ctinds); 
    [~,scores_typei,~] = pca(celldensity_typei);
    typei_score1 = scores_typei(:,1);   
    typei_score1_norm = (typei_score1 - min(typei_score1))/(max(typei_score1) - min(typei_score1));
    curinput = typei_score1_norm;
    for k = 1:length(structList)
        [~,voxinds] = ismember(structIndex{k},nonzerovox);
        allvox = nonzerovox(voxinds);
        curvox = curinput(voxinds);
        newVoxMap(allvox) = curvox;
    end
    datinput = newVoxMap;
    
    %stuff to create data goes here
    input_struct.voxthresh = vts(i);
    datinput = imresize3(datinput,[133 81 115]);
    input_struct.data = datinput;
    input_struct.regsUbins = 0;
    input_struct.img_labels = ['Zeisel_' types{i}];
%     input_struct.voxthresh = voxthreshes(i);
    input_struct.cmap = twocolor(colors1(i,:),colors2(i,:),input_struct.nbin);
    brainframe(input_struct);   
end


%% Brainframe images of connectome  
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
matdir1 = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'BrainFrame_Dependencies.mat'],'input_struct');
load([matdir filesep 'Interregional_Distances_JLT.mat']);

reggroups = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57; 
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
reggroups(tha) = 13;
reggroups = [reggroups;reggroups];
cmap = hsv(length(unique(reggroups))); %Creating colormap

D_ = cat(1,D(1:11,:),zeros(1,size(D,2)),D(12:(11+212),:),zeros(1,size(D,2)),D((12+212):end,:));
D_ = cat(2,D_(:,1:11),zeros(size(D_,1),1),D_(:,12:(11+212)),zeros(size(D_,1),1),D_(:,(12+212):end));
D = D_;

% Hippocampal connectivity, globally thresholded
% distthresh = 75; conthresh = 85;
% C = input_struct.conmat;
% C = C - diag(diag(C));
% conthresh_val = prctile(nonzeros(C(:)),conthresh);
% distthresh_val = prctile(nonzeros(D(:)),distthresh);
% thresh_inds = ((D > distthresh_val) + (C > conthresh_val)) == 2;
% C_thresh = C;
% C_thresh(~thresh_inds) = 0;
% 
% hipp_inds = [hip,(hip)+213];
% datinput = zeros(size(C,1),1);
% C_hipp_thresh = C_thresh; 
% C_hipp_thresh(:,~ismember(1:size(C,1),hipp_inds)) = 0;
% [x,y] = find(C_hipp_thresh > 0);
% datinput(y) = 1; datinput(unique(x)) = 1;
% imglabel = sprintf('HippocampalConnectivity_D%d_C%d',distthresh,conthresh);
% input_struct_hipp = brainframe_inputs_mouse(matdir1,'conmat',C_hipp_thresh,...
%                                              'region_groups',reggroups,...
%                                              'con_regiongroups',reggroups,...
%                                              'cmap',cmap,...
%                                              'con_cmap',cmap,...
%                                              'xfac',2,...
%                                              'sphere',1,...
%                                              'voxUreg',1,...
%                                              'iscon',1,...
%                                              'conarrow_WL',[1 0.75],...
%                                              'data',datinput,...
%                                              'norm_method','max',...
%                                              'bgcolor','w',...
%                                              'img_labels',imglabel);
% brainframe(input_struct_hipp);

% Neocortical incoming connectivity, globally thresholded
distthresh_high = 75; conthresh = 95;
C = input_struct.conmat;
C = C - diag(diag(C));
conthresh_val = prctile(nonzeros(C(:)),conthresh);
distthresh_val_low = prctile(nonzeros(D(:)),distthresh_high);
thresh_inds = ((D > distthresh_val_low) + (C > conthresh_val)) == 2;
C_thresh = C;
C_thresh(~thresh_inds) = 0;

neo_inds = [ncx,(ncx)+213];
datinput = zeros(size(C,1),1);
C_neo_thresh = C_thresh; 
C_neo_thresh(:,~ismember(1:size(C,1),neo_inds)) = 0;
[x,y] = find(C_neo_thresh > 0);
datinput(y) = 1; datinput(unique(x)) = 1;
imglabel = sprintf('NeocorticalInConnectivity_D%d_C%0.1f',distthresh_high,conthresh);
imglabel = strrep(imglabel,'.','_');
imgview = [-1.068505859375000e+02,-32.923828125000000];
input_struct_neo = brainframe_inputs_mouse(matdir1,'conmat',C_neo_thresh,...
                                             'region_groups',reggroups,...
                                             'con_regiongroups',reggroups,...
                                             'cmap',cmap,...
                                             'con_cmap',cmap,...
                                             'xfac',2,...
                                             'sphere',1,...
                                             'sphere_npts',15,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'iscon',1,...
                                             'conarrow_WL',[1.5 1],...
                                             'data',datinput,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'con_rescale',20,...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir,...
                                             'savenclose',1,...
                                             'con_arch',0.3);
brainframe(input_struct_neo);

% Neocortical outgoing connectivity, globally thresholded
distthresh_high = 75; conthresh = 95;
C = input_struct.conmat;
C = C - diag(diag(C));
conthresh_val = prctile(nonzeros(C(:)),conthresh);
distthresh_val_low = prctile(nonzeros(D(:)),distthresh_high);
thresh_inds = ((D > distthresh_val_low) + (C > conthresh_val)) == 2;
C_thresh = C;
C_thresh(~thresh_inds) = 0;

neo_inds = [ncx,(ncx)+213];
datinput = zeros(size(C,1),1);
C_neo_thresh = C_thresh; 
C_neo_thresh(~ismember(1:size(C,1),neo_inds),:) = 0;
[x,y] = find(C_neo_thresh > 0);
datinput(y) = 1; datinput(unique(x)) = 1;
imglabel = sprintf('NeocorticalOutConnectivity_D%d_C%0.1f',distthresh_high,conthresh);
imglabel = strrep(imglabel,'.','_');
imgview = [-1.068505859375000e+02,-32.923828125000000];
input_struct_neo = brainframe_inputs_mouse(matdir1,'conmat',C_neo_thresh,...
                                             'region_groups',reggroups,...
                                             'con_regiongroups',reggroups,...
                                             'cmap',cmap,...
                                             'con_cmap',cmap,...
                                             'xfac',2,...
                                             'sphere',1,...
                                             'sphere_npts',15,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'iscon',1,...
                                             'conarrow_WL',[1.5 1],...
                                             'data',datinput,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'con_rescale',20,...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir,...
                                             'savenclose',1,...
                                             'con_arch',0.3);
brainframe(input_struct_neo);

% Striatal outgoing connectivity, globally thresholded
distthresh_high = 75; conthresh = 50;
C = input_struct.conmat;
C = C - diag(diag(C));
conthresh_val = prctile(nonzeros(C(:)),conthresh);
distthresh_val_low = prctile(nonzeros(D(:)),distthresh_high);
thresh_inds = ((D > distthresh_val_low) + (C > conthresh_val)) == 2;
C_thresh = C;
C_thresh(~thresh_inds) = 0;

str_inds = [str,(str)+213];
datinput = zeros(size(C,1),1);
C_str_thresh = C_thresh; 
C_str_thresh(~ismember(1:size(C,1),str_inds),:) = 0;
[x,y] = find(C_str_thresh > 0);
datinput(y) = 1; datinput(unique(x)) = 1;
imglabel = sprintf('StriatumInConnectivity_D%d_C%d',distthresh_high,conthresh);
imgview = [-1.068505859375000e+02,-32.923828125000000];
input_struct_str = brainframe_inputs_mouse(matdir1,'conmat',C_str_thresh,...
                                             'region_groups',reggroups,...
                                             'con_regiongroups',reggroups,...
                                             'cmap',cmap,...
                                             'con_cmap',cmap,...
                                             'xfac',2,...
                                             'sphere',1,...
                                             'sphere_npts',15,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'iscon',1,...
                                             'conarrow_WL',[1.5 1],...
                                             'data',datinput,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'con_rescale',20,...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir,...
                                             'savenclose',1,...
                                             'con_arch',0.3);
brainframe(input_struct_str);

% Striatal incoming connectivity, globally thresholded
distthresh_high = 75; conthresh = 50;
C = input_struct.conmat;
C = C - diag(diag(C));
conthresh_val = prctile(nonzeros(C(:)),conthresh);
distthresh_val_low = prctile(nonzeros(D(:)),distthresh_high);
thresh_inds = ((D > distthresh_val_low) + (C > conthresh_val)) == 2;
C_thresh = C;
C_thresh(~thresh_inds) = 0;

str_inds = [str,(str)+213];
datinput = zeros(size(C,1),1);
C_str_thresh = C_thresh; 
C_str_thresh(:,~ismember(1:size(C,1),str_inds)) = 0;
[x,y] = find(C_str_thresh > 0);
datinput(y) = 1; datinput(unique(x)) = 1;
imglabel = sprintf('StriatumOutConnectivity_D%d_C%d',distthresh_high,conthresh);
imgview = [-1.068505859375000e+02,-32.923828125000000];
input_struct_str = brainframe_inputs_mouse(matdir1,'conmat',C_str_thresh,...
                                             'region_groups',reggroups,...
                                             'con_regiongroups',reggroups,...
                                             'cmap',cmap,...
                                             'con_cmap',cmap,...
                                             'xfac',2,...
                                             'sphere',1,...
                                             'sphere_npts',15,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'iscon',1,...
                                             'conarrow_WL',[1.5 1],...
                                             'data',datinput,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'con_rescale',20,...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figdir,...
                                             'savenclose',1,...
                                             'con_arch',0.3);
brainframe(input_struct_str);


%% Tasic class brainframe
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
matdir1 = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'Tasic_Inputs.mat'],'nonzerovox','GENGDmod');
load([matdir filesep 'CellDensity_Tasic_nG606.mat'],'outstruct');

reggroups = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57; 
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
reggroups(tha) = 13;
reggroups = [reggroups;reggroups];
cmap = hsv(length(unique(reggroups))); %Creating colormap

celldensity_voxel = outstruct.corrB;
cdv_min = min(celldensity_voxel); cdv_min = repmat(cdv_min,size(celldensity_voxel,1),1);
cdv_max = max(celldensity_voxel); cdv_max = repmat(cdv_max,size(celldensity_voxel,1),1);
celldensity_norm = (celldensity_voxel - cdv_min) ./ (cdv_max - cdv_min);

% NEOGLUT
glutneo_inds = 9:16;
celldensity_glutneo = celldensity_norm(:,glutneo_inds);
[~,scores_glutneo,~] = pca(celldensity_glutneo);
% glutneo_score1 = scores_glutneo(:,1);
glutneo_score1 = mean(celldensity_glutneo,2);
glutneo_score1_norm = (glutneo_score1 - min(glutneo_score1))/(max(glutneo_score1) - min(glutneo_score1));

newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = glutneo_score1_norm;
datinput = newVoxMap;
datinput = imresize3(datinput,[133 81 115]);
datinput(datinput < 0) = 0;

% nbin = 10;
% color1 = [1 1 0];
% color2 = [1 0.5 0];
img_views_ = [[0 0 1]; [-1 0 0]; [0 1 0]];
% cmap = twocolor(color1,color2,nbin);

input_struct_glutneo = brainframe_inputs_mouse(matdir1,'data',datinput,...
                                            'voxthresh',0.65,...
                                            'nbin',10,...
                                            'voxUreg',0,...
                                            'xfac',0.02,...
                                            'pointsize',0.1,...
                                            'bgcolor','w',...
                                            'img_format','tiffn',...
                                            'cmap',cmap,...
                                            'regsUbins',1,...
                                            'region_groups',reggroups,...
                                            'img_directory',figdir,...
                                            'img_labels','NEOGLUT_PC1',...
                                            'img_views',img_views_,...
                                            'img_renderer','painters',...
                                            'savenclose',0);
brainframe(input_struct_glutneo); 

%% Zeisel class brainframe
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
matdir1 = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'BrainFrame_Dependencies.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey','input_struct');
load([matdir filesep 'Zeisel_outstruct_nG1360.mat'],'outstruct');

reggroups = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57; 
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
reggroups(tha) = 13;
reggroups = [reggroups;reggroups];
cmap = hsv(length(unique(reggroups))); %Creating colormap

celldensity_voxel = outstruct.corrB;
cdv_min = min(celldensity_voxel); cdv_min = repmat(cdv_min,size(celldensity_voxel,1),1);
cdv_max = max(celldensity_voxel); cdv_max = repmat(cdv_max,size(celldensity_voxel,1),1);
celldensity_norm = (celldensity_voxel - cdv_min) ./ (cdv_max - cdv_min);

% TEGLU
teglu_comp = @(x) (length(x) > 5) && strcmp(x(1:5),'TEGLU');
teglu_inds = cell2mat(cellfun(teglu_comp,classkey,'UniformOutput',false));
celldensity_teglu = celldensity_norm(:,teglu_inds);
[~,scores_teglu,~] = pca(celldensity_teglu);
teglu_score1 = scores_teglu(:,1);
teglu_score1_norm = (teglu_score1 - min(teglu_score1))/(max(teglu_score1) - min(teglu_score1));

newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = teglu_score1_norm;
datinput = newVoxMap;
datinput = imresize3(datinput,[133 81 115]);
datinput(datinput < 0) = 0;

% nbin = 10;
% color1 = [1 1 0];
% color2 = [1 0.5 0];
img_views_ = [[0 0 1]; [-1 0 0]; [0 1 0]];
% cmap = twocolor(color1,color2,nbin);

input_struct_teglu = brainframe_inputs_mouse(matdir1,'data',datinput,...
                                            'voxthresh',0.65,...
                                            'nbin',10,...
                                            'voxUreg',0,...
                                            'xfac',0.02,...
                                            'pointsize',0.1,...
                                            'bgcolor','w',...
                                            'img_format','tiffn',...
                                            'cmap',cmap,...
                                            'regsUbins',1,...
                                            'region_groups',reggroups,...
                                            'img_directory',figdir,...
                                            'img_labels','TEGLU_PC1',...
                                            'img_views',img_views_,...
                                            'img_renderer','painters',...
                                            'savenclose',0);
brainframe(input_struct_teglu); 

% MSN
msn_comp = @(x) (length(x) > 3) && strcmp(x(1:3),'MSN');
msn_inds = cell2mat(cellfun(msn_comp,classkey,'UniformOutput',false));
celldensity_msn = celldensity_norm(:,msn_inds);
[~,scores_msn,~] = pca(celldensity_msn);
msn_score1 = scores_msn(:,1);
msn_score1_norm = (msn_score1 - min(msn_score1))/(max(msn_score1) - min(msn_score1));

newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = msn_score1_norm;
datinput = newVoxMap;
datinput = imresize3(datinput,[133 81 115]);
datinput(datinput < 0) = 0;

% nbin = 10;
% color1 = [1 1 0];
% color2 = [1 0.5 0];
img_views_ = [[0 0 1]; [-1 0 0]; [0 1 0]];
% cmap = twocolor(color1,color2,nbin);

input_struct_msn = brainframe_inputs_mouse(matdir1,'data',datinput,...
                                            'voxthresh',0.65,...
                                            'nbin',10,...
                                            'voxUreg',0,...
                                            'xfac',0.02,...
                                            'pointsize',0.1,...
                                            'bgcolor','w',...
                                            'img_format','tiffn',...
                                            'cmap',cmap,...
                                            'regsUbins',1,...
                                            'region_groups',reggroups,...
                                            'img_directory',figdir,...
                                            'img_labels','MSN_PC1',...
                                            'img_views',img_views_,...
                                            'img_renderer','painters',...
                                            'savenclose',0);
brainframe(input_struct_msn); 

%% Long-range proportion analysis
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'BrainFrame_Dependencies.mat'],'input_struct');
load([matdir filesep 'Interregional_Distances_JLT.mat']);
load([matdir filesep 'regionvoxels.mat']);

amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57; 
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
region_cell = {amy,cer,sub,hip,hyp,ncx,med,mid,olf,pal,pon,str,tha};
region_name_cell = {'Amygdala','Cerebellum','Cortical Subplate','Hippocampus',...
    'Hypothalamus','Neocortex','Medulla','Midbrain','Olfactory','Pallidum',...
    'Pons','Striatum','Thalamus'};

D_ = cat(1,D(1:11,:),zeros(1,size(D,2)),D(12:(11+212),:),zeros(1,size(D,2)),D((12+212):end,:));
D_ = cat(2,D_(:,1:11),zeros(size(D_,1),1),D_(:,12:(11+212)),zeros(size(D_,1),1),D_(:,(12+212):end));
D = D_;

V = [voxels; voxels];
V_row = repmat(V,1,length(V));
% V_col = repmat(V.',length(V),1);
% V_arimean = 0.5*(V_row + V_col);
% V_geomean = (V_row .* V_col).^(0.5);

distthresh_high = 75;
C = input_struct.conmat;
C = C - diag(diag(C));
C = log2(C+1);
distthresh_val_low = prctile(nonzeros(D(:)),distthresh_high);
thresh_inds = (D < distthresh_val_low);
C(~thresh_inds) = 0;
C_row = C ./ V_row;
% C_col = C ./ V_col;
% C_arimean = C ./ V_arimean;
% C_geomean = C ./ V_geomean;
C_ratios_reg = zeros(1,length(region_cell));
V_ratios_reg = C_ratios_reg;
for i = 1:length(region_cell)
    reginds = zeros(1,length(V));
    reginds([region_cell{i},(region_cell{i}+213)]) = 1;
    reginds = logical(reginds);
    C_row_reg = C_row; 
    C_row_reg(~reginds,~reginds) = 0;
    C_ratios_reg(i) = sum(C_row_reg(:))/sum(C_row(:));
    V_ratios_reg(i) = sum(V(reginds))/sum(V);
end
[C_ratios_reg_sort,sortinds] = sort(C_ratios_reg,'descend');
region_name_sort = region_name_cell(sortinds);

%% MISS rerunning on dummy genes, Tasic
rng(0);
matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles'; %define directory to draw from and save data to
load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
load([matdir filesep 'MRx3_L90_inds'],'geneinds'); %Tasic MRx3 gene indices
lambda = 90;
ng_param_list = 606; %values of nG to test and map, going through genes in MRx3 ranked order
missmethod = 'MR                                                                                                                    x3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
% savename = 'CellDensity_Tasic_nG606_random.mat'; %set name of file to be saved
niters = 100;
scrambling = 1;
Tasic_scrambleddensities_nG606 = zeros(424,length(classkey),niters);
[~,regvgene] = Voxel_To_Region(voxvgene,matdir);
clc;
for i = 1:niters
    fprintf('Iteration %d/%d\n',i,niters);
    tic
    [E_red,C_red,nGen] = GeneSelector(genevct,regvgene,gene_names,1360,90,missmethod,geneinds,scrambling);
    B = CellDensityInference_resnorm(E_red,C_red);
    B([12,12+213],:) = [];
    Tasic_scrambleddensities_nG606(:,:,i) = B;
    toc
end
save([matdir filesep 'Tasic_scrambledC_nG606_new'],'Tasic_scrambleddensities_nG606','-v7.3');

%% MISS rerunning on dummy genes, Zeisel
rng(0);
matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles'; %define directory to draw from and save data to
load([matdir filesep 'Zeisel_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
load([matdir filesep 'Zeisel_MRx3Inds.mat'],'geneinds'); %Zeisel MRx3 gene indice
missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
savename = 'CellDensity_Zeisel_nG1360_random.mat'; %set name of file to be saved
niters = 100;
scrambling = 1;
Zeisel_scrambleddensities_nG1360 = zeros(424,length(classkey),niters);      
[~,regvgene] = Voxel_To_Region(voxvgene,matdir);
clc;
for i = 1:niters
    fprintf('Iteration %d/%d\n',i,niters);
    tic
    [E_red,C_red,nGen] = GeneSelector(genevct,regvgene,gene_names,1360,90,missmethod,geneinds,scrambling);
    B = CellDensityInference_resnorm(E_red,C_red);
    B([12,12+213],:) = [];
    Zeisel_scrambleddensities_nG1360(:,:,i) = B;
    toc
end
save([matdir filesep 'Zeisel_scrambledC_nG1360_new.mat'],'Zeisel_scrambleddensities_nG1360','-v7.3');

%% Distance vs. connectivity density
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'BrainFrame_Dependencies.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey','input_struct');
load([matdir filesep 'Interregional_Distances_JLT.mat']);
load([matdir filesep 'regionvoxels.mat']);
fine = 0;
amy = 1:11; cer = 12:22; sub = 23:25; hip = 26:36; hyp = 37:56; 
ncx = 57:94; med = 95:119; mid = 120:140; olf = 141:148; pal = 149:156;
pon = 157:169; str = 170:177; tha = 178:212;
if fine
    region_cell = {amy,cer,sub,hip,hyp,ncx,med,mid,olf,pal,pon,str,tha};
    region_name_cell = {'Amygdala','Cerebellum','Cortical Subplate','Hippocampus',...
        'Hypothalamus','Neocortex','Medulla','Midbrain','Olfactory','Pallidum',...
        'Pons','Striatum','Thalamus'};
else
    region_name_cell = {'Neocortex','Other Forebrain','Di/Mesencephelon','Hindbrain'};
    region_cell = {ncx,[amy,sub,hip,olf,pal,str],[hyp,mid,tha],[cer,med,pon]};
end
cmap = hsv(length(region_cell));

% D_ = cat(1,D(1:11,:),zeros(1,size(D,2)),D(12:(11+212),:),zeros(1,size(D,2)),D((12+212):end,:));
% D_ = cat(2,D_(:,1:11),zeros(size(D_,1),1),D_(:,12:(11+212)),zeros(size(D_,1),1),D_(:,(12+212):end));
% D = D_;

V = [voxels; voxels];
V_row = repmat(V,1,length(V));
% V_col = repmat(V.',length(V),1);
% V_arimean = 0.5*(V_row + V_col);
% V_geomean = (V_row .* V_col).^(0.5);

distthresh_high = 75;
distthresh_low = 25;
C = input_struct.conmat;
C = C - diag(diag(C));
C = log2(C+1);
C_row = C ./ V_row;
C([12,213+12],:) = []; C(:,[12,213+12]) = [];
C_nonzeroinds = find(C);
D_nonzeros = D(C_nonzeroinds);
distthresh_val_high = prctile(nonzeros(D_nonzeros(:)),distthresh_high);
distthresh_val_low = prctile(nonzeros(D_nonzeros(:)),distthresh_low);
% C_col = C ./ V_col;
% C_arimean = C ./ V_arimean;
% C_geomean = C ./ V_geomean;
region_D_cell = cell(length(region_cell),2);
for i = 1:length(region_cell)
    C1 = zeros(size(C));
    C2 = C1;
    C3 = C1;
    reginds = region_cell{i};
    C1(reginds,:) = C(reginds,:);
    C2(:,reginds) = C(:,reginds);
    C3(reginds,reginds) = C(reginds,reginds);
    C_reg = C1 + C2 - C3; % remove double-counting in center block-diagonal
    C_reginds = find(C_reg);
    D_reg = D(C_reginds); 
    [f,xi] = ksdensity(D_reg(:));                               
    region_D_cell{i,1} = xi; 
    region_D_cell{i,2} = f*length(D_reg(:));
end 
figure('Units','inches','Position',[0 0 15 8.5]); hold on;
xmaxes = zeros(1,length(region_cell)); ymaxes = xmaxes;
for i = 1:length(region_cell)
    plot(region_D_cell{i,1},region_D_cell{i,2},'Color',cmap(i,:),'LineWidth',2);
    xmaxes(i) = max(region_D_cell{i,1}); ymaxes(i) = max(region_D_cell{i,2});
end
plot([distthresh_val_low,distthresh_val_low],[0,max(ymaxes)],'k--','LineWidth',1.5);
plot([distthresh_val_high,distthresh_val_high],[0,max(ymaxes)],'k:','LineWidth',1.5);
ylim([0,max(ymaxes)]); yticks([0,max(ymaxes)/2,max(ymaxes)]); 
yticklabels({'0',num2str(max(ymaxes)/2,'%1.f'),num2str(max(ymaxes),'%1.f')});
xlim([0,max(xmaxes)]); xticks([0,max(xmaxes)/2,max(xmaxes)]); 
xticklabels({'0',num2str(max(xmaxes)/2,'%.1f'),num2str(max(xmaxes),'%.1f')});
legend([region_name_cell, sprintf('Lower Quartile - %.1f mm',distthresh_val_low),...
    sprintf('Upper Quartile - %.1f mm',distthresh_val_high)],'FontName','DejaVuSans','FontSize',20,'box','off');
xlabel('Distance (mm)'); ylabel('Probability Density'); set(gca,'FontName','DejaVuSans','FontSize',24,'box','on');
print([figdir filesep 'ConnectionDistanceKSPlot'],'-dtiffn','-r300'); close;

%% Taxonomic distance hea   tmap
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'Taxonomic_Distance_Info.mat'],'taxonomic_distance_matrix')
cmap = flipud(copper(500));
figure;
imagesc(taxonomic_distance_matrix); colormap(cmap) ; colorbar;
set(gca,'FontSize',16,'FontName','DejaVuSans','TickLength',[0 0])
print([figdir filesep 'TaxonomicDistance'],'-dtiffn','-r300'); close;

%% Entropy calculation
% Zeisel
matdir = '/Users/justintorok/Documents/MATLAB/CellTypeConnectivity/MatFiles';
figdir = '/Users/justintorok/Documents/MATLAB/MISS/Figures';
load([matdir filesep 'Zeisel_outstruct_nG1360.mat'],'outstruct');
load([matdir filesep 'Zeisel_Ontology.mat'],'zeisel_ontology_cns');
load([matdir filesep 'Zeisel_FeatureImportance.mat']);
% Average cell types over broader classes in Zeisel ontology
zeisel_classes = unique(zeisel_ontology_cns(:,2));
zeisel_fi_all_class = zeros(length(zeisel_classes),size(zeisel_fi_all,2));
zeisel_fi_long_class = zeisel_fi_all_class;
zeisel_fi_short_class = zeisel_fi_all_class;
zeisel_fi_out_class = zeisel_fi_all_class;
zeisel_fi_in_class = zeisel_fi_all_class;
meth = 'CoV';

for i = 1:length(zeisel_classes)
    classinds = ismember(zeisel_ontology_cns(:,2),zeisel_classes{i});
    zeisel_fi_all_class(i,:) = mean(zeisel_fi_all(classinds,:));
    zeisel_fi_short_class(i,:) = mean(zeisel_fi_short(classinds,:));
    zeisel_fi_long_class(i,:) = mean(zeisel_fi_long(classinds,:));
    zeisel_fi_out_class(i,:) = mean(zeisel_fi_out(classinds,:));
    zeisel_fi_in_class(i,:) = mean(zeisel_fi_in(classinds,:));
end

pcascores = zeros(424,length(zeisel_classes));   
for i = 1:length(zeisel_classes)
    ctinds = ismember(zeisel_ontology_cns(:,2),zeisel_classes{i});
    celldensity_typei_reg = outstruct(1).Bmeans(:,ctinds);
    celldensity_typei_reg([12,(12+213)],:) = [];
    [~,scores_typei,~] = pca(celldensity_typei_reg);
    typei_score1 = scores_typei(:,1);   
    typei_score1_norm = (typei_score1 - min(typei_score1))/(max(typei_score1) - min(typei_score1));
    pcascores(:,i) = typei_score1_norm;
end
% pcascores = (pcascores - repmat(min(pcascores,[],2),1,size(pcascores,2)))./...
%     (repmat(max(pcascores,[],2),1,size(pcascores,2))-repmat(min(pcascores,[],2),1,size(pcascores,2)));
H_zeisel = dispersion_calc(pcascores,meth);                      
clear outstruct

% Tasic
load([matdir filesep 'CellDensity_Tasic_nG606.mat'],'outstruct');
load([matdir filesep 'Tasic_Ontology.mat'],'tasic_ontology_cns');
load([matdir filesep 'Tasic_FeatureImportance.mat']);
% Average cell types over broader classes in Tasic ontology
tasic_classes = unique(tasic_ontology_cns(:,2));
tasic_fi_all_class = zeros(length(tasic_classes),size(tasic_fi_all,2));
tasic_fi_long_class = tasic_fi_all_class;
tasic_fi_short_class = tasic_fi_all_class;
tasic_fi_classification_class = tasic_fi_all_class;
tasic_fi_out_class = tasic_fi_all_class;
tasic_fi_in_class = tasic_fi_all_class;
for i = 1:length(tasic_classes)
    classinds = ismember(tasic_ontology_cns(:,2),tasic_classes{i});
    tasic_fi_all_class(i,:) = mean(tasic_fi_all(classinds,:));
    tasic_fi_short_class(i,:) = mean(tasic_fi_short(classinds,:));
    tasic_fi_long_class(i,:) = mean(tasic_fi_long(classinds,:));
    tasic_fi_classification_class(i,:) = mean(tasic_fi_classification(classinds,:));
    tasic_fi_out_class(i,:) = mean(tasic_fi_out(classinds,:));
    tasic_fi_in_class(i,:) = mean(tasic_fi_in(classinds,:));
end

pcascores = zeros(424,length(tasic_classes));
for i = 1:length(tasic_classes) 
    ctinds = ismember(tasic_ontology_cns(:,2), tasic_classes{i});
    celldensity_typei_reg = outstruct(1).Bmeans(:,ctinds);
    celldensity_typei_reg([12,(12+213)],:) = [];
    if sum(ctinds) > 1
        [~,scores_typei,~] = pca(celldensity_typei_reg);
        typei_score1 = scores_typei(:,1);   
    else
        typei_score1 = celldensity_typei_reg;
    end
    typei_score1_norm = (typei_score1 - min(typei_score1))/(max(typei_score1) - min(typei_score1));
    pcascores(:,i) = typei_score1_norm;
end
% pcascores = (pcascores - repmat(min(pcascores,[],2),1,size(pcascores,2)))./...
%     (repmat(max(pcascores,[],2),1,size(pcascores,2))-repmat(min(pcascores,[],2),1,size(pcascores,2)));
H_tasic = dispersion_calc(pcascores,meth);          
fi_cell = [{tasic_fi_all_class,tasic_fi_short_class,tasic_fi_long_class,tasic_fi_out_class,tasic_fi_in_class};...
    {zeisel_fi_all_class,zeisel_fi_short_class, zeisel_fi_long_class,zeisel_fi_out_class,zeisel_fi_in_class}];
H_cell = {H_tasic,H_zeisel};
labelcell = {'Tasic FI','Zeisel FI'};
ficellnames = {'All Connections','Short-range','Long-range','Source','Target'};
% titlecell = ficellnames;
titlecell = {'All Connections','Short-range','Long-range'};

% Plot
cmap = [[1 0 0]; [0 1 0]; [0 0 1]];
cmap_z = hsv(17);
figure('units','inches','Position',[0 0 15 9]); 
% figure('units','inches','Posit    ion',[0 0 20 10]); 
tiledlayout(length(labelcell),length(titlecell),'TileSpacing','Compact','Padding','Compact');
fiinds = find(ismember(ficellnames,titlecell));
for i = 1:length(labelcell)
    for j = fiinds
        nexttile; hold on;
        fi_mean = mean(fi_cell{i,j},2); fi_std = std(fi_cell{i,j},[],2);
        if i == 1
             errorbar(H_cell{i},fi_mean,fi_std,'o','LineStyle','none',...
            'MarkerEdgeColor',cmap(j,:),'MarkerFaceColor',cmap(j,:),'MarkerSize',7,...
            'LineWidth',2,'Color',cmap(j,:));
            title(ficellnames{j});
            yticks([1.1*min(fi_mean),max(fi_mean)]);
            ytickformat('%.2f');
        end
        if i == 2
            for k = 1:17                 
                errorbar(H_cell{i}(k),fi_mean(k),fi_std(k),'o','LineStyle','none',...
                'MarkerEdgeColor',cmap_z(k,:),'MarkerFaceColor',cmap_z(k,:),'MarkerSize',7,...
                'LineWidth',2,'Color',cmap_z(k,:));
            end
            mdl = fitlm(H_cell{i},fi_mean);
            xi = linspace(min(H_cell{i}),max(H_cell{i})).';
            yi = predict(mdl,xi);
            plot(xi,yi,'k-','LineWidth',1.5,'Color',[0.75,0.75,0.25]);
            [R,p] = corr(H_cell{i},fi_mean);
            if j == 3
                text(0.025,0.09,sprintf("Pearson\'s R = %.2f \nP value = %.1d",...
                    R,p),'Units','normalized','FontSize',18,'FontName','DejaVuSans');
            elseif j == 1
                text(0.025,0.2,sprintf("Pearson\'s R = %.2f \nP value = %.1d",...
                    R,p),'Units','normalized','FontSize',18,'FontName','DejaVuSans');
            elseif j == 2
                text(0.025,0.2,sprintf("Pearson\'s R = %.2f \nP value = %.1d",...
                    R,p),'Units','normalized','FontSize',18,'FontName','DejaVuSans');
                xlabel(sprintf('%s',meth));
            end
%             text(0.1,0.1,sprintf("Pearson\'s R = %.2f \np = %.2f",...
%                 mdl.Rsquared.Adjusted,mdl.ModelFitVsNullModel.Pvalue),...
%                 'Units','normalized','FontSize',14,'FontName','DejaVuSans');
            yticks([1.1*min(fi_mean),max(fi_mean)]);
            ylim([min(fi_mean-fi_std),max(fi_mean+fi_std)])
            ytickformat('%.1f');
        end
        if j == 1
            ylabel('F.I.');
        end
        xticks([1.1*min(H_cell{i}),max(H_cell{i})]);
        xtickformat('%.1f');
        set(gca,'FontSize',20,'FontName','DejaVuSans','box','on')
    end
end
print([figdir filesep 'DispersionvsFIplot'],'-dtiffn','-r300'); close;


function H = dispersion_calc(X,method)
    if nargin < 2
        method = 'Entropy';
    end
    if strcmp(method,'Entropy')
        X = X ./ repmat(sum(X),size(X,1),1);
        for k = 1:size(X,2)
            ks_(:,k) = ksdensity(X(:,k),'NumPoints',1000);
%             trapz(xi,ks_(:,k))
        end
        X = ks_;
        logX = log(X);
        logX(logX == -Inf) = 0;
        H = diag(-logX.' * X);
    elseif strcmp(method,'F-stat')
        Xmean = repmat(mean(X),size(X,1),1);
        Xdisp = (X - Xmean).^2 / (size(X,1)-1);
        H = sum(Xdisp.',2);
    elseif strcmp(method,'CoV')
        H = (std(X)./mean(X)).';
    end
end

% %% Plot long vs. short FI scatterplot
% matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles';
% load([matdir filesep 'Tasic_Ontology.mat'],'tasic_ontology_cns');
% load([matdir filesep 'Zeisel_Ontology.mat'],'zeisel_ontology_cns');
% load([matdir filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_short','tasic_fi_long')
% load([matdir filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_short','zeisel_fi_long')
% savenclose = 0;
% 
% tasic_fi_short_mean = mean(tasic_fi_short,2);
% tasic_fi_short_mean_norm = (tasic_fi_short_mean - min(tasic_fi_short_mean))...
%     /(max(tasic_fi_short_mean) - min(tasic_fi_short_mean));
% tasic_fi_long_mean = mean(tasic_fi_long,2);
% tasic_fi_long_mean_norm = (tasic_fi_long_mean - min(tasic_fi_long_mean))...
%     /(max(tasic_fi_long_mean) - min(tasic_fi_long_mean));
% 
% zeisel_fi_short_mean = mean(zeisel_fi_short,2);
% zeisel_fi_short_mean_norm = (zeisel_fi_short_mean - min(zeisel_fi_short_mean))...
%     /(max(zeisel_fi_short_mean) - min(zeisel_fi_short_mean));
% zeisel_fi_long_mean = mean(zeisel_fi_long,2);
% zeisel_fi_long_mean_norm = (zeisel_fi_long_mean - min(zeisel_fi_long_mean))...
%     /(max(zeisel_fi_long_mean) - min(zeisel_fi_long_mean));
% 
% teglu_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Neo Glu')).';
% oligo_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Oligo')).';
% vasc_ind_tasic = find(ismember(tasic_ontology_cns(:,2),'Vasc')).';
% other_ind_tasic = setdiff(1:length(tasic_ontology_cns(:,2)),[teglu_ind_tasic,oligo_ind_tasic,vasc_ind_tasic]);
% teglu_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Hip Neo Glu')).';
% oligo_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Oligo')).';
% vasc_ind_zeisel = find(ismember(zeisel_ontology_cns(:,2),'Vasc')).';
% other_ind_zeisel = setdiff(1:length(zeisel_ontology_cns(:,2)),[teglu_ind_zeisel,oligo_ind_zeisel,vasc_ind_zeisel]);
% 
% scatter_vals = cell(2,4);
% scatter_c = [[1 0 0]; [0 0 1]; [1 0 1]; [0.85 0.85 0.85]];
% scatter_s = {'o','o'};
% 
% scatter_vals{1,1} = [tasic_fi_short_mean_norm(teglu_ind_tasic),...
%     tasic_fi_long_mean_norm(teglu_ind_tasic)];
% scatter_vals{1,2} = [tasic_fi_short_mean_norm(oligo_ind_tasic),...
%     tasic_fi_long_mean_norm(oligo_ind_tasic)];
% scatter_vals{1,3} = [tasic_fi_short_mean_norm(vasc_ind_tasic),...
%     tasic_fi_long_mean_norm(vasc_ind_tasic)];
% scatter_vals{1,4} = [tasic_fi_short_mean_norm(other_ind_tasic),...
%     tasic_fi_long_mean_norm(other_ind_tasic)];
% scatter_vals{2,1} = [zeisel_fi_short_mean_norm(teglu_ind_zeisel),...
%     zeisel_fi_long_mean_norm(teglu_ind_zeisel)];
% scatter_vals{2,2} = [zeisel_fi_short_mean_norm(oligo_ind_zeisel),...
%     zeisel_fi_long_mean_norm(oligo_ind_zeisel)];
% scatter_vals{2,3} = [zeisel_fi_short_mean_norm(vasc_ind_zeisel),...
%     zeisel_fi_long_mean_norm(vasc_ind_zeisel)];
% scatter_vals{2,4} = [zeisel_fi_short_mean_norm(other_ind_zeisel),...
%     zeisel_fi_long_mean_norm(other_ind_zeisel)];
% 
% figure('Position',[0 0 500 400]); hold on;
% s1 = scatter(scatter_vals{1,4}(:,1),scatter_vals{1,4}(:,2),75,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(4,:),'MarkerFaceColor',scatter_c(4,:));
% scatter(scatter_vals{2,4}(:,1),scatter_vals{2,4}(:,2),75,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(4,:),'MarkerFaceColor',scatter_c(4,:));
% s2 = scatter(scatter_vals{1,1}(:,1),scatter_vals{1,1}(:,2),75,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(1,:),'MarkerFaceColor',scatter_c(1,:));
% s3 = scatter(scatter_vals{1,2}(1),scatter_vals{1,2}(2),75,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(2,:),'MarkerFaceColor',scatter_c(2,:));
% s4 = scatter(scatter_vals{1,3}(1),scatter_vals{1,3}(2),75,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(3,:),'MarkerFaceColor',scatter_c(3,:));
% scatter(scatter_vals{2,1}(:,1),scatter_vals{2,1}(:,2),75,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(1,:),'MarkerFaceColor',scatter_c(1,:));
% scatter(scatter_vals{2,2}(:,1),scatter_vals{2,2}(:,2),75,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(2,:),'MarkerFaceColor',scatter_c(2,:));
% scatter(scatter_vals{2,3}(:,1),scatter_vals{2,3}(:,2),75,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(3,:),'MarkerFaceColor',scatter_c(3,:));
% ylim([0 1]); yticks(1); xlim([0 1]); xticks([0 1]);
% ylabel('Long-range F.I.'); xlabel('Short-range F.I.'); 
% legend([s2, s3, s4, s1], {'Tel Glu','Oligo','Vasc','Other'},'Location','north')
% % text(0.05,0.925,'Tel Glu','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(1,:));
% % text(0.65,0.5,'Oligo','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(2,:));
% % text(0.9,0.2,'Vasc','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(3,:));
% set(gca,'box','on','FontName','DejaVunSans','FontSize',20);
% 
% % testshort = []; testlong = [];
% % for i = 1:2
% %     for j = 1:4
% % %         if size(scatter_vals{i,j},1) == 1
% % %             testshort = [testshort; scatter_vals{i,j}(1)];
% % %             testlong = [testlong; scatter_vals{i,j}(2)];
% % %         else
% %             testshort = [testshort; scatter_vals{i,j}(:,1)];
% %             testlong = [testlong; scatter_vals{i,j}(:,2)];
% % %         end
% %     end
% % end
% % [r,p] = corr(testshort,testlong)
% %% Plot long vs. short FI scatterplot 3 (supertypes)
% matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles';
% load([matdir filesep 'Tasic_Ontology.mat'],'tasic_ontology_cns');
% load([matdir filesep 'Zeisel_Ontology.mat'],'zeisel_ontology_cns');
% load([matdir filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_short','tasic_fi_long')
% load([matdir filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_short','zeisel_fi_long')
% savenclose = 0;
% 
% zeisel_fi_long_class_mean_norm = (zeisel_fi_long_class_mean - min(zeisel_fi_long_class_mean))...
%     /(max(zeisel_fi_long_class_mean) - min(zeisel_fi_long_class_mean));
% zeisel_fi_short_class_mean_norm = (zeisel_fi_short_class_mean - min(zeisel_fi_short_class_mean))...
%     /(max(zeisel_fi_short_class_mean) - min(zeisel_fi_short_class_mean));
% tasic_fi_long_class_mean_norm = (tasic_fi_long_class_mean - min(tasic_fi_long_class_mean))...
%     /(max(tasic_fi_long_class_mean) - min(tasic_fi_long_class_mean));
% tasic_fi_short_class_mean_norm = (tasic_fi_short_class_mean - min(tasic_fi_short_class_mean))...
%     /(max(tasic_fi_short_class_mean) - min(tasic_fi_short_class_mean));
% 
% tasic_classes_longshort = tasic_classes;
% % tasic_long_low = (tasic_fi_long_class_mean_norm < mean(tasic_fi_long_class_mean_norm));
% % tasic_short_low = (tasic_fi_short_class_mean_norm < mean(tasic_fi_short_class_mean_norm));
% % tasic_lowfi = ((tasic_long_low + tasic_short_low) == 2);
% % tasic_classes_longshort(tasic_lowfi) = [];
% % tasic_fi_long_class_mean_norm(tasic_lowfi) = [];
% % tasic_fi_short_class_mean_norm(tasic_lowfi) = [];
% zeisel_classes_longshort = zeisel_classes;
% % zeisel_long_low = (zeisel_fi_long_mean_norm < mean(zeisel_fi_long_class_mean_norm));
% % zeisel_short_low = (zeisel_fi_short_class_mean_norm < mean(zeisel_fi_short_class_mean_norm));
% % zeisel_lowfi = ((zeisel_long_low + zeisel_short_low) == 2);
% % zeisel_classes_longshort(zeisel_lowfi) = [];
% % zeisel_fi_long_class_mean_norm(zeisel_lowfi) = [];
% % zeisel_fi_short_class_mean_norm(zeisel_lowfi) = [];
% 
% 
% scatter_vals = cell(2,4);
% scatter_c = [[1 0 0]; [0 0 1]; [1 0 1]; [0.5 0.5 0.5]];
% scatter_s = {'o','o'};
% 
% teglu_ind_tasic = find(ismember(tasic_classes_longshort,'Neo Glu'));
% oligo_ind_tasic = find(ismember(tasic_classes_longshort,'Oligo'));
% vasc_ind_tasic = find(ismember(tasic_classes_longshort,'Vasc'));
% other_ind_tasic = setdiff(1:length(tasic_classes_longshort),[teglu_ind_tasic,oligo_ind_tasic,vasc_ind_tasic]);
% teglu_ind_zeisel = find(ismember(zeisel_classes_longshort,'Hip Neo Glu'));
% oligo_ind_zeisel = find(ismember(zeisel_classes_longshort,'Oligo'));
% vasc_ind_zeisel = find(ismember(zeisel_classes_longshort,'Vasc'));
% other_ind_zeisel = setdiff(1:length(zeisel_classes_longshort),[teglu_ind_zeisel,oligo_ind_zeisel,vasc_ind_zeisel]);
% 
% scatter_vals{1,1} = [tasic_fi_short_class_mean_norm(teglu_ind_tasic),...
%     tasic_fi_short_mean_norm(teglu_ind_tasic)];
% scatter_vals{1,2} = [tasic_fi_short_class_mean_norm(oligo_ind_tasic),...
%     tasic_fi_short_mean_norm(oligo_ind_tasic)];
% scatter_vals{1,3} = [tasic_fi_short_class_mean_norm(vasc_ind_tasic),...
%     tasic_fi_long_class_mean_norm(vasc_ind_tasic)];
% scatter_vals{1,4} = [tasic_fi_short_class_mean_norm(other_ind_tasic),...
%     tasic_fi_long_class_mean_norm(other_ind_tasic)];
% scatter_vals{2,1} = [zeisel_fi_short_class_mean_norm(teglu_ind_zeisel),...
%     zeisel_fi_long_class_mean_norm(teglu_ind_zeisel)];
% scatter_vals{2,2} = [zeisel_fi_short_class_mean_norm(oligo_ind_zeisel),...
%     zeisel_fi_long_class_mean_norm(oligo_ind_zeisel)];
% scatter_vals{2,3} = [zeisel_fi_short_class_mean_norm(vasc_ind_zeisel),...
%     zeisel_fi_long_class_mean_norm(vasc_ind_zeisel)];
% scatter_vals{2,4} = [zeisel_fi_short_class_mean_norm(other_ind_zeisel),...
%     zeisel_fi_long_class_mean_norm(other_ind_zeisel)];
% 
% figure('Position',[0 0 600 500]); hold on;
% scatter(scatter_vals{1,1}(1),scatter_vals{1,1}(2),100,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(1,:),'MarkerFaceColor',scatter_c(1,:));
% scatter(scatter_vals{1,2}(1),scatter_vals{1,2}(2),100,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(2,:),'MarkerFaceColor',scatter_c(2,:));
% scatter(scatter_vals{1,3}(1),scatter_vals{1,3}(2),100,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(3,:),'MarkerFaceColor',scatter_c(3,:));
% scatter(scatter_vals{1,4}(:,1),scatter_vals{1,4}(:,2),100,scatter_s{1},'filled',...
%     'MarkerEdgeColor',scatter_c(4,:),'MarkerFaceColor',scatter_c(4,:));
% scatter(scatter_vals{2,1}(1),scatter_vals{2,1}(2),100,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(1,:),'MarkerFaceColor',scatter_c(1,:));
% scatter(scatter_vals{2,2}(1),scatter_vals{2,2}(2),100,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(2,:),'MarkerFaceColor',scatter_c(2,:));
% scatter(scatter_vals{2,3}(1),scatter_vals{2,3}(2),100,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(3,:),'MarkerFaceColor',scatter_c(3,:));
% scatter(scatter_vals{2,4}(:,1),scatter_vals{2,4}(:,2),100,scatter_s{2},'filled',...
%     'MarkerEdgeColor',scatter_c(4,:),'MarkerFaceColor',scatter_c(4,:));
% ylim([0 1]); yticks(1); xlim([0 1]); xticks([0 1]);
% ylabel('Long-range F.I.'); xlabel('Short-range F.I.'); 
% legend({'Tel Glu','Oligo','Vasc'},'Location','northeast')
% % text(0.05,0.925,'Tel Glu','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(1,:));
% % text(0.65,0.5,'Oligo','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(2,:));
% % text(0.9,0.2,'Vasc','FontSize',18,'FontName','DejaVuSans','Color',scatter_c(3,:));
% set(gca,'box','on','FontName','DejaVunSans','FontSize',20);
% if savenclose
%     print([figdir filesep 'FI_AllCTScatter2'],'-dtiffn','-r300'); close;
% end
