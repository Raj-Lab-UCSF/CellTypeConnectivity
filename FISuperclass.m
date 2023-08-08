function [fi_reorder,class_reorder,fi_reorderinds] = FISuperclass(matdir_,study_,type_)

if strcmp(study_,'Zeisel')
    load([matdir_ filesep 'Zeisel_Ontology.mat'],'zeisel_ontology_cns');
    zeisel_classes = unique(zeisel_ontology_cns(:,2));
    zeisel_fi_class = zeros(length(zeisel_classes),10); % k = 10
    switch type_
        case 'All'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_all');
            zeisel_fi = zeisel_fi_all;
        case 'Long'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_long');
            zeisel_fi = zeisel_fi_long;
        case 'Short'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_short');
            zeisel_fi = zeisel_fi_short;
        case 'Out'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_out');
            zeisel_fi = zeisel_fi_out;
        case 'In'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_in');
            zeisel_fi = zeisel_fi_in;
        case 'Out Long'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_outlong');
            zeisel_fi = zeisel_fi_outlong;
        case 'Out Short'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_outshort');
            zeisel_fi = zeisel_fi_outshort;
        case 'In Long'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_inlong');
            zeisel_fi = zeisel_fi_inlong;
        case 'In Short'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_inshort');
            zeisel_fi = zeisel_fi_inshort;
        case 'NeoToTha'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_tha');
            zeisel_fi = zeisel_fi_ncx_tha;
        case 'NeoToTha Out'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_tha_out');
            zeisel_fi = zeisel_fi_ncx_tha_out;
        case 'NeoToTha In'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_tha_in');
            zeisel_fi = zeisel_fi_ncx_tha_in;
        case 'NeoToNeo'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_ncx');
            zeisel_fi = zeisel_fi_ncx_ncx;
        case 'NeoToNeo Out'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_ncx_out');
            zeisel_fi = zeisel_fi_ncx_ncx_out;
        case 'NeoToNeo In'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_ncx_in');
            zeisel_fi = zeisel_fi_ncx_ncx_in;
        case 'NeoToAll Out'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_all_out');
            zeisel_fi = zeisel_fi_ncx_all_out;
        case 'NeoToAll In'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_ncx_all_in');
            zeisel_fi = zeisel_fi_ncx_all_in;
        case 'Classification'
            load([matdir_ filesep 'Zeisel_FeatureImportance.mat'],'zeisel_fi_classification');
            zeisel_fi = zeisel_fi_classification;
    end
    for i = 1:length(zeisel_classes)
        classinds = ismember(zeisel_ontology_cns(:,2),zeisel_classes{i});
        zeisel_fi_class(i,:) = mean(zeisel_fi(classinds,:));
    end
    zeisel_fi_class_mean = mean(zeisel_fi_class,2);
    [~,fi_reorderinds] = sort(zeisel_fi_class_mean,'descend');
    fi_reorder = zeisel_fi_class(fi_reorderinds,:);
    class_reorder = zeisel_classes(fi_reorderinds);
elseif strcmp(study_,'Tasic')
    load([matdir_ filesep 'Tasic_Ontology.mat'],'tasic_ontology_cns');
    tasic_classes = unique(tasic_ontology_cns(:,2));
    tasic_fi_class = zeros(length(tasic_classes),10); % k = 10
    switch type_
        case 'All'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_all');
            tasic_fi = tasic_fi_all;
        case 'Long'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_long');
            tasic_fi = tasic_fi_long;
        case 'Short'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_short');
            tasic_fi = tasic_fi_short;
        case 'Out'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_out');
            tasic_fi = tasic_fi_out;
        case 'In'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_in');
            tasic_fi = tasic_fi_in;
        case 'Out Long'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_outlong');
            tasic_fi = tasic_fi_outlong;
        case 'Out Short'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_outshort');
            tasic_fi = tasic_fi_outshort;
        case 'In Long'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_inlong');
            tasic_fi = tasic_fi_inlong;
        case 'In Short'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_inshort');
            tasic_fi = tasic_fi_inshort;
        case 'Classification'
            load([matdir_ filesep 'Tasic_FeatureImportance.mat'],'tasic_fi_classification');
            tasic_fi = tasic_fi_classification;
    end
    for i = 1:length(tasic_classes)
        classinds = ismember(tasic_ontology_cns(:,2),tasic_classes{i});
        tasic_fi_class(i,:) = mean(tasic_fi(classinds,:));
    end
    tasic_fi_class_mean = mean(tasic_fi_class,2);
    [~,fi_reorderinds] = sort(tasic_fi_class_mean,'descend');
    fi_reorder = tasic_fi_class(fi_reorderinds,:);
    class_reorder = tasic_classes(fi_reorderinds);
end



end