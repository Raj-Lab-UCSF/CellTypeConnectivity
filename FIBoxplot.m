function [xpos_cell,g_cell,fi_reorder_cell] = FIBoxplot(matdir_,study_,types_,savenclose_,figdir_)

rng(0)
% if ~isequal(types_,'Classification')

fi_reorder_cell = cell(1,2); class_reorder_cell = cell(1,2);
fi_reorder_inds_cell = cell(1,2); cmap_cell = cell(1,2);
if strcmp('Both',study_)
    studycell = {'Zeisel','Tasic'};
    switch types_
        case 'All'
            titlecell = {'Zeisel, et al. - All Types', 'Tasic, et al. - All Types'};
            types_ = {'All','All'};  
            figstr = [figdir_ filesep 'FI_AllConnections'];
        case 'Classification'
            titlecell = {'Zeisel, et al. - Binary Connections', 'Tasic, et al. - Binary Connections'};
            types_ = {'Classification','Classification'};  
            figstr = [figdir_ filesep 'FI_Classification'];
    end
    xoffsets = [-0.35 -0.17];
    for i = 1:2
        [fi_reorder,class_reorder,fi_reorder_inds] = FISuperclass(matdir_,studycell{i},types_{i});
        fi_reorder_cell{i} = fi_reorder; class_reorder_cell{i} = class_reorder;
        fi_reorder_inds_cell{i} = fi_reorder_inds;
        cmap_cell{i} = hsv(length(class_reorder));
    end
else
    figstr = [figdir_ filesep 'FI_' study_ types_{1} types_{2}];
    titlecell = cell(1,2);
    if strcmp(study_,'Zeisel')
        if isequal(types_,{'Out','In'})
            xoffsets = [-0.35 -0.35];
        else
            xoffsets = [-0.475 -0.475];
        end
    else
        if isequal(types_,{'Out','In'})
            xoffsets = [-0.17 -0.17];
        else
            xoffsets = [-0.25, -0.25];
        end
    end
    for i = 1:2
        type = types_{i};
        [fi_reorder,class_reorder,fi_reorder_inds] = FISuperclass(matdir_,study_,type);
        fi_reorder_cell{i} = fi_reorder; class_reorder_cell{i} = class_reorder;
        fi_reorder_inds_cell{i} = fi_reorder_inds;
        cmap_cell{i} = hsv(length(class_reorder));
        switch type
            case 'Long'
                titlecell{i} = sprintf('%s et al. - Long-range',study_);
            case 'Short'
                titlecell{i} = sprintf('%s et al. - Short-range',study_);
            case 'Out'
                titlecell{i} = sprintf('%s et al. - Source Types',study_);
            case 'In'
                titlecell{i} = sprintf('%s et al. - Target Types',study_);
            case 'Out Long'
                titlecell{i} = sprintf('%s et al. - Source Types, Long-range',study_);
            case 'Out Short'
                titlecell{i} = sprintf('%s et al. - Source Types, Short-range',study_);
            case 'In Long'
                titlecell{i} = sprintf('%s et al. - Target Types, Long-range',study_);
            case 'In Short'
                titlecell{i} = sprintf('%s et al. - Target Types, Short-range',study_);
            case 'NeoToAll Out'
                titlecell{i} = sprintf('%s et al. - Source, Neocortical/Other',study_);
            case 'NeoToAll In'
                titlecell{i} = sprintf('%s et al. - Target, Neocortical/Other',study_);
        end
    end
end

xpos_cell = cell(1,2); g_cell = cell(1,2);
xposscatter = @(x,y) 0.15 * (2*rand(length(x),1) - 1) + y;
for i = 1:length(xpos_cell)
    fi_reorder = fi_reorder_cell{i}.';
    xpos = zeros(size(fi_reorder)); xpos = xpos(:); g = xpos;
    for j = 1:size(fi_reorder,2)
        curinds = (1:size(fi_reorder,1)) + (j-1)*size(fi_reorder,1);
        g(curinds) = j;
        xpos(curinds) = xposscatter(curinds,j);
    end
    xpos_cell{i} = xpos; g_cell{i} = g;
end
if strcmp('Both',study_)
    figure('Position',[0 0 1000 500]);
    fontsize = 18;
elseif strcmp('Zeisel',study_) && isequal(types_,{'Out','In'})
    figure('Position',[0 0 1000 525]);
    fontsize = 18;
elseif strcmp('Tasic',study_) && isequal(types_,{'Out','In'})
    figure('Position',[0 0 1000 475]);
    fontsize = 18;
else
    figure('Position',[0 0 650 450]);
    fontsize = 16;
end
tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(xpos_cell)
    fi_reorder = fi_reorder_cell{i}.';
    fi_reorder_inds = fi_reorder_inds_cell{i};
    class_reorder = class_reorder_cell{i};
    cmap = cmap_cell{i}; cmap_reorder = cmap(fi_reorder_inds,:);
    nexttile; hold on;
    warning('off', 'stats:boxplot:BadObjectType');
    b = boxplot(fi_reorder,'Colors',cmap_reorder,'Symbol','');
    set(b,{'linew'},{1})
    gscatter(xpos_cell{i},fi_reorder(:),g_cell{i},cmap_reorder,[],15,'off');
    ylim_ = max(fi_reorder(:));
    ylim([0,1.1*ylim_]);
    ylabel('F.I.','FontSize',fontsize,'FontName','DejaVuSans');
    xlabel('');
    set(gca,'XTick',1:length(class_reorder),'XTickLabel',class_reorder,...
        'XTickLabelRotation',40,'FontName','DejaVuSans','FontSize',fontsize,...
        'TickLength',[0 0],'YTick',ylim_);
    ytickformat('%.2f')
    text(xoffsets(i),0,'0','FontSize',fontsize,'FontName','DejaVuSans')
    title(titlecell{i},'FontSize',20,...
        'FontName','DejaVuSans','FontWeight','bold');
end
% else
%     figstr = [figdir_ filesep 'FI_AllConnections_Classification'];
%     [fi_reorder,class_reorder,fi_reorder_inds] = FISuperclass(matdir_,study_,types_);
%     xposscatter = @(x,y) 0.15 * (2*rand(length(x),1) - 1) + y;
%     fi_reorder = fi_reorder.';
%     cmap = hsv(length(class_reorder)); cmap_reorder = cmap(fi_reorder_inds,:);
%     xpos = zeros(size(fi_reorder)); xpos = xpos(:); g = xpos;
%     for j = 1:size(fi_reorder,2)
%         curinds = (1:size(fi_reorder,1)) + (j-1)*size(fi_reorder,1);
%         g(curinds) = j;
%         xpos(curinds) = xposscatter(curinds,j);
%     end
%     figure('Position',[0 0 1000 500]);
%     nexttile; hold on;
%     warning('off', 'stats:boxplot:BadObjectType');
%     b = boxplot(fi_reorder,'Colors',cmap_reorder,'Symbol','');
%     set(b,{'linew'},{1})
%     gscatter(xpos,fi_reorder(:),g,cmap_reorder,[],15,'off');
%     ylim_ = max(fi_reorder(:));
%     ylim([0,1.1*ylim_]);
%     ylabel('Feature Importance','FontSize',18,'FontName','DejaVuSans');
%     xlabel('');
%     set(gca,'XTick',1:length(class_reorder),'XTickLabel',class_reorder,...
%         'XTickLabelRotation',35,'FontName','DejaVuSans','FontSize',18,...
%         'TickLength',[0 0],'YTick',ylim_);
%     ytickformat('%.1f')
%     text(-0.17,0,'0','FontSize',18,'FontName','DejaVuSans')
%     title('Tasic, et al. - Binary Connections','FontSize',20,...
%         'FontName','DejaVuSans','FontWeight','bold');
% end

if savenclose_
    print(figstr,'-dtiffn','-r300'); close;
end
end