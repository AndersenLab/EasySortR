#' Test for correlation among traits in a data frame and output a plot
#' 
#' @param df A dataframe output by easysorter
#' @param gropus Number of groups to split the correlation into, usually based on number of significant PCs
#' @param phenocol When spread, the column number on which traits begin, usually 11 with HTA data.
#' @return A dendrogram and a correlation plot for the data
#' @export

corplot <- function(df, groups, phenocol = NA){
    df <- df %>%
        spread(trait, phenotype) %>%
        ungroup() 
    if(is.na(phenocol)){
        df2 <- df[,11:ncol(df)]
    } else {
        df2 <- df[,phenocolo:ncol(df)]
    }
    
    reorder_cormat <- function(cormat){
        # Use correlation between variables as distance
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd) 
        cormat <-cormat[hc$order, hc$order]
    }
    
    traitcormat <- stats::cor(df2, use = "pairwise.complete.obs")
    reordered <- reorder_cormat(traitcormat)
    
    dd <- stats::as.dist((1-traitcormat)/2)
    hc <- stats::hclust(dd)
    
    dend <- ggdendro::ggdendrogram(hc, labels = FALSE, leaf_labels = FALSE) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = .5),
                       plot.margin = grid::unit(c(0,0,0,0), "cm"))
    
    #clusters1 is the Var1 trait and which cluster it is in:
    #first divide the hc tree into the number of groups you want based on how many significant PCs we have:
    clusters1 <- data.frame(level1 = stats::cutree(hc, k = groups)) %>%
        dplyr::mutate(Var1 = rownames(.))
    #then order the traits based on the way they appear in the dendrogram:
    #extract dendrogram info
    hc2 <- ggdendro::dendro_data(hc)$labels
    #factor so you can order them in the heatmap based on this order
    hc2$label <- factor(hc2$label)
    #order the traits based on their order in the dendrogram
    clusters1$Var1 <- factor(clusters1$Var1, levels = hc2$label)
    #assign cluster groups based on this order that you want the traits (they were grouped by cutree, but group number is not how we want it)
    clusters1 <- clusters1 %>%
        dplyr::arrange(Var1) %>%
        dplyr::group_by(level1)%>%
        dplyr::mutate(newlev = 1)
    for(i in 2:nrow(clusters1)){
        if(clusters1$level1[i] == clusters1$level1[i-1]){
            clusters1$newlev[i] <- clusters1$newlev[i-1]
        } else clusters1$newlev[i] <- clusters1$newlev[i-1]+1
    }
    #clusters2 is the same thing, but for clustering the y-axis traits in the heatmap
    clusters2 <- clusters1 %>%
        dplyr::rename(newlev2 = newlev, Var2 = Var1) %>%
        dplyr::ungroup()%>%
        dplyr::select(-level1)
    #add correlations and the grouping data altogether
    melted <- data.table::melt(reordered) %>%
        dplyr::left_join(clusters1)%>%
        dplyr::left_join(clusters2)
    #reorder Var2 so white line goes diagonal across whole heat map, not just within facet
    melted$Var2 <- factor(melted$Var2, levels = rev(clusters2$Var2))
    
    traitcorplot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var1, y = Var2, fill=value^2))+
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(colours=terrain.colors(6), name=expression(bold(italic("r")^2)))+
        ggplot2::scale_x_discrete(position = "top")+
        ggplot2::facet_grid(newlev2~newlev, scales = "free", space = "free", switch = "x")+
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=8, color="black",angle = 90, hjust=0),
                       axis.text.y = ggplot2::element_text(size=8, color="black"),
                       panel.background = ggplot2::element_rect(fill = '#E5E8E8'),
                       strip.text = ggplot2::element_blank(),
                       strip.text.y = ggplot2::element_blank(),
                       panel.spacing = grid::unit(1, "pt"),
                       plot.margin = grid::unit(c(-1,0,0,0),"cm"))+
        ggplot2::labs(title="", x=NULL, y=NULL) 
    
    
    finaldend <- cowplot::plot_grid(NULL, dend, NULL, nrow = 1, ncol = 3, rel_widths = c(.12, 1, .13))
    finalplot <- cowplot::plot_grid(finaldend, traitcorplot, nrow = 2, ncol = 1, rel_heights = c(.5, 1))
    return(finalplot)
}  


#' Perform principal component analysis on traits
#' 
#' @param df A regressed dataframe of HTA traits for PCA
#' @return pca_obj, a PCA object output from princomp; 
#' loadings, a vector of loadings for each HTA trait into each PC;
#' cumsums, a matrix of cumulative sums of variance explained by each PC;
#' RIAILPCphenos, a dataframe of PC values for each strain;
#' corr_PC_trait, correlation of each trait to each PC

runPCA <- function(df){
    
    traits_for_PCA <- df %>%
        dplyr::filter(!strain %in% c("N2","CB4856")) %>% #remove parents from the dataset
        dplyr::select(trait, strain, phenotype)%>% #select columns for PCA
        tidyr::spread(trait, phenotype) %>% #spread so each row is a strain
        na.omit()
    
    row.names(traits_for_PCA) <- traits_for_PCA$strain
    
    pc_traits <- traits_for_PCA %>%
        dplyr::select(-strain)
    
    #scale the phenotypes before running princomp
    scales_pc_traits <- as.data.frame(scale(pc_traits))
    
    pca_obj <- princomp(scales_pc_traits)
    
    #pull the loadings (amount of each trait in each PC)
    loadings <- loadings(pca_obj)[]
    
    #pull the total variance explained with each PC and call it "drug.cumsum"
    cumsums <- t(as.matrix(cumsum(pca_obj$sdev^2/sum(pca_obj$sdev^2))))
    
    #keep the PCs that explain at least 95% of the variation in phenotype
    keep <- data.frame("Val" = cumsums[1,], "PC" = 1:ncol(cumsums)) %>%
        dplyr::filter(Val > .95)
    keep <- as.numeric(min(keep$PC))
    
    #make df of all loadings for the first five PCs
    pcaoutput <- as.data.frame(loadings) %>%
        dplyr::mutate(trait = rownames(.)) %>%
        tidyr::gather(component, variance, -trait) %>%
        dplyr::mutate(component = as.numeric(stringr::str_split_fixed(component, "Comp.", 2)[,2])) %>%
        dplyr::filter(component <= keep) %>%
        dplyr::mutate(component = paste0("PC", component))
    
    #calculate the phenotypic values of each strain for each of the top five PCs for mapping!
    RIAIL_PCphenos <- as.data.frame(pca_obj$scores) %>%
        dplyr::mutate(strain = rownames(.), condition = "bleomycin")%>%
        tidyr::gather(trait, phenotype, -strain, -condition) %>%
        dplyr::filter(trait %in% c("Comp.1", "Comp.2","Comp.3","Comp.4","Comp.5"))
    
    ### Correlation of each PC with each of the HTA traits:
    spread_PCA <- RIAIL_PCphenos %>%
        tidyr::spread(trait, phenotype) %>%
        dplyr::select(-condition)
    comparephenos <- merge(spread_PCA, traits_for_PCA, by = "strain")%>%
        na.omit() %>%
        ungroup()%>%
        dplyr::select(-strain)
    
    corr_PC_trait <- cor(comparephenos)
    
    corr_PC_trait <- corr_PC_trait[1:keep,(keep+1):31]
    
    return(list(pca_obj, loadings, RIAIL_PCphenos, corr_PC_trait))
}