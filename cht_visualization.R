args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)






###############################################
#  Boxplot of number of egenes as a function of number of pcs
#  Each point in boxplot is a time step
egenes_as_function_of_pcs_boxplot <- function(cht_output_dir, parameter_string, output_file) {
    # Extract data
    egenes <- c()
    time_step <- c()
    num_pcs <- c()

    # loop through time steps and pcs
    for (temp_time_step in 0:15) {
        for (temp_num_pc in 0:5) {
            # egene file for this time step and this pc
            egene_file <- paste0(cht_output_dir,"cht_results_",parameter_string,"_num_pc_", temp_num_pc,"_time_",temp_time_step,"_qval_.1_significant_egenes.txt")
            egene_data <- read.table(egene_file, header=TRUE)
            # Get number of egenes at this time step and pc num
            num_egenes <- dim(egene_data)[1]
            # Now store data 
            egenes <- c(egenes,num_egenes)
            time_step <- c(time_step, temp_time_step)
            num_pcs <- c(num_pcs, temp_num_pc)
        }
    }
    # Put data in organized data frame
    df <- data.frame(egenes=as.numeric(egenes),time_step = factor(time_step), num_pcs = factor(num_pcs))
    # PLOT!!
    box_plot <- ggplot(df, aes(x=num_pcs, y=egenes)) + geom_boxplot(width=.54)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(x = "Number of PCs", y = "Number of eGenes")
    ggsave(box_plot, file=output_file,width = 20,height=10.5,units="cm")
}



read_in_summary_statistic_data <- function(file_name) {
    xx <- read.table(file_name)
    xxx <- xx[,2:dim(xx)[2]]
    return(xxx)
}

# Plot correlation histogram for summary stat
correlation_heatmap <- function(correlation_matrix, output_file, version,num_pc) {
    colnames(correlation_matrix) <- paste0(0:15)
    rownames(correlation_matrix) <- paste0(0:15)
    melted_corr <- melt(correlation_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu")
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 0, vjust=.5))
    heatmap <- heatmap + labs(x = "Time Step", y = "Time Step", title=paste0("Num PCs = ",num_pc," / Statistic = ", version), fill= "Spearman Rho")

    ggsave(heatmap, file=output_file,width = 18,height=13.5,units="cm")


}




###############################################
# Heatmap of correlation of summary statistics between time steps
# Ie. heatmap is of dimension number of time steps by number of time steps
# Do independently for each number of PCs
summary_statistic_correlation_heatmap <- function(parameter_string, pc_num, cht_output_dir, cht_visualization_dir) {
    # All summary statistic files have the following file stem
    file_stem <- paste0(cht_output_dir,"best_variant_per_egene_",parameter_string, "_num_pc_",pc_num)

    # Specific summmary statistic files
    alpha_file <- paste0(file_stem,"_alpha.txt")
    beta_file <- paste0(file_stem, "_beta.txt")
    pvalue_file <- paste0(file_stem, "_pvalues.txt")

    # Load in summary stat data
    # Of dimension (number of eGenes) X (number of time steps)
    alpha <- read_in_summary_statistic_data(alpha_file)
    beta <- read_in_summary_statistic_data(beta_file)
    pvalue <- read_in_summary_statistic_data(pvalue_file)
    # Compute allelic fraction (or p as it is referred to in the WASP paper)
    allelic_fraction <- alpha/(alpha + beta)

    # All output files have the following stem
    output_file_stem <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_correlation_heatmap_")

    # Plot correlation histogram for pvalue summary stat
    correlation_heatmap(cor(pvalue,method="spearman"), paste0(output_file_stem, "pvalues.png"), "Pvalue", pc_num)

    # Plot correlation histogram for pvalue summary stat
    correlation_heatmap(cor(alpha,method="spearman"), paste0(output_file_stem, "alpha.png"), "Alpha", pc_num)

    # Plot correlation histogram for pvalue summary stat
    correlation_heatmap(cor(beta,method="spearman"), paste0(output_file_stem, "beta.png"), "Beta", pc_num)

    # Plot correlation histogram for pvalue summary stat
    correlation_heatmap(cor(allelic_fraction,method="spearman"), paste0(output_file_stem, "allelic_fraction.png"), "Allelic Fraction",pc_num)

}
sample_non_significant_hits <- function(pvalues, fraction_kept=.01,fraction_sampled=.001) {
    index <- floor(length(pvalues)*fraction_kept)
    to_keep <- pvalues[1:index]
    to_filter <- pvalues[(index+1):length(pvalues)]
    filtered <- sort(sample(to_filter,floor(length(to_filter)*fraction_sampled)))
    return(c(to_keep,filtered))
}

###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to uniform
## 2. Permuted-data pvalues compared to uniform
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)

qq_plot_vs_uniform <- function(input_stem, null_stem, output_file) {
    # Make qq-plot for each of 16 time steps

    time_step <- 0
    p0 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 1
    p1 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 2
    p2 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
        
    time_step <- 4
    p3 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 4
    p4 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
        
    time_step <- 5
    p5 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 6
    p6 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
        
    time_step <- 7
    p7 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 8
    p8 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 9
    p9 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 10
    p10 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 11
    p11 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 12
    p12 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 13
    p13 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 14
    p14 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 15
    p15 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    # Merge those 16 plots with cowplot!

    pdf(output_file)
    gg <- plot_grid(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,nrow=4,ncol=4,label_size=8)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()
}

# Helper method to qq_plot_vs_uniform
# Makes a qqplot for a specific time step
make_qq_plot_vs_uniform_one_time_step <- function(real_eqtl_file, null_eqtl_file, time_step) {
    # Read in real data
    all_eqtl_nominal_pvalues <- read.table(real_eqtl_file, header=TRUE)
    # Extract pvalues
    pvalues <- sort(all_eqtl_nominal_pvalues$pvalue)
    # Simulate uniform distribution
    uniform_1 <- sort(runif(length(pvalues)))
    # Sample points (because impossible to plot ALL hits)
    pvalues <- sample_non_significant_hits(pvalues)
    
    # Read in null data
    null_data <- read.table(null_eqtl_file,header=TRUE)
    # Extract pvalues
    null_pvalues <- sort(null_data$pvalue)
    # Simulate uniform distribution
    uniform_2 <- sort(runif(length(null_pvalues)))
    # Sample points (because impossible to plot ALL hits)
    null_pvalues <- sample_non_significant_hits(null_pvalues)


    # Organize into data frame
    all_pvalues <- c(pvalues, null_pvalues)
    type <- c(rep("real",length(pvalues)),rep("null",length(null_pvalues)))

    uniform <- c(sample_non_significant_hits(uniform_1), sample_non_significant_hits(uniform_2))

    df <- data.frame(pvalues=-log10(all_pvalues + .000000000001), expected_pvalues=-log10(uniform + .000000000001), type=factor(type))

    # PLOT!
    max_val <-max(max(-log10(uniform + .000000000001)), max(-log10(all_pvalues + .000000000001)))
    #PLOT!
    scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues, colour = type)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(colour="data type",x = "Uniform", y = "Real", title = paste0("Time step ", time_step))
    scatter <- scatter + geom_abline() +  theme(legend.position="none")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    return(scatter)
}

###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to permuted
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)
qq_plot_vs_permuted <- function(input_stem, null_stem, output_file) {
    # Make qq-plot for each of 16 time steps
    time_step <- 0
    p0 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 1
    p1 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 2
    p2 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 3
    p3 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 4
    p4 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 5
    p5 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 6
    p6 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 7
    p7 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 8
    p8 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 9
    p9 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 10
    p10 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 11
    p11 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 12
    p12 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 13
    p13 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 14
    p14 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 15
    p15 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    # Merge those 16 plots with cowplot!
    pdf(output_file)
    gg <- plot_grid(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,nrow=4,ncol=4,label_size=8)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()

}

# Helper method to qq_plot_vs_permuted
# Makes a qqplot for a specific time step
make_qq_plot_vs_permuted_one_time_step <- function(real_eqtl_file, null_eqtl_file, time_step) {
    # Extract pvalues from real data
    all_eqtl_nominal_pvalues <- read.table(real_eqtl_file, header=TRUE)
    # Subselect hits (can't plot all hits)
    pvalues <- sample_non_significant_hits(sort(all_eqtl_nominal_pvalues$pvalue))

    # Extract pvalues from permuted data
    null_data <- read.table(null_eqtl_file,header=TRUE)
    # Subselect hits (can't plot all hits)
    null_pvalues <- sample_non_significant_hits(sort(null_data$pvalue))

    # put into convenient data frame
    df <- data.frame(real_pvalues=-log10(pvalues + .000000000001), expected_pvalues=-log10(null_pvalues + .000000000001))


    max_val <-max(max(-log10(pvalues + .000000000001)), max(-log10(null_pvalues + .000000000001)))
    #PLOT!
    scatter <- ggplot(df, aes(x = expected_pvalues, y = real_pvalues)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Permuted", y = "Real", title = paste0("Time step ", time_step))
    scatter <- scatter + geom_abline() +  theme(legend.position="none")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    return(scatter)
}




eqtl_comparison_to_reference_bar_plot <- function(cht_enrichment_dir,parameter_string, data_set_name, pc_num, plot_file) {
    # First extract data. And get into nice data format
    pvalues <- c()
    version <- c()
    time_step <- c()
    # Loop through time steps
    for (temp_time_step in 0:15) {
        # Get and parse enrichment file for this time step
        ipsc_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num,"_time_",temp_time_step,"_", data_set_name,"_real_v_matched_controls.txt")
        data <- read.table(ipsc_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("real", length(data$real_pvalue))))
        pvalues <- c(pvalues,data$matched_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$matched_pvalue)))
        version <- c(version, as.character(rep("matched", length(data$matched_pvalue))))
        #print(temp_time_step)
        #print(wilcox.test(data$real_pvalue,data$matched_pvalue))
    }
    df <- data.frame(pvalues = as.numeric(pvalues), version = factor(version,c("real","matched")), time_step = factor(time_step))

    box_plot <- ggplot(df, aes(x=time_step, y=pvalues, fill=version)) + geom_boxplot(width=.54)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Test version",x = "time step", y = "pvalue")
    ggsave(box_plot, file=plot_file,width = 20,height=10.5,units="cm")
}

###############################################
# Boxplot showing pvalues found in our data for only the eqtls in a specific data set
# 1 plot per data set
# 1 plot for pc
eqtl_comparison_to_background_shell <- function(pc_num, cht_visualization_dir, parameter_string, cht_enrichment_dir, eqtl_data_set_file) {
    # Load in data set data
    data_sets <- read.table(eqtl_data_set_file)
    num_data_sets <- dim(data_sets)[1]
    # Loop through each of data sets
    for (data_set_num in 1:num_data_sets) {
        # Name of data set of this line
        data_set_name <- paste0(data_sets[data_set_num, 1])
        # Output file
        plot_file <- paste0(cht_visualization_dir,parameter_string,"_num_pc_",pc_num,"_data_set_comparison_",data_set_name,"_boxplot.png")
        # Make boxplot
        eqtl_comparison_to_reference_bar_plot(cht_enrichment_dir, parameter_string, data_set_name, pc_num, plot_file)
    }
}

###############################################
# Boxplot showing pvalues found in our data based on eqtls found in GTEx v7:
## 1. HLV
## 2. skin_not_sun_exposed
## 3. stomach
## 4. Breast
# 1 plot for pc
gtex_tissue_comparison_bar_plot <- function(tissue_comparison_plot_file, pc_num, parameter_string, cht_enrichment_dir) {
    # First extract data. And get into nice data format
    pvalues <- c()
    version <- c()
    time_step <- c()
    for (temp_time_step in 0:15) {
        # HLV
        hlv_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_hlv_beta_filter_real_v_matched_controls.txt")
        data <- read.table(hlv_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("heart_left", length(data$real_pvalue))))
        # stomach
        stomach_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_stomach_beta_filter_real_v_matched_controls.txt")
        data <- read.table(stomach_file, header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("stomach", length(data$real_pvalue))))
        # breast
        breast_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_breast_mammary_beta_filter_real_v_matched_controls.txt")
        data <- read.table(breast_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("breast", length(data$real_pvalue))))
        # Skin not sun exposed
        skin_no_sun_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_skin_not_sun_exposed_beta_filter_real_v_matched_controls.txt")
        data <- read.table(skin_no_sun_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("skin_not_sun", length(data$real_pvalue))))  
        
    }
    # Put everything into data frame
    df <- data.frame(pvalues = -log10(as.numeric(pvalues + .000001)), version = factor(version,c("heart_left","breast","skin_not_sun","stomach")), time_step = factor(time_step))

    # PLOT!
    box_plot <- ggplot(df, aes(x=time_step, y=pvalues, fill=version)) + geom_boxplot(width=.7,outlier.size = 0.1)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Test version",x = "time step", y = "-log10(pvalue)")
    box_plot <- box_plot + theme(legend.position="bottom") 
    ggsave(box_plot, file=tissue_comparison_plot_file,width = 20,height=10.5,units="cm")


}


visualize_trajectories_for_one_cluster <- function(filtered_statistic_matrix, cluster_center, statistic_type, k) {
    filtered_statistic_matrix <- t(scale(t(filtered_statistic_matrix)))
    #filtered_statistic_matrix <- t(apply(filtered_statistic_matrix, 1, function(x)(x-min(x))/(max(x)-min(x))))

    group <- c()
    time_step <- c()
    statistics <- c()
    nrow <- dim(filtered_statistic_matrix)[1]
    for (row_num in 1:nrow) {
        time_step <- c(time_step,1:16)
        statistics <- c(statistics, as.numeric(filtered_statistic_matrix[row_num,]))
        group <- c(group, rep(row_num, 16))
    }
    df <- data.frame(statistic = statistics, time_step = time_step, group = factor(group))

    line_plot <- ggplot(df,aes(time_step,statistic,group=group))
    line_plot <- line_plot + geom_line(alpha=0.15)
    line_plot <- line_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    line_plot <- line_plot + labs(x = "time step", y = statistic_type, title = paste0("Cluster ", k))


    return(line_plot)
}

kmeans_cluster_of_summary_statistics <- function(statistic_matrix, output_file, pc_num, k, statistic_type) {
    # Run Kmeans clustering
    kmeans_obj <- kmeans(statistic_matrix, k,nstart=30,iter.max=40)
    # Get centers of kmeans clusters
    # Matrix is of dimension k X (time_steps)
    centers <- kmeans_obj$centers
    # Get assignments of each sample to 1 of k clusters
    cluster_assignments <- kmeans_obj$cluster




    # Loop through clusters
    cluster_num <- 1
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p1 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 2
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p2 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 3
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p3 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 4
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p4 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 5
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p5 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 6
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p6 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 7
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p7 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 8
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p8 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 9
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p9 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)


    pdf(output_file)
    gg <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=3,ncol=3,label_size=8)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()

}

visualize_number_of_genome_wide_significant_egenes <- function(input_stem, output_file) {
    num_genes <- c()
    for (temp_time_step in 0:15) {
        sig_egene_file <- paste0(input_stem, temp_time_step, "_qval_.1_significant_egenes.txt")
        data <- read.table(sig_egene_file,header=TRUE)
        time_step_num_egenes <- dim(data)[1]
        num_genes <- c(num_genes, time_step_num_egenes)
    }
    df <- data.frame(time_step=0:15, num_egenes=num_genes)
    p <- ggplot(df, aes(time_step, num_genes)) 
    p <- p + geom_bar(stat = "identity",aes(fill=time_step)) 
    p <- p + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p <- p + labs(x = "time step", y = "Number of egenes (FDR <= .1)", title="egenes found with WASP-CHT")
    p <- p + scale_fill_gradient(low="pink",high="blue")
    p <- p + theme(legend.position="none")
    ggsave(p, file=output_file,width = 20,height=10.5,units="cm")

}



parameter_string = args[1]  # string used to keep track of files used with specified parameter settting
cht_output_dir = args[2]  # input directory with cht test results
cht_visualization_dir = args[3]  # output directory to save images
cht_enrichment_dir = args[4]   # Input directory that has files showing our pvalues of eQTLs (according to various data sets)
eqtl_data_set_file = args[5]  # Input file containing info (1 line) for each of the data sets

###############################################
# Boxplot of number of egenes as a function of number of pcs
#  Each point in boxplot is a time step
output_file <- paste0(cht_visualization_dir, parameter_string, "egenes_as_function_of_pcs_boxplot.png")
# egenes_as_function_of_pcs_boxplot(cht_output_dir, parameter_string, output_file)






###############################################
# Bar plot showing number of genome wide significant egenes at each of the 16 time steps
# Do this for each of the pc_nums in dependently
for (pc_num in 0:5) {
    output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_number_of_significant_egenes_per_time_step_bar_plot.png")
    input_stem <- paste0(cht_output_dir, "cht_results_",parameter_string,"_num_pc_",pc_num,"_time_")
    #visualize_number_of_genome_wide_significant_egenes(input_stem, output_file)
}


###############################################
# Heatmap of correlation of summary statistics between time steps
# Ie. heatmap is of dimension number of time steps by number of time steps
# Do independently for each number of PCs
for (pc_num in 0:5) {
    #summary_statistic_correlation_heatmap(parameter_string, pc_num, cht_output_dir, cht_visualization_dir)
}


###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to uniform
## 2. Permuted-data pvalues compared to uniform
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)
for (pc_num in 0:5) {
    # Output file
    output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_qq_plot_vs_uniform.pdf")
    # Input file stems
    input_stem <- paste0(cht_output_dir, "cht_results_",parameter_string,"_num_pc_",pc_num,"_time_")
    null_stem <- paste0(cht_output_dir,"cht_perm_results_",parameter_string,"_num_pc_",pc_num,"_time_")
    #qq_plot_vs_uniform(input_stem,null_stem,output_file)
}


###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to permuted
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)
for (pc_num in 0:5) {
    # Output file
    output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_qq_plot_vs_permuted.pdf")
    # Input file stems
    input_stem <- paste0(cht_output_dir, "cht_results_",parameter_string,"_num_pc_",pc_num,"_time_")
    null_stem <- paste0(cht_output_dir,"cht_perm_results_",parameter_string,"_num_pc_",pc_num,"_time_")
    # qq_plot_vs_permuted(input_stem,null_stem,output_file)
}





###############################################
# Boxplot showing pvalues found in our data for only the eqtls in a specific data set
# 1 plot per data set
# 1 plot for pc
for (pc_num in 0:5) {
    #eqtl_comparison_to_background_shell(pc_num, cht_visualization_dir, parameter_string, cht_enrichment_dir, eqtl_data_set_file)
}



###############################################
# Boxplot showing pvalues found in our data based on eqtls found in GTEx v7:
## 1. HLV
## 2. skin_not_sun_exposed
## 3. stomach
## 4. Breast
# 1 plot for pc
for (pc_num in 0:5) {
    tissue_comparison_plot_file <- paste0(cht_visualization_dir,parameter_string,"_num_pc_",pc_num, "_gtex_tissue_beta_filter_comparison.png")
    #gtex_tissue_comparison_bar_plot(tissue_comparison_plot_file, pc_num, parameter_string, cht_enrichment_dir)
}




# Clustering stuff (not quite done)



k <- 9 # Number of cluster centers for kmeans

pc_num <- 0
    alpha_file <- paste0(cht_output_dir, "best_variant_per_egene_", parameter_string, "_num_pc_", pc_num, "_alpha.txt")
    beta_file <- paste0(cht_output_dir, "best_variant_per_egene_", parameter_string, "_num_pc_", pc_num, "_beta.txt")
    pvalue_file <- paste0(cht_output_dir, "best_variant_per_egene_", parameter_string, "_num_pc_", pc_num, "_pvalues.txt")


    # Load in summary statistic data
    # Each data structure is of dimension (Num_samples) X (Num_time steps)
    alpha <- read_in_summary_statistic_data(alpha_file)
    beta <- read_in_summary_statistic_data(beta_file)
    pvalue <- read_in_summary_statistic_data(pvalue_file)
    # Compute allelic fraction (or p as it is referred to in the WASP paper)
    allelic_fraction <- alpha/(alpha + beta)


    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_kmeans_clustering_", k, "_allelic_fraction.pdf")
    #kmeans_cluster_of_summary_statistics(allelic_fraction, output_file, pc_num, k, "allelic fraction")

    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_kmeans_clustering_", k, "_pvalues.pdf")
    #kmeans_cluster_of_summary_statistics(-log10(pvalue + .00000001), output_file, pc_num, k, "-log10(pvalue)")


