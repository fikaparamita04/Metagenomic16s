# Use the example data in ggpicrust2 package

# Analyze the KO data
data(ko_abundance)
data(metadata)
results_data_input <- ggpicrust2(data = ko_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")
results_data_input[[1]]$plot
KO_results <- results_data_input[[1]]$results

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
ko_daa_results_df <- pathway_daa(abundance = ko_abundance %>% column_to_rownames("#NAME"), metadata = metadata, group = "Environment", daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
ko_daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = ko_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
ko_feature_with_p_0.05 <- ko_daa_results_df %>% filter(p_adjust < 0.05) %>% top_n(10, p_adjust)
pathway_heatmap(abundance = ko_abundance %>% right_join(ko_feature_with_p_0.05 %>% select(all_of(c("feature"))),by = c("#NAME" = "feature")) %>% column_to_rownames("#NAME"), metadata = metadata, group = "Environment")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = ko_abundance %>% column_to_rownames("#NAME"), metadata = metadata, group = "Environment")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
ko_daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = ko_abundance %>% column_to_rownames("#NAME"), metadata = metadata, group = "Environment", daa_method = method)
})

# Compare results across different methods
ko_comparison_results <- compare_daa_results(daa_results_list = ko_daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))


#############################################################################

# Analyze the EC or MetaCyc pathway
data(metacyc_abundance)
results_file_input <- ggpicrust2(data = metacyc_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "MetaCyc",
                                 daa_method = "LinDA",
                                 ko_to_kegg = FALSE,
                                 order = "group",
                                 p_values_bar = TRUE,
                                 x_lab = "description")
results_file_input[[1]]$plot
metacyc_results <- results_file_input[[1]]$results


# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05) %>% top_n(10, p_adjust)
pathway_heatmap(abundance = metacyc_abundance %>% right_join(metacyc_feature_with_p_0.05 %>% select(all_of(c("feature"))),by = c("pathway" = "feature")) %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))
