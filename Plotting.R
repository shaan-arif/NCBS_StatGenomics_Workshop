library(tidyverse)
library(readxl)
if (!requireNamespace("treemapify", quietly = TRUE)) install.packages("treemapify")
library(treemapify)

# ---- file path ----
file_path <- "~/SGW project/DNA_Signature_Data_PCAWG_and_HMF.xlsx"

# ---- Read sheet 'a' ----
all_sheets <- excel_sheets(file_path)
sheet_a <- all_sheets[tolower(all_sheets) == "a"]
if (length(sheet_a) == 0) stop("Sheet 'a' not found.")
sheet_a <- sheet_a[1]

df <- read_excel(file_path, sheet = sheet_a)

# ---- Clean column name " Tumor_location " ----
# Find exact column
old_colname <- colnames(df)[trimws(colnames(df)) == "Tumor_location"]

if (length(old_colname) == 0) {
  stop("Column named ' Tumor_location ' (with spaces) not found. Check column names with colnames(df).")
}

df <- df %>%
  rename(Tumor_location = !!old_colname)  # rename to clean name

# ---- Detect sample ID column ----
possible_ids <- colnames(df)[str_detect(tolower(colnames(df)), "sample|id")]
chosen_id <- NULL
for (p in c("Sample_ID","sample_id","sampleID","sample","ID","id")) {
  if (p %in% colnames(df)) { chosen_id <- p; break }
}
if (is.null(chosen_id) && length(possible_ids) > 0) chosen_id <- possible_ids[1]
if (is.null(chosen_id)) {
  message("No Sample_ID detected; cohort set to Unknown.")
  df <- df %>% mutate(Sample_ID = NA_character_)
} else {
  df <- df %>% rename(Sample_ID = !!sym(chosen_id))
}

# ---- Cohort classification ----
df <- df %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    cohort = case_when(
      str_detect(Sample_ID, regex("^\\s*HMF", ignore_case = TRUE)) ~ "HMF",
      str_detect(Sample_ID, regex("^\\s*SA",  ignore_case = TRUE)) ~ "PCAWG",
      TRUE ~ "Other"
    ),
    Tumor_location = str_squish(as.character(Tumor_location)),
    Tumor_location = if_else(is.na(Tumor_location) | Tumor_location == "", 
                             "Unknown", Tumor_location)
  )

# ---- Aggregate counts ----
agg <- df %>%
  count(Tumor_location, cohort, name = "n") %>%
  group_by(Tumor_location) %>%
  mutate(total = sum(n)) %>%
  ungroup()

# Choose top 20 Tumor_locations
top_n <- 20
top_locs <- agg %>%
  distinct(Tumor_location, total) %>%
  arrange(desc(total)) %>%
  slice_head(n = top_n) %>%
  pull(Tumor_location)

agg_top <- agg %>% filter(Tumor_location %in% top_locs)

# =============================
#  Plot 1: TREEMAP
# =============================
treemap_df <- agg_top %>%
  distinct(Tumor_location, total) %>%
  arrange(desc(total))

p_treemap <- ggplot(treemap_df, aes(area = total, label = Tumor_location, fill = total)) +
  geom_treemap() +
  geom_treemap_text(
    aes(label = paste0(Tumor_location, "\n(", total, ")")),
    colour = "white", place = "centre", grow = TRUE, reflow = TRUE
  ) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(title = paste("Treemap — Top", top_n, "Tumor Locations"),
       fill = "Total Samples") +
  theme_minimal()

print(p_treemap)


# =============================
#  Plot 2: FACETED DOT + LOLLIPOP
# =============================
agg_top <- agg_top %>%
  mutate(Tumor_location = fct_reorder(Tumor_location, total))

p_facet <- ggplot(agg_top, aes(x = n, y = Tumor_location, color = cohort)) +
  geom_segment(aes(x = 0, xend = n, y = Tumor_location, yend = Tumor_location),
               size = 0.6, alpha = 0.6) +
  geom_point(aes(size = n)) +
  scale_size(range = c(2, 8)) +
  facet_wrap(~ cohort, scales = "free_x") +
  labs(title = "Counts by Tumor Location — Split by Cohort",
       x = "Samples", y = "Tumor Location") +
  theme_minimal()

print(p_facet)


# =============================
# Plot 3: STACKED PROPORTION PLOT
# =============================
stacked <- agg_top %>%
  group_by(Tumor_location) %>%
  mutate(pct = n / total) %>%
  ungroup()

p_stack <- ggplot(stacked, aes(x = Tumor_location, y = pct, fill = cohort)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cohort Composition for Top Tumor Locations",
       x = "Tumor Location", y = "Percentage") +
  theme_minimal()

print(p_stack)


# ===============================
# CONTINUATION: Sankey, Circular bar, Chord, Plotly treemap, Grouped bar
# Assumes `df` exists with columns:
#   - Tumor_location (clean, character)
#   - cohort (values e.g. "HMF", "SA", "Other")
# ===============================
library(tidyverse)
# networkD3 for Sankey
if (!requireNamespace("networkD3", quietly = TRUE)) install.packages("networkD3")
library(networkD3)
# circlize for chord
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
library(circlize)
# plotly for interactive treemap
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
library(plotly)
# scales for formatting
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
library(scales)

# Safety check: df must exist
if (!exists("df")) stop("Data frame `df` not found. Run the previous code that reads & cleans the data (creates df with Tumor_location and cohort).")

# Choose top tumor locations for most plots (to keep visuals readable)
top_n <- 30
top_locations <- df %>%
  count(Tumor_location, name = "total") %>%
  arrange(desc(total)) %>%
  slice_head(n = top_n) %>%
  pull(Tumor_location)

# Filter to top locations where appropriate
df_top <- df %>% filter(Tumor_location %in% top_locations)

# -----------------------
# 1) SANKEY (networkD3)
# -----------------------
# Prepare links: counts from Tumor_location -> cohort
links_df <- df_top %>%
  count(Tumor_location, cohort, name = "value") %>%
  arrange(desc(value))

# Create nodes: unique tumor locations + unique cohorts
nodes <- tibble(name = c(unique(links_df$Tumor_location), unique(links_df$cohort)))

# helper to get node index
get_index <- function(x) match(x, nodes$name) - 1  # networkD3 uses zero-based indices

links <- links_df %>%
  mutate(
    source = get_index(Tumor_location),
    target = get_index(cohort)
  ) %>%
  select(source, target, value)

# Plot Sankey
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes,
                             Source = "source", Target = "target", Value = "value",
                             NodeID = "name", fontSize = 12, nodeWidth = 30, sinksRight = FALSE)
# Print interactive htmlwidget
sankey_plot

# Optionally save as html (requires htmlwidgets)
if (requireNamespace("htmlwidgets", quietly = TRUE)) {
  htmlwidgets::saveWidget(sankey_plot, "sankey_tumor_cohort.html", selfcontained = TRUE)
  message("Sankey saved: sankey_tumor_cohort.html")
}

# -----------------------
# 2) CIRCULAR BAR PLOT (polar)
# -----------------------
# Aggregate totals by Tumor_location for top_n
circ_df <- df_top %>%
  count(Tumor_location, name = "total") %>%
  arrange(desc(total)) %>%
  mutate(Tumor_location = fct_reorder(Tumor_location, total))

# Create circular bar (using ggplot2)
p_circular <- ggplot(circ_df, aes(x = Tumor_location, y = total)) +
  geom_col(fill = "steelblue", width = 0.8, show.legend = FALSE) +
  coord_polar(start = 0) +
  labs(title = paste0("Circular bar plot — top ", top_n, " Tumor locations"),
       x = "", y = "Number of samples") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

print(p_circular)
ggsave("circular_bar_top_tumors.png", p_circular, width = 10, height = 8, dpi = 300)

# -----------------------
# 3) CHORD DIAGRAM (circlize)
# -----------------------
# Prepare a dataframe: from = Tumor_location, to = cohort, value = count
chord_df <- df_top %>%
  count(Tumor_location, cohort, name = "value") %>%
  filter(value > 0)

# circlize wants a matrix or a data frame with columns from, to, value
# Choose a color palette for cohorts vs tumor locations
# Build a vector of colors matching unique categories
tumors_unique <- unique(chord_df$Tumor_location)
cohorts_unique <- unique(chord_df$cohort)
all_nodes <- c(tumors_unique, cohorts_unique)

# Color mapping: tumors get a palette, cohorts get distinct colors
tumor_cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(tumors_unique))
cohort_cols <- c("HMF" = "#D55E00", "SA" = "#0072B2", "Other" = "#009E73")
node_colors <- c(tumor_cols, cohort_cols[cohorts_unique])

# Plot chord diagram
# Use a temporary plot device to avoid overwriting graphics
circos.clear()
# Use directional chord (from tumor -> cohort)
# chordDiagram can accept data frame
circos.par(gap.degree = 2, start.degree = 90)
p_pdf <- tryCatch({
  chordDiagram(x = chord_df, grid.col = node_colors, transparency = 0.35,
               annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05))
}, error = function(e) { message("Error in chordDiagram: ", e$message); NULL })
# Add labels (works automatically)
title("Chord diagram — Tumor_location to Cohort (top locations)")

# Save chord as PNG by opening a device first (optional)
png("chord_tumor_cohort.png", width = 1600, height = 1200, res = 150)
circos.clear()
circos.par(gap.degree = 2, start.degree = 90)
chordDiagram(x = chord_df, grid.col = node_colors, transparency = 0.35,
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05))
title("Chord diagram — Tumor_location to Cohort (top locations)")
dev.off()
message("Chord diagram saved: chord_tumor_cohort.png")
circos.clear()

# -----------------------
# 4) INTERACTIVE PLOTLY TREEMAP
# -----------------------
treemap_df <- df_top %>%
  count(Tumor_location, name = "total") %>%
  arrange(desc(total)) %>%
  slice_head(n = top_n)

# plotly treemap
p_plotly_tm <- plot_ly(
  treemap_df,
  type = "treemap",
  values = ~total,
  labels = ~Tumor_location,
  textinfo = "label+value",
  hoverinfo = "label+value",
  marker = list(colors = treemap_df$total, colorscale = "Viridis")
)
p_plotly_tm <- p_plotly_tm %>% layout(title = paste("Interactive treemap — top", top_n, "Tumor locations"))
p_plotly_tm

# Optionally save as an html widget
htmlwidgets::saveWidget(p_plotly_tm, "interactive_treemap_top_tumors.html", selfcontained = TRUE)
message("Interactive treemap saved: interactive_treemap_top_tumors.html")

# -----------------------
# 5) GROUPED / DODGED BAR PLOT (Tumor types vs HMF & SA)
# -----------------------
# Prepare counts for HMF and SA specifically (and Other optionally)
group_df <- df %>%
  filter(Tumor_location %in% top_locations) %>%
  filter(cohort %in% c("HMF", "SA")) %>%   # include only HMF and SA for clear comparison
  count(Tumor_location, cohort, name = "n") %>%
  complete(Tumor_location = top_locations, cohort = c("HMF", "SA"), fill = list(n = 0)) %>%
  mutate(Tumor_location = fct_reorder(Tumor_location, -n, .fun = sum))

# Plot grouped (dodged) bar chart so HMF and SA are side-by-side
p_grouped <- ggplot(group_df, aes(x = Tumor_location, y = n, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("HMF" = "#D55E00", "SA" = "#0072B2")) +
  labs(title = "Comparison of sample counts by Tumor_location: HMF vs SA (top locations)",
       x = "Tumor location", y = "Number of samples", fill = "Cohort") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))

print(p_grouped)
ggsave("grouped_bar_hmf_sa_top_tumors.png", p_grouped, width = 12, height = 10, dpi = 300)

# If you want One plot that shows HMF and SA stacked as proportions:
prop_df <- df %>%
  filter(Tumor_location %in% top_locations) %>%
  count(Tumor_location, cohort, name = "n") %>%
  group_by(Tumor_location) %>%
  mutate(pct = n / sum(n)) %>% ungroup() %>%
  mutate(Tumor_location = fct_reorder(Tumor_location, -pct, .fun = sum))

p_stacked_prop <- ggplot(prop_df %>% filter(cohort %in% c("HMF", "SA")), aes(x = Tumor_location, y = pct, fill = cohort)) +
  geom_col(position = "fill") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("HMF" = "#D55E00", "SA" = "#0072B2")) +
  labs(title = "Proportion of HMF vs SA within each Tumor location (top locations)",
       x = "Tumor location", y = "Proportion", fill = "Cohort") +
  theme_minimal()

print(p_stacked_prop)
ggsave("stacked_prop_hmf_sa_top_tumors.png", p_stacked_prop, width = 12, height = 10, dpi = 300)

# -----------------------

# -----------------------
pkgs <- c("tidyverse","readxl","uwot","Rtsne","pheatmap","cluster","factoextra","mclust","scales")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse)
library(readxl)
library(uwot)
library(Rtsne)
library(pheatmap)
library(cluster)      # silhouette
library(factoextra)   # visualization helpers (optional)
library(mclust)       # adjustedRandIndex
library(scales)

set.seed(42)

# ---------- 1) Load & prepare data (adjust path if needed) ----------
file_path <- "~/SGW project/DNA_Signature_Data_PCAWG_and_HMF.xlsx"
sheets_all <- excel_sheets(file_path)
sheets_to_load <- sheets_all[!tolower(sheets_all) %in% tolower("C")]
data_list <- map(set_names(sheets_to_load), ~ read_excel(file_path, sheet = .x))
combined_df <- bind_rows(data_list, .id = "source_sheet")

# Clean the " Tumor_location " if present and create Tumor_location column
col_trimmed <- tibble(original = colnames(combined_df),
                      trimmed = trimws(colnames(combined_df)))
if ("Tumor_location" %in% col_trimmed$trimmed) {
  exact_old <- col_trimmed$original[col_trimmed$trimmed == "Tumor_location"][1]
  combined_df <- combined_df %>% rename(Tumor_location = !!sym(exact_old))
}

# Ensure Sample_ID and cohort exist (HMF and PCAWG)
possible_id_cols <- colnames(combined_df)[str_detect(tolower(colnames(combined_df)), "sample|^id$|^sid$")]
chosen_id <- possible_id_cols[1] %||% NA_character_
if (!is.na(chosen_id)) combined_df <- combined_df %>% rename(Sample_ID = !!sym(chosen_id)) else combined_df$Sample_ID <- NA_character_
combined_df <- combined_df %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    cohort = case_when(
      str_detect(Sample_ID, regex("^\\s*HMF", ignore_case = TRUE)) ~ "HMF",
      str_detect(Sample_ID, regex("^\\s*SA",  ignore_case = TRUE)) ~ "PCAWG",
      TRUE ~ "Other"
    ),
    Tumor_location = ifelse(is.na(Tumor_location) | trimws(Tumor_location) == "", "Unknown", as.character(Tumor_location))
  )

# ---------- 2) Select numeric mutation-feature columns ----------
# Exclude common metadata columns (case-insensitive)
meta_patterns <- c("sample","id","source_sheet","tumor_location","tumor","cohort","cancer","type","site")
meta_cols <- colnames(combined_df)[map_lgl(colnames(combined_df), ~ any(str_detect(tolower(.x), meta_patterns)))]
# numeric candidates: try converting to numeric, keep columns that become numeric (and not metadata)
candidate_df <- combined_df %>% select(-any_of(meta_cols))
# Determine numeric columns robustly
numeric_cols <- candidate_df %>% select(where(~ all(is.na(.x)) || is.numeric(.x) || suppressWarnings(!any(is.na(as.numeric(as.character(.x))))))) %>% colnames()
# A safer filter: keep cols that have at least some non-NA numeric values and not all unique strings
numeric_cols <- numeric_cols[ sapply(numeric_cols, function(cn){
  x <- candidate_df[[cn]]
  suppressWarnings(sum(!is.na(as.numeric(as.character(x)))) ) > max(5, 0.01 * nrow(candidate_df))
}) ]

if (length(numeric_cols) < 2) stop("Not enough numeric feature columns detected. Inspect your data and adjust numeric_cols selection.")

message("Using ", length(numeric_cols), " numeric features for unsupervised analysis.")

# Build numeric matrix (scaled)
X <- combined_df %>% select(all_of(numeric_cols)) %>% mutate(across(everything(), ~ as.numeric(as.character(.x))))
# Impute simple: column median for NA (you can swap for better methods)
X_imp <- X %>% mutate(across(everything(), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))
X_scaled <- scale(X_imp) %>% as.matrix()

# Keep metadata for plotting (Tumor_location, cohort, Sample_ID)
meta <- combined_df %>% select(Sample_ID, Tumor_location, cohort) %>% mutate(Tumor_location = as.factor(Tumor_location), cohort = as.factor(cohort))

# For plotting readability, optionally restrict to top tumor locations (by count)
top_n_loc <- 25
top_locs <- meta %>% count(Tumor_location, sort = TRUE) %>% slice_head(n = top_n_loc) %>% pull(Tumor_location)
keep_idx <- which(meta$Tumor_location %in% top_locs)
# You can comment out the next line to use all samples
#X_plot <- X_scaled[keep_idx, , drop = FALSE]; meta_plot <- meta[keep_idx, ]
X_plot <- X_scaled; meta_plot <- meta  # using all samples by default

# ---------- 3) PCA ----------
pca_res <- prcomp(X_plot, center = TRUE, scale. = TRUE)
pca_df <- as_tibble(pca_res$x[,1:6]) %>% bind_cols(meta_plot)

# Percent variance
var_expl <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_var <- percent(var_expl[1])
pc2_var <- percent(var_expl[2])

# PCA plot colored by Tumor_location (top N) and by cohort
library(ggplot2)
p_pca_loc <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tumor_location)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "PCA: PC1 vs PC2 (colored by Tumor_location)", subtitle = paste0("Var PC1: ", pc1_var, ", PC2: ", pc2_var),
       x = paste0("PC1 (", round(100 * var_expl[1],1), "%)"), y = paste0("PC2 (", round(100 * var_expl[2],1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")

p_pca_cohort <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cohort, shape = cohort)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA: colored by cohort") +
  theme_minimal() +
  scale_color_manual(values = c("HMF" = "#D55E00","PCAWG" = "#0072B2","Other" = "grey"))

print(p_pca_loc)
print(p_pca_cohort)
ggsave("PCA_PC1_PC2_Tumor_location.png", p_pca_loc, width = 10, height = 8, dpi = 300)
ggsave("PCA_PC1_PC2_cohort.png", p_pca_cohort, width = 8, height = 6, dpi = 300)

# ---------- 4) UMAP ----------
set.seed(42)
umap_res <- uwot::umap(X_plot, n_components = 2, n_neighbors = 15, metric = "cosine")
umap_df <- as_tibble(umap_res) %>% setNames(c("UMAP1","UMAP2")) %>% bind_cols(meta_plot)

p_umap_loc <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Tumor_location)) +
  geom_point(alpha = 0.7, size = 1.8) +
  labs(title = "UMAP projection (colored by Tumor_location)") +
  theme_minimal()

p_umap_cohort <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cohort)) +
  geom_point(alpha = 0.8, size = 1.8) +
  labs(title = "UMAP: colored by cohort") +
  theme_minimal() +
  scale_color_manual(values = c("HMF" = "#D55E00","PCAWG" = "#0072B2","Other" = "grey"))

print(p_umap_loc)
print(p_umap_cohort)
ggsave("UMAP_Tumor_location.png", p_umap_loc, width = 10, height = 8, dpi = 300)

# ---------- 5) t-SNE ----------
set.seed(42)
tsne_res <- Rtsne(X_plot, perplexity = 30, pca = TRUE, check_duplicates = FALSE) # adjust perplexity by sample size
tsne_df <- as_tibble(tsne_res$Y) %>% setNames(c("tSNE1","tSNE2")) %>% bind_cols(meta_plot)

p_tsne_loc <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tumor_location)) +
  geom_point(alpha = 0.7, size = 1.8) +
  labs(title = "t-SNE (colored by Tumor_location)") +
  theme_minimal()

print(p_tsne_loc)
ggsave("tSNE_Tumor_location.png", p_tsne_loc, width = 10, height = 8, dpi = 300)

# ---------- 6) Hierarchical clustering + heatmap of top features ----------
# For heatmap, showing all features is heavy; choose features with highest variance
feat_var <- apply(X_scaled, 2, var, na.rm = TRUE)
top_feat_n <- 100   # adjust if needed
top_feats <- names(sort(feat_var, decreasing = TRUE))[1:min(top_feat_n, length(feat_var))]

# build matrix for heatmap (samples x top_feats), cluster rows & columns
mat_heat <- X_scaled[, top_feats]
rownames(mat_heat) <- meta$Sample_ID %||% seq_len(nrow(mat_heat))

# annotation for rows
anno_row <- data.frame(Tumor_location = meta$Tumor_location, cohort = meta$cohort)
rownames(anno_row) <- rownames(mat_heat)

# pheatmap
pheatmap(mat_heat, show_rownames = FALSE, show_colnames = TRUE, annotation_row = anno_row,
         clustering_method = "ward.D2", main = "Heatmap (top variable features)")

# Save heatmap as PNG
png("heatmap_top_features.png", width = 1400, height = 1000, res = 150)
pheatmap(mat_heat, show_rownames = FALSE, show_colnames = TRUE, annotation_row = anno_row,
         clustering_method = "ward.D2", main = "Heatmap (top variable features)")
dev.off()

# ---------- 7) k-means clustering: choose k by silhouette and WSS ----------
# limit to reasonable sample size for speed (you can use X_plot which may be filtered)
max_k <- 10
wss <- sapply(1:max_k, function(k) {
  kmeans(X_plot, centers = k, nstart = 10, iter.max = 100)$tot.withinss
})

# silhouette average width per k (compute for k >= 2)
avg_sil <- sapply(2:max_k, function(k) {
  km <- kmeans(X_plot, centers = k, nstart = 10)
  ss <- silhouette(km$cluster, dist(X_plot))
  mean(ss[, 3])
})

# Plot WSS and silhouette
df_k <- tibble(k = 1:max_k, wss = wss)
p_wss <- ggplot(df_k, aes(x = k, y = wss)) + geom_line() + geom_point() + labs(title = "Elbow: WSS by k")
print(p_wss)
ggsave("k_elbow_wss.png", p_wss, width = 7, height = 5, dpi = 200)

df_s <- tibble(k = 2:max_k, avg_sil = avg_sil)
p_sil <- ggplot(df_s, aes(x = k, y = avg_sil)) + geom_line() + geom_point() + labs(title = "Average silhouette by k")
print(p_sil)
ggsave("k_silhouette.png", p_sil, width = 7, height = 5, dpi = 200)

# Choose k with max silhouette (or manual)
best_k <- which.max(c(NA, avg_sil))  # adds NA for k=1 index alignment; returns index
if (best_k < 2) best_k <- 2
message("Selected k = ", best_k)

# Run kmeans with best_k and project clusters onto PCA
km_final <- kmeans(X_plot, centers = best_k, nstart = 20)
pca_df$cluster_kmeans <- factor(km_final$cluster[rownames(pca_df) %||% seq_len(nrow(pca_df))])

p_km_on_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_kmeans)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = paste0("K-means (k=", best_k, ") clusters on PCA")) +
  theme_minimal()
print(p_km_on_pca)
ggsave("kmeans_on_pca.png", p_km_on_pca, width = 8, height = 6, dpi = 200)

# ---------- 8) Compare clustering to true Tumor_location: confusion table & ARI ----------
# Use only samples with non-Unknown Tumor_location
valid_idx <- which(!is.na(meta_plot$Tumor_location) & meta_plot$Tumor_location != "Unknown")
if (length(valid_idx) > 1) {
  ari_val <- adjustedRandIndex(km_final$cluster[valid_idx], as.numeric(as.factor(meta_plot$Tumor_location[valid_idx])))
  message("Adjusted Rand Index (kmeans vs Tumor_location): ", round(ari_val, 4))
  
  # show confusion for top tumor types
  true_fac <- factor(meta_plot$Tumor_location[valid_idx])
  tab <- table(cluster = km_final$cluster[valid_idx], truth = true_fac)
  print(tab)
} else {
  message("Not enough labeled Tumor_location data to compute ARI.")
}

pkgs <- c("tidyverse", "readxl", "ggpubr", "patchwork", "viridis", "pheatmap", "scales")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(viridis)
library(pheatmap)
library(scales)

set.seed(42)

# ---------------------------
# 1) Load + basic cleaning (exclude sheet C)
# ---------------------------
file_path <- "~/SGW project/DNA_Signature_Data_PCAWG_and_HMF.xlsx"
sheets_all <- excel_sheets(file_path)
sheets_to_load <- sheets_all[!tolower(sheets_all) %in% tolower("C")]
data_list <- map(set_names(sheets_to_load), ~ read_excel(file_path, sheet = .x))
combined_df <- bind_rows(data_list, .id = "source_sheet")

# Clean " Tumor_location " variations -> Tumor_location
col_trimmed <- tibble(original = colnames(combined_df), trimmed = trimws(colnames(combined_df)))
if ("Tumor_location" %in% col_trimmed$trimmed) {
  exact_old <- col_trimmed$original[col_trimmed$trimmed == "Tumor_location"][1]
  combined_df <- combined_df %>% rename(Tumor_location = !!sym(exact_old))
}

# Detect sample ID column and create Sample_ID (robust)
possible_id_cols <- colnames(combined_df)[str_detect(tolower(colnames(combined_df)), "sample|sid|sampleid|sample_id|id$")]
chosen_id <- NULL
for (p in c("Sample_ID","sample_id","SampleID","sampleid","sample","Sample","ID","id")) {
  if (p %in% colnames(combined_df)) { chosen_id <- p; break }
}
if (!is.null(chosen_id)) {
  combined_df <- combined_df %>% rename(Sample_ID = !!sym(chosen_id))
} else {
  combined_df <- combined_df %>% mutate(Sample_ID = NA_character_)
}

# Create cohort: HMF or PCAWG (SA -> PCAWG); else Other
combined_df <- combined_df %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    cohort = case_when(
      str_detect(Sample_ID, regex("^\\s*HMF", ignore_case = TRUE)) ~ "HMF",
      str_detect(Sample_ID, regex("^\\s*SA",  ignore_case = TRUE)) ~ "PCAWG",
      TRUE ~ "Other"
    ),
    Tumor_location = if_else(is.na(Tumor_location) | trimws(as.character(Tumor_location)) == "", "Unknown", as.character(Tumor_location))
  )

# ---------------------------
# 2) Identify numeric mutation-feature columns (exclude metadata)
# ---------------------------
meta_patterns <- c("sample","id","source_sheet","tumor_location","tumor","cohort","cancer","type","site","age","sex","gender")
meta_cols <- colnames(combined_df)[map_lgl(colnames(combined_df), ~ any(str_detect(tolower(.x), meta_patterns)))]

# Candidate columns excluding metadata
candidate_cols <- setdiff(colnames(combined_df), meta_cols)

# Convert candidate columns to numeric where possible and keep those with enough numeric entries
numeric_cols <- candidate_cols[sapply(candidate_cols, function(cn){
  x <- combined_df[[cn]]
  # can coerce to numeric? count non-NA numeric conversions
  ok <- suppressWarnings(!all(is.na(as.numeric(as.character(x)))))
  # also require at least some numeric observations (>= 5 or >=1% of samples)
  num_ok <- suppressWarnings(sum(!is.na(as.numeric(as.character(x))))) 
  ok && (num_ok >= max(5, 0.01 * nrow(combined_df)))
})]

# Build numeric feature dataframe (coerce to numeric)
num_df <- combined_df %>%
  select(all_of(numeric_cols)) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.x))))

# If no numeric features found, stop
if (ncol(num_df) < 2) stop("No numeric mutation-feature columns were detected. Inspect column names and metadata patterns.")

# ---------------------------
# 3) Which features have highest variance?
# ---------------------------
feature_var <- apply(num_df, 2, var, na.rm = TRUE) %>% sort(decreasing = TRUE)
top_n_var <- 20
top_var_feats <- names(feature_var)[1:min(top_n_var, length(feature_var))]
var_df <- tibble(feature = names(feature_var), variance = as.numeric(feature_var)) %>% arrange(desc(variance))

# Save CSV
write_csv(var_df, "feature_variances.csv")

# Plot 1: Barplot of top variance features
p_varbar <- var_df %>%
  slice_head(n = top_n_var) %>%
  mutate(feature = fct_reorder(feature, variance)) %>%
  ggplot(aes(x = feature, y = variance, fill = variance)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "magma") +
  labs(title = paste("Top", top_n_var, "mutational features by variance"),
       x = "Feature", y = "Variance") +
  theme_minimal()

ggsave("top_variance_features_bar.png", p_varbar, width = 10, height = 7, dpi = 300)
print(p_varbar)

# Plot 2: Violin + boxplots for top features to see distribution across all samples
# Pick top 6 high-variance features to visualize
top_show <- top_var_feats[1:min(6, length(top_var_feats))]  

# Build long format data directly from combined_df so row counts match
long_top <- combined_df %>%
  select(Sample_ID, Tumor_location, cohort, all_of(top_show)) %>%
  pivot_longer(cols = all_of(top_show), names_to = "feature", values_to = "value")
# adds metadata

p_violin_top <- ggplot(long_top, aes(x = feature, y = value)) +
  geom_violin(fill = "lightgray") +
  geom_boxplot(width = 0.12, outlier.size = 0.7) +
  coord_flip() +
  labs(title = "Distribution of top mutational features (violin + box)", y = "Value", x = "") +
  theme_minimal()

ggsave("top_features_violin.png", p_violin_top, width = 10, height = 6, dpi = 300)
print(p_violin_top)

# ---------------------------
# 4) Mutational burden: compute and compare HMF vs PCAWG
# ---------------------------
# Define mutational burden as sum of numeric mutation features per sample (simple proxy)
mutational_burden <- rowSums(num_df, na.rm = TRUE)
combined_df <- combined_df %>% mutate(mutational_burden = mutational_burden)

# Keep only samples with cohort HMF or PCAWG
df_burden <- combined_df %>% filter(cohort %in% c("HMF", "PCAWG")) %>% 
  mutate(cohort = factor(cohort, levels = c("PCAWG", "HMF")))

# Basic summary
burden_summary <- df_burden %>% group_by(cohort) %>%
  summarise(n = n(), median = median(mutational_burden, na.rm = TRUE), mean = mean(mutational_burden, na.rm = TRUE), sd = sd(mutational_burden, na.rm = TRUE))
write_csv(burden_summary, "mutational_burden_summary_by_cohort.csv")
print(burden_summary)

# Wilcoxon test (non-parametric)
wil <- wilcox.test(mutational_burden ~ cohort, data = df_burden, exact = FALSE)
# also compute effect size (rank-biserial) using common formula: r = Z / sqrt(n)
# but wilcox.test in R doesn't return Z. We'll compute Cliff's delta as effect size:
if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")
library(effsize)
cliff <- cliff.delta(df_burden$mutational_burden[df_burden$cohort=="HMF"], df_burden$mutational_burden[df_burden$cohort=="PCAWG"])

# Plot: violin + boxplot + jitter + p-value annotation
p_burden <- ggplot(df_burden, aes(x = cohort, y = mutational_burden, fill = cohort)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  scale_y_continuous(trans = "log10", labels = comma) +   # log scale often helps visualizing counts
  scale_fill_manual(values = c("PCAWG" = "#0072B2", "HMF" = "#D55E00")) +
  labs(title = "Mutational burden (sum of numeric features): PCAWG vs HMF",
       subtitle = paste0("Wilcoxon p = ", signif(wil$p.value, 3), "; Cliff's delta (HMF vs PCAWG) = ", round(cliff$estimate,3)),
       x = "", y = "Mutational burden (log10 scale)") +
  theme_minimal()

ggsave("mutational_burden_cohort_violin.png", p_burden, width = 8, height = 6, dpi = 300)
print(p_burden)

# Optional: density by cohort
p_burden_den <- ggplot(df_burden, aes(x = mutational_burden, color = cohort)) +
  geom_density(adjust = 1.5) +
  scale_x_continuous(trans = "log10") +
  labs(title = "Density of mutational burden (log scale)", x = "Mutational burden (log10)") +
  theme_minimal()
ggsave("mutational_burden_density.png", p_burden_den, width = 8, height = 5, dpi = 300)
print(p_burden_den)

# Save test results
sink("mutational_burden_test_results.txt")
print(wil)
print(cliff)
sink()

# ---------------------------
# 5) Correlations with clinical variables (age numeric, sex categorical)
# ---------------------------
# Detect age column (common names) and sex/gender
age_col_candidates <- colnames(combined_df)[str_detect(tolower(colnames(combined_df)), "age|patient_age|age_at")]
sex_col_candidates <- colnames(combined_df)[str_detect(tolower(colnames(combined_df)), "sex|gender")]

age_col <- if (length(age_col_candidates) > 0) age_col_candidates[1] else NA_character_
sex_col <- if (length(sex_col_candidates) > 0) sex_col_candidates[1] else NA_character_

# Convert age to numeric if present
if (!is.na(age_col)) {
  combined_df <- combined_df %>% mutate(age_numeric = as.numeric(as.character(.data[[age_col]])))
  age_present <- TRUE
} else age_present <- FALSE

# Clean sex if present
if (!is.na(sex_col)) {
  combined_df <- combined_df %>% mutate(sex_clean = as.character(.data[[sex_col]]) %>% str_to_title() %>% str_trim())
  sex_present <- TRUE
} else sex_present <- FALSE

# 5a) If age present: compute Spearman correlation between features and age
if (age_present) {
  # Use features in num_df
  cor_age <- sapply(colnames(num_df), function(feat) {
    x <- num_df[[feat]]
    y <- combined_df$age_numeric
    # use Spearman; require enough non-NA pairs
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 10) return(NA_real_)
    cor_val <- cor(x[ok], y[ok], method = "spearman")
    return(cor_val)
  })
  cor_age_df <- tibble(feature = names(cor_age), spearman_r = as.numeric(cor_age)) %>%
    drop_na() %>%
    arrange(desc(abs(spearman_r)))
  write_csv(cor_age_df, "feature_age_correlations.csv")
  
  # Heatmap of top features correlated with age (top 20 by absolute rho)
  top_age_feats <- cor_age_df %>% slice_max(order_by = abs(spearman_r), n = min(20, nrow(cor_age_df))) %>% pull(feature)
  mat_age <- num_df %>% select(all_of(top_age_feats)) %>% mutate(across(everything(), ~ scale(.x))) %>% as.matrix()
  rownames(mat_age) <- combined_df$Sample_ID %||% seq_len(nrow(mat_age))
  # order samples by age
  sample_age <- combined_df$age_numeric
  annotation_row <- data.frame(age = sample_age)
  rownames(annotation_row) <- rownames(mat_age)
  pheatmap(mat_age, show_rownames = FALSE, main = "Top features correlated with age (scaled)", annotation_row = annotation_row)
  png("age_feature_heatmap.png", width = 1200, height = 900, res = 150)
  pheatmap(mat_age, show_rownames = FALSE, main = "Top features correlated with age (scaled)", annotation_row = annotation_row)
  dev.off()
  
  # Scatter plot for top single feature vs age
  top_feat_age <- cor_age_df$feature[1]
  scatter_df_age <- tibble(age = combined_df$age_numeric, feature_val = num_df[[top_feat_age]])
  p_scatter_age <- ggplot(scatter_df_age, aes(x = age, y = feature_val)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    labs(title = paste("Top feature correlated with age:", top_feat_age),
         x = "Age", y = top_feat_age) +
    theme_minimal()
  ggsave("top_feature_vs_age_scatter.png", p_scatter_age, width = 7, height = 5, dpi = 300)
  print(p_scatter_age)
} else {
  message("No age column detected. Skipping age correlation analysis.")
}

# 5b) If sex present: compare feature distributions between sexes
if (sex_present) {
  # Keep only two-level sexes ideally (Male/Female)
  combined_df <- combined_df %>% mutate(sex_clean = if_else(sex_clean == "" | is.na(sex_clean), "Unknown", sex_clean))
  sexes_present <- unique(combined_df$sex_clean)
  message("Detected sex/gender categories: ", paste(head(sexes_present,10), collapse = ", "))
  
  # For each feature compute test (t-test or Wilcoxon depending on distribution)
  sex_stats <- tibble(feature = character(), mean_M = numeric(), mean_F = numeric(), p_value = numeric(), method = character())
  for (feat in colnames(num_df)) {
    x <- num_df[[feat]]
    df_tmp <- tibble(x = x, sex = combined_df$sex_clean) %>% filter(!is.na(x) & sex != "Unknown")
    # require at least 10 samples in two groups
    counts <- table(df_tmp$sex)
    if (length(counts) < 2 || any(counts < 5)) {
      next
    }
    # pick two largest groups (commonly Male/Female)
    top_groups <- names(sort(counts, decreasing = TRUE))[1:2]
    df2 <- df_tmp %>% filter(sex %in% top_groups)
    # normality check: small sample -> wilcox; else t-test with var.test?
    # we'll use Wilcoxon (robust)
    wt <- wilcox.test(x ~ sex, data = df2)
    pval <- wt$p.value
    means <- df2 %>% group_by(sex) %>% summarise(m = mean(x, na.rm = TRUE)) %>% pivot_wider(names_from = sex, values_from = m)
    # add a row
    row <- tibble(feature = feat, p_value = pval, method = "wilcox", mean_M = means[[2]][1], mean_F = means[[1]][1]) # approximate positions
    sex_stats <- bind_rows(sex_stats, row)
  }
  sex_stats <- sex_stats %>% arrange(p_value)
  write_csv(sex_stats, "feature_sex_tests.csv")
  
  # plot top 6 features by significance
  top_sex_feats <- sex_stats %>% slice_head(n = 6) %>% pull(feature)
  if (length(top_sex_feats) > 0) {
    long_sex <- num_df %>% select(all_of(top_sex_feats)) %>% pivot_longer(everything(), names_to = "feature", values_to = "value") %>%
      bind_cols(combined_df %>% select(sex_clean))
    p_sex_box <- ggplot(long_sex, aes(x = sex_clean, y = value, fill = sex_clean)) +
      geom_boxplot() +
      facet_wrap(~ feature, scales = "free_y", ncol = 2) +
      theme_minimal() +
      labs(title = "Top features differing by sex (boxplots)")
    ggsave("top_features_by_sex_boxplots.png", p_sex_box, width = 10, height = 8, dpi = 300)
    print(p_sex_box)
  } else {
    message("No features passed sex-test requirements (insufficient counts), or no sex differences detected.")
  }
} else {
  message("No sex/gender column detected. Skipping sex-related tests.")
}

# ---------------------------
# 6) Save useful outputs and final messages
# ---------------------------
# Save top variance list and feature variance CSV already saved earlier
write_csv(var_df, "feature_variances.csv")   # repeated to ensure saved
write_csv(tibble(feature = top_var_feats), "top_variance_features_list.csv")

# Save mutational burden per sample
write_csv(combined_df %>% select(Sample_ID, cohort, Tumor_location, mutational_burden), "mutational_burden_per_sample.csv")

pkgs <- c("tidyverse","readxl","pheatmap","uwot","RColorBrewer","igraph","ggraph","networkD3","plotly","scales","broom","ggrepel")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl); library(pheatmap); library(uwot)
library(RColorBrewer); library(igraph); library(ggraph); library(networkD3)
library(plotly); library(scales); library(broom); library(ggrepel)

set.seed(42)

# ---- user params ----
file_path <- "~/SGW project/DNA_Signature_Data_PCAWG_and_HMF.xlsx"
exclude_sheet <- "C"
top_n <- 25            # used for plots showing top tumor types/signatures
sig_corr_threshold <- 0.6  # absolute correlation threshold for signature network
cluster_k <- 6         # a default k for cluster visualizations (can try other values)
# ----------------------------------------------------------------------

# ---- 1) Read & basic cleaning ----
sheets_all <- excel_sheets(file_path)
sheets_to_load <- sheets_all[!tolower(sheets_all) %in% tolower(exclude_sheet)]
data_list <- map(set_names(sheets_to_load), ~ read_excel(file_path, sheet = .x))
combined_df <- bind_rows(data_list, .id = "source_sheet")

# Clean " Tumor_location " if present (leading/trailing spaces)
col_trimmed <- tibble(original = colnames(combined_df), trimmed = trimws(colnames(combined_df)))
if ("Tumor_location" %in% col_trimmed$trimmed) {
  old <- col_trimmed$original[col_trimmed$trimmed == "Tumor_location"][1]
  combined_df <- combined_df %>% rename(Tumor_location = !!sym(old))
}

# Detect Sample ID and add Sample_ID column if not present
possible_id_cols <- colnames(combined_df)[str_detect(tolower(colnames(combined_df)), "sample|^id$|^sid$")]
chosen_id <- NULL
for (p in c("Sample_ID","sample_id","SampleID","sampleid","sample","Sample","ID","id")) {
  if (p %in% colnames(combined_df)) { chosen_id <- p; break }
}
if (!is.null(chosen_id)) {
  combined_df <- combined_df %>% rename(Sample_ID = !!sym(chosen_id))
} else { combined_df <- combined_df %>% mutate(Sample_ID = paste0("s", row_number())) }


# Cohort detection: HMF -> HMF, SA -> PCAWG
combined_df <- combined_df %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    cohort = case_when(
      str_detect(Sample_ID, regex("^\\s*HMF", ignore_case = TRUE)) ~ "HMF",
      str_detect(Sample_ID, regex("^\\s*SA",  ignore_case = TRUE)) ~ "PCAWG",
      TRUE ~ "Other"
    ),
    Tumor_location = if_else(is.na(Tumor_location) | trimws(as.character(Tumor_location)) == "", "Unknown", as.character(Tumor_location))
  )

# ---- 2) Detect signature columns and SV / indel columns ----
all_cols <- colnames(combined_df)

# Signature patterns: typical names include SBS, DBS, ID, "signature", "Sig"
sig_patterns <- c("^SBS", "^DBS", "^ID", "signature", "Signature", "Sig")
sig_cols <- unique(unlist(map(sig_patterns, ~ grep(.x, all_cols, value = TRUE, ignore.case = TRUE))))

# Also include columns named like "SBS1","SBS2"... if not matched earlier
if (length(sig_cols) == 0) {
  sig_cols <- all_cols[str_detect(all_cols, regex("SBS|DBS|ID|signature|Sig", ignore_case = TRUE))]
}

# SV / indel patterns
sv_patterns <- c("sv_", "^sv", "structural", "indel", "indels", "DEL_", "INS_", "SV", "svcount", "sv_count", "indel_count")
sv_cols <- unique(unlist(map(sv_patterns, ~ grep(.x, all_cols, value = TRUE, ignore.case = TRUE))))

# If no sv_cols found, try common names
if (length(sv_cols) == 0) {
  sv_cols <- all_cols[str_detect(all_cols, regex("indel|del|ins|structural|sv|SV", ignore_case = TRUE))]
}


# Numeric filter function: keep only columns that can be coerced to numeric and have > 5 numeric values
keep_numeric <- function(cols) {
  cols[sapply(cols, function(cn) {
    x <- combined_df[[cn]]
    num_ok <- suppressWarnings(sum(!is.na(as.numeric(as.character(x)))))
    !is.null(x) && num_ok >= max(5, 0.01 * nrow(combined_df))
  })]
}

sig_cols <- keep_numeric(sig_cols)
sv_cols <- keep_numeric(sv_cols)

message("Detected ", length(sig_cols), " signature columns, ", length(sv_cols), " SV/indel columns.")

# Create numeric matrices
sig_mat <- combined_df %>% select(all_of(sig_cols)) %>% mutate(across(everything(), ~ as.numeric(as.character(.x))))
sv_mat <- combined_df %>% select(all_of(sv_cols)) %>% mutate(across(everything(), ~ as.numeric(as.character(.x))))

# ---- Q4: Which cancer types exhibit the highest SVs/indels? ----
# If multiple sv cols, compute SV_burden as rowSum; else use single col
if (ncol(sv_mat) == 0) {
  message("No SV/indel columns detected. Q4 and Q5 require SV columns. Skipping Q4/Q5.")
} else {
  combined_df <- combined_df %>% mutate(SV_burden = rowSums(sv_mat, na.rm = TRUE))
  
  # Summary by Tumor_location
  sv_by_tumor <- combined_df %>% group_by(Tumor_location) %>% summarise(n = n(), median_SV = median(SV_burden, na.rm = TRUE), mean_SV = mean(SV_burden, na.rm = TRUE)) %>% arrange(desc(mean_SV))
  write_csv(sv_by_tumor, "sv_burden_by_tumor.csv")
  
  # Plot: top tumor types by mean SV burden
  top_tumors <- sv_by_tumor %>% slice_max(mean_SV, n = top_n) %>% pull(Tumor_location)
  plot_df <- combined_df %>% filter(Tumor_location %in% top_tumors)
  
  p_sv_box <- ggplot(plot_df, aes(x = fct_reorder(Tumor_location, SV_burden, .fun = median), y = SV_burden)) +
    geom_violin(fill = "lightblue", alpha = 0.6) +
    geom_boxplot(width = 0.12, outlier.size = 0.7) +
    coord_flip() +
    labs(title = paste("SV / indel burden by Tumor_location (top", length(top_tumors), ")"), x = "", y = "SV burden (sum)") +
    theme_minimal()
  ggsave("SV_burden_by_tumor_violin.png", p_sv_box, width = 10, height = 8, dpi = 300)
  print(p_sv_box)
  
  # Statistical test across tumor types: Kruskal-Wallis (non-parametric)
  kw <- kruskal.test(SV_burden ~ Tumor_location, data = combined_df %>% filter(Tumor_location %in% top_tumors))
  sink("SV_burden_kruskal.txt"); print(kw); sink()
  
  # pairwise Wilcoxon for the top few tumor types (optional)
  # We'll compute pairwise Wilcoxon and output matrix of p-values
  pairwise_p <- combined_df %>% filter(Tumor_location %in% top_tumors) %>%
    group_by(Tumor_location) %>% summarise(vals = list(SV_burden)) %>% ungroup()
  # simple pairwise matrix
  tumor_levels <- pairwise_p$Tumor_location
  pair_mat <- matrix(NA, nrow = length(tumor_levels), ncol = length(tumor_levels), dimnames = list(tumor_levels, tumor_levels))
  for (i in seq_along(tumor_levels)) {
    for (j in seq_along(tumor_levels)) {
      if (i < j) {
        xi <- unlist(pairwise_p$vals[i])
        xj <- unlist(pairwise_p$vals[j])
        if (length(xi) >= 5 & length(xj) >=5) {
          pair_mat[i,j] <- wilcox.test(xi, xj)$p.value
        }
      }
    }
  }
  write_csv(as_tibble(pair_mat, rownames = "tumor"), "sv_pairwise_pvalues_top_tumors.csv")
}

# ---- Q3: Which mutational signatures co-occur most strongly? Do they form modules? ----
if (ncol(sig_mat) < 2) {
  message("Not enough signature columns detected for co-occurrence analysis. Skipping Q3.")
} else {
  # compute correlation matrix (Spearman preferred for counts/distributions)
  sig_mat_num <- sig_mat %>% mutate(across(everything(), ~ as.numeric(as.character(.x))))
  cor_sig <- cor(sig_mat_num, use = "pairwise.complete.obs", method = "spearman")
  cor_sig[is.na(cor_sig)] <- 0
  
  # find strongly correlated pairs
  th <- sig_corr_threshold
  upper_idx <- which(abs(cor_sig) >= th & upper.tri(cor_sig), arr.ind = TRUE)
  edges_sig <- tibble(
    feature1 = rownames(cor_sig)[upper_idx[,1]],
    feature2 = colnames(cor_sig)[upper_idx[,2]],
    cor = cor_sig[upper_idx]
  ) %>% arrange(desc(abs(cor)))
  write_csv(edges_sig, "signature_strong_pairs.csv")
  message(nrow(edges_sig), " signature pairs with |rho| >= ", th)
  
  # Build network graph and plot with ggraph
  if (nrow(edges_sig) > 0) {
    # nodes
    nodes <- tibble(name = unique(c(edges_sig$feature1, edges_sig$feature2)))
    edges_igraph <- edges_sig %>% rename(from = feature1, to = feature2, weight = cor)
    g_sig <- graph_from_data_frame(edges_igraph, vertices = nodes, directed = FALSE)
    
    # community detection
    comm <- cluster_louvain(g_sig, weights = abs(E(g_sig)$weight))
    V(g_sig)$module <- as.factor(membership(comm))
    V(g_sig)$degree <- degree(g_sig)
    
    # plot
    p_sig_net <- ggraph(g_sig, layout = "fr") +
      geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.8) +
      scale_edge_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0) +
      geom_node_point(aes(size = degree, color = module)) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      labs(title = paste0("Signature co-occurrence network (|rho| >= ", th, ")")) +
      theme_void()
    ggsave("signature_cooccurrence_network.png", p_sig_net, width = 12, height = 10, dpi = 300)
    print(p_sig_net)
  } else {
    message("No signature pairs pass threshold ", th)
  }
  
  # Also produce clustered heatmap of signature correlations (use all signatures)
  # reorder by hierarchical clustering
  d <- as.dist(1 - abs(cor_sig))
  hc <- hclust(d, method = "average")
  cor_sig_ordered <- cor_sig[hc$order, hc$order]
  pheatmap(cor_sig_ordered, main = "Signature correlation (Spearman)", 
           color = colorRampPalette(c("blue","white","red"))(50), show_rownames = TRUE, show_colnames = TRUE,
           filename = "signature_correlation_heatmap.png", width = 10, height = 10)
}


# ---- Q1: Are certain signatures enriched in HMF vs PCAWG? ----
if (ncol(sig_mat) < 1) {
  message("No signature columns for enrichment tests. Skipping Q1.")
} else {
  sig_df <- sig_mat_num %>% as_tibble() %>% bind_cols(combined_df %>% select(Sample_ID, cohort, Tumor_location))
  # filter to PCAWG and HMF
  sig_df_hp <- sig_df %>% filter(cohort %in% c("HMF","PCAWG"))
  # for each signature compute log2 fold-change (HMF / PCAWG) and Wilcoxon p-value
  res_list <- map_dfr(colnames(sig_mat_num), function(feat) {
    xh <- sig_df_hp %>% filter(cohort == "HMF") %>% pull(feat) %>% as.numeric()
    xp <- sig_df_hp %>% filter(cohort == "PCAWG") %>% pull(feat) %>% as.numeric()
    # compute medians; avoid zero divide
    med_h <- median(xh, na.rm = TRUE); med_p <- median(xp, na.rm = TRUE)
    lf <- log2((med_h + 1e-6) / (med_p + 1e-6))
    # wilcox
    pval <- tryCatch(wilcox.test(xh, xp)$p.value, error = function(e) NA_real_)
    tibble(feature = feat, log2FC = lf, pval = pval)
  })
  res_list <- res_list %>% mutate(padj = p.adjust(pval, method = "BH")) %>% arrange(padj)
  write_csv(res_list, "signature_enrichment_HMF_vs_PCAWG.csv")
  
  # Volcano plot: log2FC vs -log10(padj)
  res_plot <- res_list %>% mutate(signif = ifelse(padj < 0.05, "yes", "no"))
  p_volcano <- ggplot(res_plot, aes(x = log2FC, y = -log10(padj+1e-12), label = feature)) +
    geom_point(aes(color = signif), alpha = 0.7) +
    scale_color_manual(values = c("yes" = "red", "no" = "grey60")) +
    geom_text_repel(data = res_plot %>% filter(padj < 0.01 & abs(log2FC) > 0.5), size = 3) +
    labs(title = "Signature enrichment HMF vs PCAWG (log2FC vs -log10 padj)", x = "log2 median(HMF / PCAWG)", y = "-log10(adj p-value)") +
    theme_minimal()
  ggsave("signature_enrichment_volcano.png", p_volcano, width = 9, height = 6, dpi = 300)
  print(p_volcano)
  
  # Boxplots for top enriched features (by padj)
  top_enriched <- res_list %>% filter(!is.na(padj)) %>% slice_min(padj, n = 8) %>% pull(feature)
  if (length(top_enriched) > 0) {
    long_enriched <- sig_df_hp %>% pivot_longer(cols = all_of(top_enriched), names_to = "feature", values_to = "value")
    p_enriched_box <- ggplot(long_enriched, aes(x = cohort, y = value, fill = cohort)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.12) +
      facet_wrap(~ feature, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = c("PCAWG" = "#0072B2","HMF" = "#D55E00")) +
      theme_minimal() +
      labs(title = "Top signature differences HMF vs PCAWG")
    ggsave("top_signature_differences_boxplots.png", p_enriched_box, width = 12, height = 8, dpi = 300)
    print(p_enriched_box)
  }
}

# ---- Q2: Do cancer types cluster into mechanistic groups based on signatures? ----
if (ncol(sig_mat) < 2) {
  message("Not enough signature columns to cluster cancer types. Skipping Q2.")
} else {
  # Aggregate mean signature per tumor type (Tumor_location)
  sig_long <- combined_df %>% select(Sample_ID, Tumor_location, cohort) %>% bind_cols(sig_mat_num)
  mean_by_tumor <- sig_long %>% group_by(Tumor_location) %>% summarise(across(all_of(colnames(sig_mat_num)), ~ mean(.x, na.rm = TRUE)))
  # scale signatures per tumor for visualization (rows=Tumor_location)
  mat_tumor <- mean_by_tumor %>% column_to_rownames("Tumor_location") %>% as.matrix()
  mat_tumor_scaled <- t(scale(t(mat_tumor)))  # center/scale across tumors per signature
  mat_tumor_scaled[is.na(mat_tumor_scaled)] <- 0
  
  # Heatmap with clustering
  pheatmap(mat_tumor_scaled, show_rownames = TRUE, show_colnames = TRUE, 
           clustering_method = "ward.D2", main = "Mean signature profile by Tumor_location (scaled)",
           filename = "mean_signature_by_tumor_heatmap.png", width = 12, height = 10)
  
  # UMAP on sample-level signatures colored by Tumor_location (for sample-level clustering visual)
  umap_res <- uwot::umap(sig_mat_num, n_neighbors = 15, min_dist = 0.3, metric = "cosine")
  umap_df <- as_tibble(umap_res) %>% setNames(c("UMAP1","UMAP2")) %>% bind_cols(combined_df %>% select(Sample_ID, Tumor_location, cohort))
  
  # Plot UMAP but show only top tumor types for readability
  top_tumors_by_n <- combined_df %>% count(Tumor_location, sort = TRUE) %>% slice_head(n = top_n) %>% pull(Tumor_location)
  umap_plot_df <- umap_df %>% mutate(Tumor_location_filtered = if_else(Tumor_location %in% top_tumors_by_n, Tumor_location, "Other"))
  
  p_umap <- ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2, color = Tumor_location_filtered)) +
    geom_point(size = 1.2, alpha = 0.7) +
    labs(title = "UMAP of samples by signature (top tumor labels shown)") +
    theme_minimal() +
    theme(legend.position = "right")
  ggsave("UMAP_signatures_by_tumor.png", p_umap, width = 10, height = 8, dpi = 300)
  print(p_umap)
  
  # Optionally run hierarchical clustering of mean_by_tumor and cut into groups
  hc <- hclust(dist(mat_tumor_scaled), method = "ward.D2")
  tumor_clusters <- cutree(hc, k = cluster_k)
  mean_by_tumor <- mean_by_tumor %>% mutate(cluster = tumor_clusters[match(Tumor_location, names(tumor_clusters))])
  write_csv(mean_by_tumor %>% rowwise() %>% mutate(Tumor_location = rownames(mat_tumor)), "mean_signature_by_tumor_with_cluster.csv")
}

# ---- Q5: Does SV complexity differ between PCAWG and HMF? ----
# Define SV complexity metrics:
# - SV_burden (already computed if sv columns exist)
# - number of unique SV types if columns for specific SV types (DEL, INS, DUP, INV) exist
if (exists("sv_mat") && ncol(sv_mat) > 0) {
  # attempt to detect type-specific columns
  type_patterns <- c("del","ins","dup","inv","tra","trans")
  type_cols <- tibble(col = colnames(sv_mat)) %>% mutate(type = map_chr(col, ~ {
    t <- NA_character_
    for (p in type_patterns) if (str_detect(tolower(.x), p)) { t <- p; break }
    t
  }))
  # compute per-sample type counts if any identified
  if (any(!is.na(type_cols$type))) {
    type_map <- type_cols %>% filter(!is.na(type)) %>% group_by(type) %>% summarise(cols = list(col))
    # add counts
    for (r in seq_len(nrow(type_map))) {
      t <- type_map$type[r]
      cs <- unlist(type_map$cols[r])
      newname <- paste0("SV_type_count_", t)
      combined_df[[newname]] <- rowSums(as.matrix(combined_df %>% select(all_of(cs))), na.rm = TRUE)
    }
  }
  
  # Compute complexity metric: SV_burden * number_of_SV_types_present (simple heuristic)
  # number_of_SV_types_present: count type-specific columns >0 across DEL/INS/DUP/INV
  type_count_cols <- colnames(combined_df)[str_detect(colnames(combined_df), "^SV_type_count_")]
  if (length(type_count_cols) > 0) {
    combined_df <- combined_df %>% mutate(SV_type_present = rowSums(across(all_of(type_count_cols), ~ .x > 0), na.rm = TRUE))
    combined_df <- combined_df %>% mutate(SV_complexity = SV_burden * SV_type_present)
  } else {
    combined_df <- combined_df %>% mutate(SV_complexity = SV_burden)  # fallback
  }
  
  # Compare SV_complexity between cohorts
  df_svcomp <- combined_df %>% filter(cohort %in% c("HMF","PCAWG")) %>% mutate(cohort = factor(cohort, levels = c("PCAWG","HMF")))
  # Wilcoxon test
  test_svcomp <- wilcox.test(SV_complexity ~ cohort, data = df_svcomp)
  sink("SV_complexity_test.txt"); print(test_svcomp); sink()
  
  # Plot violin/box
  p_svcomp <- ggplot(df_svcomp, aes(x = cohort, y = SV_complexity, fill = cohort)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    scale_y_continuous(trans = "log10", labels = comma) +
    scale_fill_manual(values = c("PCAWG" = "#0072B2","HMF" = "#D55E00")) +
    labs(title = "SV complexity by cohort (PCAWG vs HMF)", y = "SV complexity (log10 scale)") +
    theme_minimal()
  ggsave("SV_complexity_cohort_violin.png", p_svcomp, width = 8, height = 6, dpi = 300)
  print(p_svcomp)
} else {
  message("SV columns not detected; Q5 skipped.")
}


