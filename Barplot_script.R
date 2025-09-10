# ============================
# GO Enrichment Analysis - Top 20 Terms
# ============================

# --- Load GO results ---
# Expecting a table with columns: GO, Log.q.value., LogP, FDR
go_df <- read.csv("GO.csv", header = TRUE, stringsAsFactors = FALSE) 

# --- Filter terms ---
# Keep only significant terms (Log.q.value. < -2 means q < 0.01 on log scale)
filtered_genes <- go_df[go_df$Log.q.value. < -2, ]

# --- Select top 20 terms by -log10(p-value) ---
top_20 <- filtered_genes %>%
  mutate(LogP = abs(LogP)) %>%              # ensure positivity for plotting
  arrange(desc(LogP)) %>%                   # sort descending by significance
  slice_head(n = 20) %>%                    # keep top 20
  mutate(
    Description = factor(GO, levels = rev(GO)),  # reverse for coord_flip()
    FDR_Category = case_when(
      FDR < 0.01 ~ "FDR < 0.01",
      FDR < 0.05 ~ "FDR < 0.05",            
      TRUE ~ "Not Significant"
    )
  )

# --- Plot ---
ggplot(top_20, aes(x = Description, y = LogP, fill = FDR_Category)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) -x,  # if LogP was originally negative, flip sign
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = c(
    "FDR < 0.01" = "#F4A261",   # orange
    "FDR < 0.05" = "#e75480",   # pink
    "Not Significant" = "#8DA399" # muted green-grey
  )) +
  labs(
    #title = "Top Enriched GO Terms",    # optional title
    x = "",
    y = "-log10(p-value)",
    fill = "FDR Significance"
  ) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.y  = element_text(size = 15, hjust = 0, lineheight = 1.2), # more spacing
    legend.position = c(0.65, 0.2),  # inside plot
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "grey", size = 0.3),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )