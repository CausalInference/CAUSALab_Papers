source("scripts/load_globals.R")


# load and process models -----------------------------------------------------------------------------------------------------

models <- list()

# load models
for (analysis in c("singleton", "plevel", "restricted")) {
  if (analysis %in% c("singleton")) {
    multiple_stub <- "_probable_multiple"
  } else {
    multiple_stub <- ""
  }
  for (comparison in c("_artiui_wcalyear", "_nonuse_wcalyear")) {
    for (outcome in c("malformation", "nicu_infant")) {
      if (comparison %in% c("_artiui_wcalyear")) {
        exposure <- "art_only"
      } else {
        exposure <- "art"
      }
      if (analysis %in% c("restricted")) {
        obtain_estimate <- get_model_estimate
      } else {
        obtain_estimate <- get_robust_estimate
      }
      models[[glue("{analysis}_crude_{outcome}{comparison}")]] <- readRDS(file=glue("models/{analysis}_crude_pois_model{multiple_stub}_{outcome}{comparison}.rds")) %>% 
        obtain_estimate(exposure) %>% mutate(type = "Crude", comparison=comparison, outcome=outcome, analysis=analysis)
      models[[glue("{analysis}_adjusted_{outcome}{comparison}")]] <- readRDS(file=glue("models/{analysis}_adj_pois_model{multiple_stub}_{outcome}{comparison}.rds")) %>% 
        obtain_estimate(exposure) %>% mutate(type = "Adjusted", comparison=comparison, outcome=outcome, analysis=analysis)
    }
  }
}

# combine model outputs, exponentiate RRs, tidy variable values
model_outputs <- models %>% bind_rows() %>% 
  mutate(across(c(mean, lower, upper), ~ exp(.))) %>%
  mutate(rr = as.character(glue("{scales::number(mean, accuracy=0.01)} ({scales::number(lower, accuracy=0.01)}-{scales::number(upper, accuracy=0.01)})")))  %>%
  mutate(comparison = case_when(comparison == "_artiui_wcalyear" ~ "IUI", comparison == "_nonuse_wcalyear" ~ "Non-use", .default="Error"),
         outcome = case_when(outcome == "malformation" ~ "Malformation", outcome == "nicu_infant" ~ "NICU admission"),
         analysis = case_when(analysis=="singleton" ~ "Singletons", analysis=="plevel" ~ "Pregnancy-level", analysis=="restricted" ~ "Infant-level IPSW")) 

# arrange outputs
model_outputs <- model_outputs %>% arrange(desc(comparison), outcome, desc(analysis))

# prepare forest
extra_prep_non_use_adj <- prepare_forest(model_outputs, "Adjusted", "Non-use")

# specify forest header
header <- tibble(
  outcome = c("Outcome"),
  analysis = c("Analysis"),
  rr = c("Risk ratio (95% CI)"),
  summary=TRUE,
)

# plot forest plot non-use adjusted 
fp_non_use_adjusted <- bind_rows(header, extra_prep_non_use_adj) %>% forestplot(labeltext=c(outcome, analysis, rr),
                                                                                is.summary=c(TRUE, rep(FALSE, nrow(model_outputs)/2)),
                                                                                title="ART vs. Non-use",
                                                                                graph.pos=4,
                                                                                xlog=TRUE,
                                                                                xticks = c(log(0.71), log(1), log(2.5)),
                                                                                boxsize = 0.1,
                                                                                txt_gp = fpTxtGp(ticks=gpar(cex=1), label=gpar(cex=1), summary=gpar(cex=1)),
                                                                                col = fpColors(box = "royalblue",
                                                                                               line = "darkblue", 
                                                                                               summary = "royalblue",
                                                                                               hrz_lines = "#444444"),
                                                                                colgap = unit(0.01,"npc"),
                                                                                align=c("l", "l", "l")) 

# prepare forest plot IUI adjusted 
extra_prep_iui_adj <- prepare_forest(model_outputs, "Adjusted", "IUI")

# specify simplified header for forest plot
header_simple <- tibble(
  rr = c("Risk ratio (95% CI)"),
  summary=TRUE,
)

# plot forest plot IUI adjusted
fp_iui_adjusted <- bind_rows(header_simple, extra_prep_iui_adj) %>% forestplot(labeltext=c(rr),
                                                                                       is.summary=c(TRUE, rep(FALSE, nrow(model_outputs)/2)),
                                                                                       title="ART vs. IUI",
                                                                                       graph.pos=2,
                                                                                       xlog=TRUE,
                                                                                       xticks = c(log(0.71), log(1), log(2.5)),
                                                                                       boxsize = 0.1,
                                                                                       txt_gp = fpTxtGp(ticks=gpar(cex=1), label=gpar(cex=1), summary=gpar(cex=1)),
                                                                                       col = fpColors(box = "royalblue",
                                                                                                      line = "darkblue", 
                                                                                                      summary = "royalblue",
                                                                                                      hrz_lines = "#444444"),
                                                                                       colgap = unit(0.01,"npc"),
                                                                                       align=c("l")) 

# output as JPEG file (FIGURE 3)
jpeg(glue("output/rr_adjusted_combined_forest_plot.jpeg"), units="mm", width=200, height=80, res=300)
plot_forests_side_by_side(fp_non_use_adjusted, fp_iui_adjusted, 0.66) 
dev.off()

# prepare forest non-use crude
extra_prep_non_use_crude <- prepare_forest(model_outputs, "Crude", "Non-use")

# plot forest non-use crude
fp_non_use_crude <- bind_rows(header, extra_prep_non_use_crude) %>% forestplot(labeltext=c(outcome, analysis, rr),
                                                                                                                  is.summary=c(TRUE, rep(FALSE, nrow(model_outputs)/2)),
                                                                                                                  title="ART vs. Non-use",
                                                                                                                  graph.pos=4,
                                                                                                                  xlog=TRUE,
                                                                                                                  xticks = c(log(0.71), log(1), log(2.5)),
                                                                                                                  boxsize = 0.1,
                                                                                                                  txt_gp = fpTxtGp(ticks=gpar(cex=1), label=gpar(cex=1), summary=gpar(cex=1)),
                                                                                                                  col = fpColors(box = "royalblue",
                                                                                                                                 line = "darkblue", 
                                                                                                                                 summary = "royalblue",
                                                                                                                                 hrz_lines = "#444444"),
                                                                                                            colgap = unit(0.01,"npc"),
                                                                                                            align=c("l", "l", "l")) 

# prepare forest IUI crude
extra_prep_iui_crude <- prepare_forest(model_outputs, "Crude", "IUI")

# plot forest IUI crude
fp_iui_crude <- bind_rows(header_simple, extra_prep_iui_crude) %>% forestplot(labeltext=c(rr),
                                                                                               is.summary=c(TRUE, rep(FALSE, nrow(model_outputs)/2)),
                                                                                               title="ART vs. IUI",
                                                                                               graph.pos=2,
                                                                                               xlog=TRUE,
                                                                                               xticks = c(log(0.71), log(1), log(2.5)),
                                                                                               boxsize = 0.1,
                                                                                               txt_gp = fpTxtGp(ticks=gpar(cex=1), label=gpar(cex=1), summary=gpar(cex=1)),
                                                                                               col = fpColors(box = "royalblue",
                                                                                                              line = "darkblue", 
                                                                                                              summary = "royalblue",
                                                                                                              hrz_lines = "#444444"),
                                                                                         colgap = unit(0.01,"npc"),
                                                                                         align=c("l"))


# output as JPEG (APPENDIX FIGURE 2)
jpeg(glue("output/rr_crude_combined_forest_plot.jpeg"), units="mm", width=200, height=80, res=300)
plot_forests_side_by_side(fp_non_use_crude, fp_iui_crude, 0.66) 
dev.off()

# plot forests side by side
plot_forests_side_by_side(fp_non_use_crude, fp_iui_crude, 0.665) 
crude_plot <- grid.grab(wrap.grobs = TRUE)
plot_forests_side_by_side(fp_non_use_adjusted, fp_iui_adjusted, 0.665) 
adj_plot <- grid.grab(wrap.grobs = TRUE)

# output combined plot as JPEG 
jpeg(glue("output/rr_both_combined_forest_plot.jpeg"), units="mm", width=200, height=160, res=300)
ggarrange(crude_plot, adj_plot, ncol=1, nrow =2, labels=c("A","B"), hjust=-0.2)
dev.off()
