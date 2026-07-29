# ============================================================================
#  Quick GLMM Review  —  R Shiny app (glmmTMB / DHARMa / emmeans)
# ----------------------------------------------------------------------------
#  Extends the "Quick Mixed-Model Review" (lmerTest) app to generalized linear
#  mixed models via glmmTMB. Same overall structure and UF/IFAS look, but:
#
#    * family        : dropdown of glmmTMB distribution families
#    * ziformula      : optional zero-inflation model (checkbox + fixed-effect
#                       picker); built as ~0 (off) / ~1 (intercept only) /
#                       ~var1 + var2 ...
#    * dispformula    : optional dispersion model, 0, 1, or 2 predictors
#                       (can reuse any created interaction/combo variable)
#    * Residuals      : DHARMa simulated residuals (QQ + resid-vs-fitted),
#                       DHARMa overdispersion test, zero-inflation test,
#                       outlier test  (glmmTMB residuals are NOT plain
#                       Pearson/deviance residuals in the lmer sense, so
#                       DHARMa's simulation-based approach is used instead)
#    * Binary outcomes : handled on their OWN tab, because 0/1 (Bernoulli)
#                       data need different guardrails than general GLMMs
#                       (no dispersion model, restricted family list, link
#                       choice, and different DHARMa expectations). Keeping
#                       this separate avoids the "it silently did something
#                       odd with my 0/1 column" failure mode. Supports both
#                       individual 0/1 rows AND grouped successes/total-trials
#                       counts (failures = total - successes, computed for you).
#    * Ordered beta    : "Ordered beta (logit link)" on the General GLMM tab
#                       covers proportions that include exact 0s and/or 1s,
#                       which the standard Beta family cannot (it requires
#                       values strictly inside (0,1)).
#
#  Install once:
#    install.packages(c("shiny","shinyjs","bslib","glmmTMB","DHARMa","emmeans",
#                       "multcomp","multcompView","ggplot2","DT","car"))
#  Branding:
#    Uses the shared UF/IFAS components.R if present (same as the lmer app),
#    otherwise falls back to an inline copy so the app still runs standalone.
#  Run:
#    shiny::runApp("glmmTMB_app.R")
# ============================================================================

library(shiny)
library(shinyjs)
library(bslib)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)
library(DT)

# --- Shared UF/IFAS look-and-feel (same convention as the lmer app) --------
if (file.exists("components.R")) {
  source("components.R", local = FALSE)
} else {
  UF_BLUE   <- "#003087"
  UF_ORANGE <- "#FA4616"
  UF_COLORS <- c(UF_BLUE, UF_ORANGE, "#2ca25f", "#8856a7")
  uf_theme  <- function() {
    bslib::bs_theme(bootswatch = "flatly", primary = UF_BLUE,
                    secondary = UF_ORANGE, font_scale = 0.95,
                    "navbar-bg" = UF_BLUE) |>
      bslib::bs_add_rules(".shiny-output-error-validation { white-space: nowrap; }")
  }
  uf_logo_uri <- function() NULL
  uf_title <- function(text, logo = uf_logo_uri())
    tags$span(if (!is.null(logo)) tags$img(src = logo, height = "51px", alt = "UF/IFAS",
                                           style = "margin-right:12px; vertical-align:middle;"),
              tags$span(text, style = "vertical-align:middle;"))
  info_tip <- function(..., placement = "right")
    bslib::tooltip(tags$span(icon("circle-question"),
                             style = "color:#aaa;cursor:help;margin-left:5px;font-size:0.82em;"),
                   ..., placement = placement)
  copy_js <- "
function DEcopy(id, btn){
  var el = document.getElementById(id);
  if(!el) return;
  var txt = el.innerText || el.textContent || '';
  navigator.clipboard.writeText(txt).then(function(){
    if(btn){ var o = btn.innerHTML; btn.innerHTML = 'Copied!'; setTimeout(function(){ btn.innerHTML = o; }, 1200); }
  });
}
"
}
HAS_CAR <- requireNamespace("car", quietly = TRUE)
UF_CHARCOAL <- "#222222"
UF_BLUE_LT  <- "#cdd8ea"
MAX_FIXED   <- 3L

# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && !nzchar(a))) b else a
bq <- function(x) sprintf("`%s`", x)

is_categorical <- function(df, vars)
  vapply(vars, function(v) is.factor(df[[v]]) || is.character(df[[v]]), logical(1))

round_df <- function(d, digits = 4) {
  num  <- vapply(d, is.numeric, logical(1))
  is_p <- grepl("p\\.value|p_value|Pr\\(|^p$", names(d), ignore.case = TRUE)
  for (j in which(num)) d[[j]] <- if (is_p[j]) signif(d[[j]], 3) else round(d[[j]], digits)
  d
}

build_random_part <- function(random, slope = NULL, slope_group = NULL) {
  if (length(random) == 0) return("1")   # glmmTMB needs some formula; caller guards
  if (!is.null(slope) && nzchar(slope)) {
    if (is.null(slope_group) || !length(slope_group)) slope_group <- random
    terms <- vapply(random, function(g)
      if (g %in% slope_group) sprintf("(1 + %s | %s)", bq(slope), bq(g))
      else sprintf("(1 | %s)", bq(g)), character(1))
    paste(terms, collapse = " + ")
  } else paste(sprintf("(1 | %s)", bq(random)), collapse = " + ")
}

build_fixed_part <- function(fixed, interactions) {
  if (length(fixed) == 0) return("1")
  if (isTRUE(interactions) && length(fixed) > 1) paste(bq(fixed), collapse = " * ")
  else paste(bq(fixed), collapse = " + ")
}

# Main conditional (fixed + random) formula, e.g.  y ~ A * B + (1 | Block)
# `response` is normally a bare column name (quoted internally via bq()), but
# for grouped-binomial fits the caller passes an already-built expression like
# "cbind(`succ`, `tot` - `succ`)" and sets quote_response = FALSE so it is
# used as-is instead of being wrapped in backticks.
build_formula_string <- function(response, fixed, random, interactions,
                                 slope = NULL, slope_group = NULL,
                                 quote_response = TRUE) {
  resp_part   <- if (quote_response) bq(response) else response
  fixed_part  <- build_fixed_part(fixed, interactions)
  random_part <- build_random_part(random, slope, slope_group)
  paste(resp_part, "~", fixed_part, "+", random_part)
}

# Zero-inflation / dispersion side formulas: character vector of predictor
# names -> "~0" (off), "~1" (intercept only), or "~v1 + v2"
build_side_formula <- function(vars, enabled = TRUE) {
  if (!enabled) return("~0")
  if (length(vars) == 0) return("~1")
  paste("~", paste(bq(vars), collapse = " + "))
}

# Family constructor lookup: label shown in UI -> glmmTMB family object + a
# note on when it's appropriate
FAMILY_CHOICES <- c(
  "Gaussian (identity link)"          = "gaussian_identity",
  "Gaussian (log link)"               = "gaussian_log",
  "Gamma (log link)  \u2013 right-skewed positive" = "Gamma_log",
  "Poisson (log link)  \u2013 counts" = "poisson_log",
  "Negative binomial, nbinom2 (log link)  \u2013 overdispersed counts" = "nbinom2_log",
  "Negative binomial, nbinom1 (log link)"           = "nbinom1_log",
  "Tweedie (log link)  \u2013 zero-heavy continuous" = "tweedie_log",
  "Beta (logit link)  \u2013 proportions in (0,1)"  = "beta_logit",
  "Ordered beta (logit link)  \u2013 proportions in [0,1], INCLUDES 0 and 1" = "ordbeta_logit"
)
family_from_key <- function(key) {
  switch(key,
         "gaussian_identity" = gaussian(link = "identity"),
         "gaussian_log"       = gaussian(link = "log"),
         "Gamma_log"          = Gamma(link = "log"),
         "poisson_log"        = poisson(link = "log"),
         "nbinom2_log"        = glmmTMB::nbinom2(link = "log"),
         "nbinom1_log"        = glmmTMB::nbinom1(link = "log"),
         "tweedie_log"        = glmmTMB::tweedie(link = "log"),
         "beta_logit"         = glmmTMB::beta_family(link = "logit"),
         "ordbeta_logit"      = glmmTMB::ordbeta(link = "logit"),
         gaussian())
}
# Matching valid-R constructor text for the "R code" / formula display panels,
# so the code a user copies actually runs standalone (rather than showing the
# internal dropdown key).
FAMILY_CODE_TEXT <- c(
  "gaussian_identity" = "gaussian(link = 'identity')",
  "gaussian_log"       = "gaussian(link = 'log')",
  "Gamma_log"          = "Gamma(link = 'log')",
  "poisson_log"        = "poisson(link = 'log')",
  "nbinom2_log"        = "glmmTMB::nbinom2(link = 'log')",
  "nbinom1_log"        = "glmmTMB::nbinom1(link = 'log')",
  "tweedie_log"        = "glmmTMB::tweedie(link = 'log')",
  "beta_logit"         = "glmmTMB::beta_family(link = 'logit')",
  "ordbeta_logit"      = "glmmTMB::ordbeta(link = 'logit')"
)
BINARY_LINK_CHOICES <- c("logit" = "logit", "probit" = "probit",
                         "cloglog" = "cloglog", "cauchit" = "cauchit")

# Plain-language valid-domain note per family, shown live under the family
# picker so a mismatched response (e.g. negative numbers with Gamma, or
# proportions with Poisson) is caught before fitting rather than after.
FAMILY_DOMAIN_NOTE <- c(
  "gaussian_identity" = "Any real number \u2014 positive, negative, or zero.",
  "gaussian_log"       = "Response must be strictly greater than 0 (the log link is undefined at or below 0).",
  "Gamma_log"          = "Response must be strictly greater than 0. Gamma is undefined at exactly 0 \u2014 if your data includes true zeros, this is the wrong family.",
  "poisson_log"        = "Non-negative integer counts: 0, 1, 2, 3, ... No negative values, no decimals.",
  "nbinom2_log"        = "Non-negative integer counts (0, 1, 2, ...). Use this over Poisson when counts are overdispersed (variance > mean).",
  "nbinom1_log"        = "Non-negative integer counts (0, 1, 2, ...). An alternative overdispersion structure to nbinom2 (variance scales linearly with the mean).",
  "tweedie_log"        = "Non-negative continuous values, and CAN include exact zeros (e.g. rainfall, biomass with true absences).",
  "beta_logit"         = "Strictly between 0 and 1 \u2014 EXCLUDES exact 0 and exact 1. Rescale data that touches the boundaries, e.g. y' = (y*(n-1)+0.5)/n.",
  "ordbeta_logit"      = "The CLOSED interval [0,1] \u2014 CAN include exact 0s and exact 1s (that's the whole point of this family vs. standard Beta). Use this instead of rescaling boundary values away."
)

# ---------------------------------------------------------------------------
#  Built-in example data: counts + a zero-inflated count + a binary outcome
# ---------------------------------------------------------------------------
make_example_data <- function() {
  set.seed(2024)
  d <- expand.grid(Site = factor(1:8), Treatment = c("Control", "Fertilized", "Grazed"),
                   Season = c("Spring", "Fall"), KEEP.OUT.ATTRS = FALSE,
                   stringsAsFactors = FALSE)
  d <- d[rep(seq_len(nrow(d)), each = 3), ]; rownames(d) <- NULL
  n <- nrow(d)
  site_eff <- setNames(rnorm(8, 0, 0.4), levels(d$Site))
  b <- site_eff[as.character(d$Site)]
  trtC <- c(Control = 0, Fertilized = 0.6, Grazed = -0.4)[d$Treatment]
  seaC <- c(Spring = 0, Fall = -0.3)[d$Season]
  
  # overdispersed counts -> nbinom2 target
  mu_count <- exp(1.2 + trtC + seaC + b)
  d$insect_count <- rnbinom(n, mu = mu_count, size = 2)
  
  # zero-heavy counts -> good candidate for ziformula
  zi_p <- plogis(-1 + 0.8 * (d$Treatment == "Grazed"))
  lam  <- exp(1.0 + trtC + b)
  d$seedling_count <- ifelse(rbinom(n, 1, zi_p) == 1, 0, rpois(n, lam))
  
  # proportion in (0,1) -> beta family
  eta <- 0.5 + trtC * 0.5 + seaC * 0.3 + b
  d$cover_prop <- pmin(pmax(plogis(eta) + rnorm(n, 0, 0.03), 0.001), 0.999)
  
  # proportion in [0,1], CAN include exact 0s/1s -> ordered beta family target
  d$cover_prop_ord <- pmin(pmax(plogis(eta) + rnorm(n, 0, 0.12), 0), 1)
  
  # binary presence/absence -> binomial (own tab, individual 0/1 rows)
  peta <- -0.3 + 0.9 * (d$Treatment == "Fertilized") - 0.6 * (d$Treatment == "Grazed") + b
  d$present <- rbinom(n, 1, plogis(peta))
  
  # grouped binomial counts -> successes out of total trials (own tab,
  # "grouped counts" mode; failures = trials - successes computed for you)
  d$trials    <- sample(8:15, n, replace = TRUE)
  d$successes <- rbinom(n, d$trials, plogis(peta))
  
  d
}

code_panel <- function(out_id, label = "R code") {
  tags$details(class = "uf-codewrap",
               tags$summary(tags$b(sprintf("\u25b6 %s (copy & run in R)", label))),
               tags$button("Copy code", class = "btn btn-default uf-copy",
                           onclick = sprintf("DEcopy('%s', this)", out_id)),
               verbatimTextOutput(out_id))
}

css_block <- tags$head(
  tags$style(HTML(sprintf("
  .uf-header { background: %1$s; color: #fff; padding: 14px 20px;
    margin: -12px -12px 16px -12px; border-bottom: 5px solid %2$s; }
  .uf-header .uf-title { font-size: 24px; font-weight: 700; }
  .uf-header .uf-sub { color: #fff; opacity: .85; font-weight: 600; font-size: 13px; }
  pre { border-left: 4px solid %1$s; }
  .alert-warning { border-left: 5px solid %2$s; }
  .alert-info    { border-left: 5px solid %1$s; }
  .uf-dl { margin: 6px 6px 12px 0; }
  .uf-copy { margin: 4px 0 0 0; font-size: 12px; padding: 3px 10px; }
  .uf-codewrap { margin: 6px 0 14px 0; }
  summary { cursor: pointer; user-select: none; padding: 2px 0; }
  summary:hover { color: %2$s; }
  ", UF_BLUE, UF_ORANGE))),
  tags$script(HTML(copy_js)))

# A block of controls reused by BOTH the general and the binary tab (data
# type overrides + combined/interaction-variable creator). Kept as a function
# so the two tabs stay independent but don't duplicate logic by hand.
setup_controls <- function(prefix) {
  tagList(
    tags$details(open = NA,
                 tags$summary(tags$b("\u25be Variable types \u2014 treat numbers as categories")),
                 helpText("If a treatment level or ID is stored as a number, read it as a ",
                          "factor so its levels are compared as groups, not fit as a slope."),
                 uiOutput(paste0(prefix, "_force_factor_ui"))),
    tags$details(open = NA,
                 tags$summary(tags$b("\u25be Create combined variable (interaction)")),
                 helpText("Paste 2+ columns into one new factor \u2014 useful as a ",
                          "dispersion-model predictor (e.g. Treatment.Season)."),
                 uiOutput(paste0(prefix, "_combo_vars_ui")),
                 textInput(paste0(prefix, "_combo_sep"), "Separator", value = "."),
                 actionButton(paste0(prefix, "_combo_add"), "Create variable", class = "btn-primary"),
                 uiOutput(paste0(prefix, "_combo_list_ui")))
  )
}

# ===========================================================================
#  UI
# ===========================================================================
ui <- fluidPage(
  useShinyjs(),
  theme = uf_theme(),
  title = "Quick GLMM Review",
  css_block,
  
  div(class = "uf-header",
      div(class = "uf-title", uf_title("Quick GLMM Review")),
      div(class = "uf-sub", "glmmTMB \u00b7 DHARMa \u00b7 emmeans  |  University of Florida")),
  
  radioButtons("data_source", "Data source",
               choices = c("Example data (counts + proportion + binary)" = "example",
                           "Upload a CSV file" = "upload"), selected = "example", inline = TRUE),
  conditionalPanel("input.data_source == 'upload'",
                   fileInput("file", "Choose CSV file", accept = c(".csv", ".txt")),
                   checkboxInput("header", "Header row", TRUE)),
  tags$hr(),
  
  tabsetPanel(
    id = "main_tabs",
    
    # =======================================================================
    # TAB 1 : General GLMM (counts, proportions, continuous, right-skewed...)
    # =======================================================================
    tabPanel("General GLMM",
             sidebarLayout(
               sidebarPanel(width = 4,
                            setup_controls("g"),
                            tags$hr(),
                            uiOutput("g_response_ui"),
                            selectInput("g_family",
                                        tagList("Family (distribution + link)",
                                                info_tip("Choose the distribution that matches your response's ",
                                                         "valid range (e.g. counts, positive continuous, or a ",
                                                         "proportion). The note below the box updates for whatever ",
                                                         "family is selected and lists exactly what values are valid. ",
                                                         "If your proportion data includes exact 0s or 1s, use ",
                                                         "Ordered beta instead of standard Beta.")),
                                        choices = FAMILY_CHOICES, selected = "nbinom2_log"),
                            uiOutput("g_family_note"),
                            uiOutput("g_fixed_ui"),
                            checkboxInput("g_interactions", "Include interactions among fixed effects", FALSE),
                            uiOutput("g_random_ui"),
                            tags$details(tags$summary(tags$b("\u25b8 Random slope (optional)")),
                                         uiOutput("g_slope_ui"), uiOutput("g_slope_group_ui")),
                            tags$hr(),
                            tags$details(open = NA, tags$summary(tags$b("\u25be Zero-inflation model (ziformula)")),
                                         helpText("Models the probability of a structural extra zero. ",
                                                  "Leave the picker empty for an intercept-only ",
                                                  "zero-inflation model (~1)."),
                                         checkboxInput("g_zi_on", "Include a zero-inflation term", FALSE),
                                         uiOutput("g_zi_vars_ui")),
                            tags$details(open = NA, tags$summary(tags$b("\u25be Dispersion model (dispformula)")),
                                         helpText("0, 1, or 2 predictors for how residual dispersion ",
                                                  "varies (e.g. by Treatment). Created interaction ",
                                                  "variables can be used here too."),
                                         uiOutput("g_disp_vars_ui")),
                            tags$hr(),
                            uiOutput("g_emm_ui"), uiOutput("g_emm_by_ui"),
                            selectInput("g_adjust", "Post-hoc p-value adjustment",
                                        choices = c("Tukey" = "tukey", "Sidak" = "sidak",
                                                    "Bonferroni" = "bonferroni", "Holm" = "holm", "None" = "none"),
                                        selected = "tukey"),
                            sliderInput("g_conf", "Confidence level", min = 0.80, max = 0.99, value = 0.95, step = 0.01),
                            tags$hr(),
                            actionButton("g_run", "Fit model", class = "btn-primary"),
                            actionButton("g_reset", "Reset", class = "btn-default")
               ),
               mainPanel(width = 8,
                         uiOutput("g_message_box"),
                         tabsetPanel(id = "g_result_tabs",
                                     tabPanel("Data", h4("Preview"), DT::dataTableOutput("g_data_head"),
                                              h4("Structure"), verbatimTextOutput("g_data_str")),
                                     tabPanel("Model & code", h4("Formula used"), verbatimTextOutput("g_formula_txt"),
                                              h4("Model summary"), verbatimTextOutput("g_model_summary"),
                                              h4("AIC / BIC / logLik"), DT::dataTableOutput("g_fit_table"),
                                              h4("R code"), verbatimTextOutput("g_code_block"),
                                              tags$button("Copy code", class = "btn btn-default uf-copy",
                                                          onclick = "DEcopy('g_code_block', this)"),
                                              downloadButton("g_dl_code", "Download .R", class = "btn-default uf-dl")),
                                     tabPanel("Wald ANOVA", helpText("Type III Wald chi-square tests (car::Anova). ",
                                                                     "A likelihood-ratio alternative is noted if car is unavailable."),
                                              verbatimTextOutput("g_anova_table")),
                                     tabPanel("DHARMa residuals",
                                              helpText("Simulation-based residuals \u2014 the right tool for glmmTMB ",
                                                       "(ordinary Pearson/deviance residuals from a GLMM are not ",
                                                       "reliably uniform, so this replaces the lmer residual panel)."),
                                              plotOutput("g_dharma_plot", height = "420px"),
                                              downloadButton("g_dl_dharma", "Download PNG", class = "btn-default uf-dl"),
                                              h4("Overdispersion test"), verbatimTextOutput("g_disp_test"),
                                              h4("Zero-inflation test"), verbatimTextOutput("g_zi_test"),
                                              h4("Outlier test"), verbatimTextOutput("g_outlier_test"),
                                              code_panel("g_code_dharma", "DHARMa code")),
                                     tabPanel("EMMeans & post-hoc",
                                              h4("Estimated marginal means (response scale)"), DT::dataTableOutput("g_emm_table"),
                                              downloadButton("g_dl_emm", "Download CSV", class = "btn-default uf-dl"),
                                              h4("Pairwise comparisons"), DT::dataTableOutput("g_pairs_table"),
                                              downloadButton("g_dl_pairs", "Download CSV", class = "btn-default uf-dl"),
                                              h4("Compact letter display"), DT::dataTableOutput("g_cld_table"),
                                              downloadButton("g_dl_cld", "Download CSV", class = "btn-default uf-dl"),
                                              h4("EMMeans plot"), plotOutput("g_emm_plot", height = "460px"),
                                              downloadButton("g_dl_emmplot", "Download PNG", class = "btn-default uf-dl"),
                                              code_panel("g_code_emm", "EMMeans code"))
                         )
               )
             )
    ),
    
    # =======================================================================
    # TAB 2 : Binary / binomial-count GLMM  \u2014 kept fully separate from the
    #         general tab. Supports individual 0/1 rows OR grouped
    #         successes/total-trials counts (failures computed automatically).
    # =======================================================================
    tabPanel("Binary (0/1) GLMM",
             sidebarLayout(
               sidebarPanel(width = 4,
                            helpText(tags$b("Binary & binomial-count outcomes are handled separately."),
                                     " Use individual 0/1 rows, or grouped successes/total-trials counts ",
                                     "below; family is fixed to a Bernoulli/binomial distribution ",
                                     "(dispersion is not estimated for binomial data, so no dispersion-model ",
                                     "picker is shown here). If you instead have a continuous PROPORTION that ",
                                     "can equal exactly 0 or 1 (not a successes/trials count), use the ",
                                     tags$b("Ordered beta"), " family on the General GLMM tab."),
                            setup_controls("b"),
                            tags$hr(),
                            radioButtons("b_resp_type", "Response format",
                                        choices = c("Individual rows (0/1 outcome per row)" = "binary",
                                                    "Grouped counts (successes out of total trials)" = "grouped"),
                                        selected = "binary"),
                            conditionalPanel(
                              "input.b_resp_type == 'binary'",
                              uiOutput("b_response_ui"),
                              helpText(tags$b("Valid values: "), "exactly 0 or 1 (or a two-level factor, ",
                                       "coded as the 2nd level = \"success\")."),
                              uiOutput("b_level_note")
                            ),
                            conditionalPanel(
                              "input.b_resp_type == 'grouped'",
                              helpText("Pick a column of success counts and a column of total trials. ",
                                       "The number of failures is calculated automatically as ",
                                       tags$b("total \u2212 successes"), " \u2014 you don't need to create a ",
                                       "separate failures column yourself."),
                              uiOutput("b_success_ui"),
                              uiOutput("b_total_ui"),
                              uiOutput("b_grouped_note")
                            ),
                            selectInput("b_link",
                                        tagList("Link function",
                                                info_tip("logit: symmetric, most common, coefficients are log-odds. ",
                                                         "probit: symmetric, coefficients on a normal-CDF scale. ",
                                                         "cloglog: asymmetric, useful when 1s are rare. ",
                                                         "cauchit: symmetric with heavier tails, robust to extreme predictors.")),
                                        choices = BINARY_LINK_CHOICES, selected = "logit"),
                            uiOutput("b_fixed_ui"),
                            checkboxInput("b_interactions", "Include interactions among fixed effects", FALSE),
                            uiOutput("b_random_ui"),
                            tags$details(tags$summary(tags$b("\u25b8 Random slope (optional)")),
                                         uiOutput("b_slope_ui"), uiOutput("b_slope_group_ui")),
                            tags$details(open = NA, tags$summary(tags$b("\u25be Zero-inflation model (rare for true 0/1 data)")),
                                         checkboxInput("b_zi_on", "Include a zero-inflation term", FALSE),
                                         uiOutput("b_zi_vars_ui")),
                            tags$hr(),
                            uiOutput("b_emm_ui"), uiOutput("b_emm_by_ui"),
                            selectInput("b_adjust", "Post-hoc p-value adjustment",
                                        choices = c("Tukey" = "tukey", "Sidak" = "sidak",
                                                    "Bonferroni" = "bonferroni", "Holm" = "holm", "None" = "none"),
                                        selected = "tukey"),
                            sliderInput("b_conf", "Confidence level", min = 0.80, max = 0.99, value = 0.95, step = 0.01),
                            tags$hr(),
                            actionButton("b_run", "Fit model", class = "btn-primary"),
                            actionButton("b_reset", "Reset", class = "btn-default")
               ),
               mainPanel(width = 8,
                         uiOutput("b_message_box"),
                         tabsetPanel(id = "b_result_tabs",
                                     tabPanel("Data", h4("Preview"), DT::dataTableOutput("b_data_head"),
                                              h4("Structure"), verbatimTextOutput("b_data_str")),
                                     tabPanel("Model & code", h4("Formula used"), verbatimTextOutput("b_formula_txt"),
                                              h4("Model summary"), verbatimTextOutput("b_model_summary"),
                                              h4("AIC / BIC / logLik"), DT::dataTableOutput("b_fit_table"),
                                              h4("R code"), verbatimTextOutput("b_code_block"),
                                              tags$button("Copy code", class = "btn btn-default uf-copy",
                                                          onclick = "DEcopy('b_code_block', this)"),
                                              downloadButton("b_dl_code", "Download .R", class = "btn-default uf-dl")),
                                     tabPanel("Wald ANOVA", verbatimTextOutput("b_anova_table")),
                                     tabPanel("DHARMa residuals",
                                              helpText("For 0/1 data the DHARMa QQ/residual plot is naturally ",
                                                       "grainy; focus on the formal tests below rather than the ",
                                                       "visual pattern."),
                                              plotOutput("b_dharma_plot", height = "420px"),
                                              downloadButton("b_dl_dharma", "Download PNG", class = "btn-default uf-dl"),
                                              h4("Overdispersion test"), verbatimTextOutput("b_disp_test"),
                                              h4("Outlier test"), verbatimTextOutput("b_outlier_test"),
                                              code_panel("b_code_dharma", "DHARMa code")),
                                     tabPanel("EMMeans & post-hoc",
                                              helpText("Shown on the probability (response) scale by default."),
                                              h4("Estimated marginal means (probability scale)"), DT::dataTableOutput("b_emm_table"),
                                              downloadButton("b_dl_emm", "Download CSV", class = "btn-default uf-dl"),
                                              h4("Pairwise comparisons (odds-ratio / link scale)"), DT::dataTableOutput("b_pairs_table"),
                                              downloadButton("b_dl_pairs", "Download CSV", class = "btn-default uf-dl"),
                                              h4("Compact letter display"), DT::dataTableOutput("b_cld_table"),
                                              downloadButton("b_dl_cld", "Download CSV", class = "btn-default uf-dl"),
                                              h4("EMMeans plot"), plotOutput("b_emm_plot", height = "460px"),
                                              downloadButton("b_dl_emmplot", "Download PNG", class = "btn-default uf-dl"),
                                              code_panel("b_code_emm", "EMMeans code"))
                         )
               )
             )
    )
  )
)

# ===========================================================================
#  Server
# ===========================================================================
server <- function(input, output, session) {
  
  # -- Shared raw data --------------------------------------------------------
  raw_data <- reactive({
    if (input$data_source == "example") {
      df <- make_example_data()
    } else {
      req(input$file)
      df <- tryCatch(read.csv(input$file$datapath, header = input$header, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      validate(need(!is.null(df) && ncol(df) > 0, "Could not read that file as CSV."))
    }
    df[] <- lapply(df, function(x) {
      if (is.character(x)) factor(x) else if (is.ordered(x)) factor(x, ordered = FALSE) else x
    })
    df
  })
  
  # Build one self-contained "module" of reactive logic per prefix ("g" or "b")
  # so the two tabs are genuinely independent (separate state, separate fit).
  make_module <- function(prefix, binary = FALSE) {
    P <- function(id) paste0(prefix, "_", id)
    rv <- reactiveValues(reset = 0, fit = NULL, message = NULL, combos = list())
    
    dataset <- reactive({
      df <- raw_data(); rv$reset
      ff <- input[[P("force_factor")]]
      if (length(ff)) for (v in ff) if (v %in% names(df)) df[[v]] <- factor(df[[v]])
      for (cm in rv$combos) {
        parts <- lapply(cm$vars, function(v) as.character(df[[v]]))
        df[[cm$name]] <- factor(do.call(paste, c(parts, sep = cm$sep)))
      }
      df
    })
    
    output[[P("force_factor_ui")]] <- renderUI({
      df <- dataset(); req(df)
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      selectizeInput(P("force_factor"), "Treat these numeric columns as factors",
                     choices = num_cols, multiple = TRUE)
    })
    
    output[[P("combo_vars_ui")]] <- renderUI({
      df <- dataset(); req(df)
      selectizeInput(P("combo_vars"), "Columns to combine (choose 2+)",
                     choices = names(df), multiple = TRUE)
    })
    observeEvent(input[[P("combo_add")]], {
      vars <- input[[P("combo_vars")]]; sep <- input[[P("combo_sep")]] %||% "."
      if (length(vars) < 2) {
        rv$message <- list(type = "warning", text = "Pick 2 or more columns to combine.")
        return()
      }
      nm <- paste(vars, collapse = sep)
      rv$combos <- c(rv$combos, list(list(name = nm, vars = vars, sep = sep)))
      rv$message <- list(type = "info", text = sprintf("Created variable '%s'.", nm))
    })
    output[[P("combo_list_ui")]] <- renderUI({
      if (!length(rv$combos)) return(helpText("No combined variables yet."))
      helpText("Created: ", paste(vapply(rv$combos, `[[`, character(1), "name"), collapse = ", "))
    })
    
    output[[P("response_ui")]] <- renderUI({
      df <- dataset(); req(df); rv$reset
      if (binary) {
        cand <- names(df)[vapply(df, function(x)
          (is.numeric(x) && all(x %in% c(0, 1, NA))) ||
            (is.factor(x) && nlevels(x) == 2), logical(1))]
        validate(need(length(cand) > 0, "No 0/1 or two-level factor column found for a binary response."))
        selectInput(P("response"), "Response variable (0/1 or 2-level factor)", choices = cand)
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
        validate(need(length(num_cols) > 0, "No numeric columns available."))
        selectInput(P("response"), "Response variable", choices = num_cols)
      }
    })
    if (binary) output$b_level_note <- renderUI({
      req(identical(input$b_resp_type, "binary"))
      df <- dataset(); r <- input$b_response; req(r)
      if (is.factor(df[[r]])) helpText(sprintf("Modeling P(%s = \"%s\").", r, levels(df[[r]])[2]))
    })
    
    output[[P("fixed_ui")]] <- renderUI({
      df <- dataset(); req(df); rv$reset
      selectizeInput(P("fixed"), sprintf("Fixed effects (up to %d)", MAX_FIXED),
                     choices = names(df), multiple = TRUE,
                     options = list(maxItems = MAX_FIXED, placeholder = "select 1-3 variables"))
    })
    output[[P("random_ui")]] <- renderUI({
      df <- dataset(); req(df); rv$reset
      selectInput(P("random"), "Random effects (grouping factors)", choices = names(df), multiple = TRUE)
    })
    output[[P("slope_ui")]] <- renderUI({
      fx <- input[[P("fixed")]]
      selectInput(P("ranslope"), "Random slope variable", choices = c("(intercept only)" = "", fx), selected = "")
    })
    output[[P("slope_group_ui")]] <- renderUI({
      rs <- input[[P("ranslope")]]; rnd <- input[[P("random")]]
      if (is.null(rs) || !nzchar(rs) || !length(rnd)) return(NULL)
      selectInput(P("slope_group"), "Apply slope to which grouping factor(s)?",
                  choices = rnd, selected = rnd[1], multiple = TRUE)
    })
    
    if (!binary) {
      # Live note for whatever family is currently selected, plus an automatic
      # check of the chosen response's actual range against that family's
      # valid domain (e.g. flags negative values under Gamma/log, or values
      # outside [0,1] under Beta/Ordered-beta) BEFORE the user hits "Fit model".
      output$g_family_note <- renderUI({
        req(input$g_family)
        note <- FAMILY_DOMAIN_NOTE[[input$g_family]] %||% ""
        df <- dataset(); resp <- input$g_response
        warn <- NULL
        if (!is.null(resp) && resp %in% names(df) && is.numeric(df[[resp]])) {
          x <- df[[resp]][is.finite(df[[resp]])]
          out_of_range <- switch(input$g_family,
                                 "gaussian_log" = ,
                                 "Gamma_log"    = any(x <= 0),
                                 "poisson_log"  = ,
                                 "nbinom2_log"  = ,
                                 "nbinom1_log"  = any(x < 0 | x != round(x)),
                                 "tweedie_log"  = any(x < 0),
                                 "beta_logit"   = any(x <= 0 | x >= 1),
                                 "ordbeta_logit" = any(x < 0 | x > 1),
                                 FALSE)
          if (isTRUE(out_of_range))
            warn <- sprintf("\u26a0 '%s' contains values outside this family's valid range \u2014 check the note above before fitting.", resp)
        }
        tagList(helpText(note),
                if (!is.null(warn)) div(class = "alert alert-warning", style = "padding:6px 10px;", warn))
      })
      output$g_zi_vars_ui <- renderUI({
        req(input$g_zi_on)
        df <- dataset(); fx <- input$g_fixed
        selectizeInput("g_zi_vars", "Zero-inflation predictors (blank = intercept only)",
                       choices = names(df), selected = NULL, multiple = TRUE)
      })
      output$g_disp_vars_ui <- renderUI({
        df <- dataset()
        selectizeInput("g_disp_vars", "Dispersion predictors (0, 1, or 2)",
                       choices = names(df), multiple = TRUE,
                       options = list(maxItems = 2, placeholder = "none = constant dispersion"))
      })
    } else {
      output$b_zi_vars_ui <- renderUI({
        req(input$b_zi_on)
        df <- dataset()
        selectizeInput("b_zi_vars", "Zero-inflation predictors (blank = intercept only)",
                       choices = names(df), multiple = TRUE)
      })
      
      # Grouped-binomial support: user picks a successes column and a total-
      # trials column; failures = total - successes are computed automatically
      # (as an arithmetic expression inside the model formula, not a stored
      # column), so no separate "failures" column is ever required.
      int_like <- function(x) is.numeric(x) && all(is.na(x) | (x >= 0 & x == round(x)))
      
      output$b_success_ui <- renderUI({
        df <- dataset(); req(df)
        cand <- names(df)[vapply(df, int_like, logical(1))]
        validate(need(length(cand) > 0,
                      "No non-negative whole-number columns found to use as a successes count."))
        selectInput("b_success", "Successes column (count)", choices = cand)
      })
      output$b_total_ui <- renderUI({
        df <- dataset(); req(df)
        cand <- names(df)[vapply(df, int_like, logical(1))]
        choices <- setdiff(cand, input$b_success %||% "")
        validate(need(length(choices) > 0,
                      "Need another whole-number column, distinct from the successes column, for total trials."))
        selectInput("b_total", "Total trials column (successes + failures)", choices = choices)
      })
      output$b_grouped_note <- renderUI({
        df <- dataset(); succ <- input$b_success; tot <- input$b_total
        req(succ, tot)
        failures <- df[[tot]] - df[[succ]]
        bad <- any(failures < 0, na.rm = TRUE)
        tagList(
          helpText(sprintf("Failures are computed automatically as %s \u2212 %s. The model's response becomes cbind(%s, %s \u2212 %s).",
                            tot, succ, succ, tot, succ)),
          if (isTRUE(bad))
            div(class = "alert alert-warning", style = "padding:6px 10px;",
                sprintf("\u26a0 '%s' is greater than '%s' in at least one row, which would give negative failures. Check the data before fitting.", succ, tot))
        )
      })
    }
    
    output[[P("emm_ui")]] <- renderUI({
      df <- dataset(); fx <- input[[P("fixed")]]
      cat_fx <- if (length(fx)) fx[is_categorical(df, fx)] else character(0)
      selectizeInput(P("emmvars"), "EMMeans factor(s)", choices = cat_fx, multiple = TRUE,
                     options = list(maxItems = 3))
    })
    output[[P("emm_by_ui")]] <- renderUI({
      ev <- input[[P("emmvars")]]
      if (length(ev) < 2) return(NULL)
      selectInput(P("emm_by"), "Compare within (optional, simple effects)",
                  choices = c("(none)" = "", ev), selected = "")
    })
    
    output[[P("data_head")]] <- DT::renderDataTable({
      df <- dataset(); req(df)
      DT::datatable(head(df, 10), options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
    })
    output[[P("data_str")]] <- renderPrint({ df <- dataset(); req(df); str(df) })
    
    output[[P("message_box")]] <- renderUI({
      m <- rv$message; req(m)
      cls <- switch(m$type, warning = "alert-warning", danger = "alert-danger", "alert-info")
      div(class = paste("alert", cls), m$text)
    })
    
    observeEvent(input[[P("reset")]], {
      rv$reset <- rv$reset + 1; rv$fit <- NULL; rv$combos <- list()
      rv$message <- list(type = "info", text = "Reset. Choose variables and fit again.")
    })
    
    # -- Fit -------------------------------------------------------------------
    observeEvent(input[[P("run")]], {
      df  <- dataset()
      fx <- input[[P("fixed")]]; rnd <- input[[P("random")]]
      
      resp_expr <- NULL; resp_label <- NULL; quote_resp <- TRUE
      
      if (binary && identical(input$b_resp_type, "grouped")) {
        succ <- input$b_success; tot <- input$b_total
        validate(need(!is.null(succ) && !is.null(tot) && nzchar(succ) && nzchar(tot),
                      "Choose both a successes column and a total-trials column."))
        if (identical(succ, tot)) {
          rv$message <- list(type = "warning", text = "Successes and total-trials columns must be different.")
          return()
        }
        bad <- df[[succ]] > df[[tot]] | df[[succ]] < 0 | df[[tot]] < 0
        if (any(bad, na.rm = TRUE)) {
          rv$message <- list(type = "danger",
                             text = sprintf("'%s' exceeds '%s' (or a negative value is present) in %d row(s) \u2014 fix the data before fitting.",
                                            succ, tot, sum(bad, na.rm = TRUE)))
          return()
        }
        resp_expr  <- sprintf("cbind(%s, %s - %s)", bq(succ), bq(tot), bq(succ))
        resp_label <- sprintf("%s successes / %s total", succ, tot)
        quote_resp <- FALSE
      } else {
        resp <- input[[P("response")]]
        validate(need(!is.null(resp), "Choose a response variable."))
        if (binary) df[[resp]] <- if (is.factor(df[[resp]])) as.integer(df[[resp]]) - 1L else df[[resp]]
        resp_expr  <- resp
        resp_label <- resp
        quote_resp <- TRUE
      }
      
      if (length(fx) == 0 && length(rnd) == 0) {
        rv$message <- list(type = "warning", text = "Choose at least one fixed or random effect."); return()
      }
      if (length(rnd) == 0) {
        rv$message <- list(type = "warning",
                           text = "No random effect selected \u2014 glmmTMB will fit an ordinary GLM (no grouping structure).")
      }
      
      main_f <- build_formula_string(resp_expr, fx, rnd, input[[P("interactions")]],
                                     input[[P("ranslope")]], input[[P("slope_group")]],
                                     quote_response = quote_resp)
      zi_on  <- if (binary) isTRUE(input$b_zi_on) else isTRUE(input$g_zi_on)
      zi_f   <- build_side_formula(if (binary) input$b_zi_vars else input$g_zi_vars, zi_on)
      disp_f <- if (binary) "~1" else build_side_formula(input$g_disp_vars, TRUE)
      
      fam <- if (binary) binomial(link = input$b_link) else family_from_key(input$g_family)
      
      mod <- tryCatch(
        glmmTMB::glmmTMB(formula = stats::as.formula(main_f),
                         ziformula = stats::as.formula(zi_f),
                         dispformula = stats::as.formula(disp_f),
                         family = fam, data = df),
        error = function(e) e)
      if (inherits(mod, "error")) {
        rv$message <- list(type = "danger", text = paste("Model failed to fit:", conditionMessage(mod)))
        rv$fit <- NULL; return()
      }
      rv$fit <- list(mod = mod, response = resp_label, fixed = fx, random = rnd,
                     main_f = main_f, zi_f = zi_f, disp_f = disp_f,
                     family_key = if (binary) input$b_link else input$g_family, data = df,
                     emmvars = input[[P("emmvars")]], emm_by = input[[P("emm_by")]],
                     adjust = input[[P("adjust")]], conf = input[[P("conf")]])
      rv$message <- list(type = "info", text = "Model fitted.")
    })
    
    output[[P("formula_txt")]] <- renderPrint({
      req(rv$fit)
      cat("Conditional: ", rv$fit$main_f, "\n")
      cat("ziformula:   ", rv$fit$zi_f, "\n")
      cat("dispformula: ", rv$fit$disp_f, "\n")
      fam_text <- if (binary) sprintf("binomial(link='%s')", rv$fit$family_key)
                  else (FAMILY_CODE_TEXT[[rv$fit$family_key]] %||% rv$fit$family_key)
      cat("family:      ", fam_text, "\n")
    })
    output[[P("model_summary")]] <- renderPrint({ req(rv$fit); print(summary(rv$fit$mod)) })
    output[[P("fit_table")]] <- DT::renderDataTable({
      req(rv$fit); m <- rv$fit$mod
      s <- summary(m)
      tab <- data.frame(AIC = round(AIC(m), 2), BIC = round(BIC(m), 2),
                        logLik = round(as.numeric(logLik(m)), 2),
                        df.resid = stats::df.residual(m))
      DT::datatable(tab, options = list(dom = "t"), rownames = FALSE)
    })
    
    build_code <- function() {
      req(rv$fit)
      fam_text <- if (binary) sprintf("binomial(link = '%s')", rv$fit$family_key)
                  else (FAMILY_CODE_TEXT[[rv$fit$family_key]] %||% rv$fit$family_key)
      sprintf(paste0(
        "library(glmmTMB); library(DHARMa); library(emmeans)\n\n",
        "mod <- glmmTMB(\n  %s,\n  ziformula = %s,\n  dispformula = %s,\n  family = %s,\n  data = your_data)\n\n",
        "summary(mod)\n",
        "sim <- simulateResiduals(mod, n = 250, plot = TRUE)\n",
        "testDispersion(sim)\n%s",
        "car::Anova(mod, type = 3)\n"),
        rv$fit$main_f, rv$fit$zi_f, rv$fit$disp_f, fam_text,
        if (!binary) "testZeroInflation(sim)\n" else "")
    }
    output[[P("code_block")]] <- renderPrint(cat(build_code()))
    output[[P("dl_code")]] <- downloadHandler(
      filename = function() "glmmTMB_model_code.R",
      content = function(file) writeLines(build_code(), file))
    
    # -- Wald ANOVA --------------------------------------------------------
    output[[P("anova_table")]] <- renderPrint({
      req(rv$fit)
      if (!HAS_CAR) { cat("Package 'car' not installed \u2014 install it for Type III Wald tests."); return() }
      res <- tryCatch(car::Anova(rv$fit$mod, type = 3), error = function(e) e)
      if (inherits(res, "error")) cat("ANOVA failed:", conditionMessage(res)) else print(res)
    })
    
    # -- DHARMa --------------------------------------------------------------
    dharma_sim <- reactive({ req(rv$fit); DHARMa::simulateResiduals(rv$fit$mod, n = 250) })
    draw_dharma <- function() plot(dharma_sim())
    output[[P("dharma_plot")]] <- renderPlot({ draw_dharma() })
    output[[P("dl_dharma")]] <- downloadHandler(
      filename = function() "dharma_residuals.png",
      content = function(file) { png(file, width = 900, height = 600); draw_dharma(); dev.off() })
    # Two independent overdispersion checks, shown together: (1) DHARMa's
    # simulation-based test (the primary, requested method), and (2) a classic
    # Pearson chi-square / residual-df ratio as a quick cross-check. Ratios
    # near 1 support the chosen family; well above 1 suggests overdispersion
    # (try nbinom2/nbinom1 instead of Poisson, or add a dispformula predictor);
    # well below 1 suggests underdispersion.
    output[[P("disp_test")]] <- renderPrint({
      req(rv$fit)
      ratio <- tryCatch({
        pr <- stats::residuals(rv$fit$mod, type = "pearson")
        rdf <- stats::df.residual(rv$fit$mod)
        sum(pr^2, na.rm = TRUE) / rdf
      }, error = function(e) NA)
      cat("Pearson chi-sq / residual df ratio:", if (is.na(ratio)) "unavailable" else round(ratio, 3), "\n")
      cat("  (~1 = fine; >>1 = overdispersion; <<1 = underdispersion)\n\n")
      cat("DHARMa simulation-based dispersion test:\n")
      print(tryCatch(DHARMa::testDispersion(dharma_sim()), error = function(e) e))
    })
    if (!binary) output$g_zi_test <- renderPrint({
      req(rv$fit)
      print(tryCatch(DHARMa::testZeroInflation(dharma_sim()), error = function(e) e))
    })
    output[[P("outlier_test")]] <- renderPrint({
      req(rv$fit)
      print(tryCatch(DHARMa::testOutliers(dharma_sim()), error = function(e) e))
    })
    
    # -- EMMeans / post-hoc ---------------------------------------------------
    emm_result <- reactive({
      req(rv$fit); ev <- rv$fit$emmvars
      validate(need(length(ev) >= 1, "Select at least one EMMeans factor to see post-hoc results."))
      by <- if (nzchar(rv$fit$emm_by %||% "")) rv$fit$emm_by else NULL
      spec <- if (!is.null(by))
        stats::as.formula(paste("~", paste(bq(setdiff(ev, by)), collapse = " * "), "|", bq(by)))
      else stats::as.formula(paste("~", paste(bq(ev), collapse = " * ")))
      emm <- tryCatch(emmeans::emmeans(rv$fit$mod, spec, type = "response", level = rv$fit$conf),
                      error = function(e) e)
      validate(need(!inherits(emm, "error"), paste("EMMeans failed:", if (inherits(emm, "error")) conditionMessage(emm) else "")))
      pairs_tab <- tryCatch(as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = rv$fit$adjust)),
                            error = function(e) e)
      cld_tab <- tryCatch(as.data.frame(multcomp::cld(emm, adjust = rv$fit$adjust, Letters = letters)),
                          error = function(e) e)
      list(emm = as.data.frame(emm), pairs = pairs_tab, cld = cld_tab, ev = ev, by = by)
    })
    output[[P("emm_table")]] <- DT::renderDataTable({
      DT::datatable(round_df(emm_result()$emm), options = list(dom = "tp"), rownames = FALSE) })
    output[[P("pairs_table")]] <- DT::renderDataTable({
      p <- emm_result()$pairs
      validate(need(!inherits(p, "error"), "Pairwise comparisons failed."))
      DT::datatable(round_df(p), options = list(dom = "tp"), rownames = FALSE) })
    output[[P("cld_table")]] <- DT::renderDataTable({
      c <- emm_result()$cld
      validate(need(!inherits(c, "error"), "Compact letter display failed."))
      DT::datatable(round_df(c), options = list(dom = "tp"), rownames = FALSE) })
    output[[P("dl_emm")]] <- downloadHandler(filename = function() "emmeans.csv",
                                             content = function(file) write.csv(emm_result()$emm, file, row.names = FALSE))
    output[[P("dl_pairs")]] <- downloadHandler(filename = function() "pairwise.csv",
                                               content = function(file) write.csv(emm_result()$pairs, file, row.names = FALSE))
    output[[P("dl_cld")]] <- downloadHandler(filename = function() "cld.csv",
                                             content = function(file) write.csv(emm_result()$cld, file, row.names = FALSE))
    
    build_emm_plot <- function() {
      er <- emm_result(); d <- er$cld; ev <- er$ev
      valcol <- if ("response" %in% names(d)) "response" else "emmean"
      lcol <- intersect(c("asymp.LCL", "lower.CL"), names(d))[1]
      ucol <- intersect(c("asymp.UCL", "upper.CL"), names(d))[1]
      x <- ev[1]; colr <- if (length(ev) >= 2) ev[2] else NA
      d[[x]] <- factor(d[[x]])
      p <- if (is.na(colr)) ggplot(d, aes(x = .data[[x]], y = .data[[valcol]])) + geom_point(size = 3, colour = UF_BLUE)
      else { d[[colr]] <- factor(d[[colr]])
      ggplot(d, aes(x = .data[[x]], y = .data[[valcol]], colour = .data[[colr]])) +
        geom_point(size = 3, position = position_dodge(0.5)) + scale_colour_manual(values = rep(UF_COLORS, 5)) }
      if (!is.na(lcol) && !is.na(ucol)) p <- p + geom_errorbar(aes(ymin = .data[[lcol]], ymax = .data[[ucol]]), width = 0.15)
      if (".group" %in% names(d)) p <- p + geom_text(aes(label = .group), vjust = -0.8, fontface = "bold")
      p + labs(x = x, y = "Estimated marginal mean (response scale)", title = "EMMeans \u00b1 CI") +
        theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", colour = UF_BLUE))
    }
    output[[P("emm_plot")]] <- renderPlot({ print(build_emm_plot()) })
    output[[P("dl_emmplot")]] <- downloadHandler(filename = function() "emmeans_plot.png",
                                                 content = function(file) ggsave(file, build_emm_plot(), width = 9, height = 6, dpi = 150))
    
    build_dharma_code <- function() paste0(
      "# quick manual overdispersion check (Pearson chi-sq / residual df; ~1 is fine)\n",
      "sum(residuals(mod, type = 'pearson')^2) / df.residual(mod)\n\n",
      "sim <- simulateResiduals(mod, n = 250, plot = TRUE)\n",
      "testDispersion(sim)\n", if (!binary) "testZeroInflation(sim)\n" else "", "testOutliers(sim)\n")
    output[[P("code_dharma")]] <- renderPrint(cat(build_dharma_code()))
    output[[P("code_emm")]] <- renderPrint(cat(
      "emm <- emmeans(mod, ~ your_factor, type = 'response')\n",
      "pairs(emm, adjust = 'tukey')\n", "multcomp::cld(emm, adjust = 'tukey', Letters = letters)\n"))
    
    invisible(NULL)
  }
  
  make_module("g", binary = FALSE)
  make_module("b", binary = TRUE)
}

shinyApp(ui, server)
