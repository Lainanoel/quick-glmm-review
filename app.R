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
#                       odd with my 0/1 column" failure mode.
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
build_formula_string <- function(response, fixed, random, interactions,
                                 slope = NULL, slope_group = NULL) {
  fixed_part  <- build_fixed_part(fixed, interactions)
  random_part <- build_random_part(random, slope, slope_group)
  paste(bq(response), "~", fixed_part, "+", random_part)
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
  "Beta (logit link)  \u2013 proportions in (0,1)"  = "beta_logit"
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
    gaussian())
}
BINARY_LINK_CHOICES <- c("logit" = "logit", "probit" = "probit",
                         "cloglog" = "cloglog", "cauchit" = "cauchit")

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

  # binary presence/absence -> binomial (own tab)
  peta <- -0.3 + 0.9 * (d$Treatment == "Fertilized") - 0.6 * (d$Treatment == "Grazed") + b
  d$present <- rbinom(n, 1, plogis(peta))

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
          selectInput("g_family", "Family (distribution + link)",
                      choices = FAMILY_CHOICES, selected = "nbinom2_log"),
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
    # TAB 2 : Binary (0/1) GLMM  \u2014 kept fully separate from the general tab
    # =======================================================================
    tabPanel("Binary (0/1) GLMM",
      sidebarLayout(
        sidebarPanel(width = 4,
          helpText(tags$b("Binary outcomes are handled separately."), " The response ",
                   "must be numeric 0/1 or a two-level factor; family is fixed to a ",
                   "Bernoulli/binomial distribution (dispersion is not estimated for ",
                   "true 0/1 data, so no dispersion-model picker is shown here)."),
          setup_controls("b"),
          tags$hr(),
          uiOutput("b_response_ui"),
          uiOutput("b_level_note"),
          selectInput("b_link", "Link function", choices = BINARY_LINK_CHOICES, selected = "logit"),
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
      resp <- input[[P("response")]]; fx <- input[[P("fixed")]]; rnd <- input[[P("random")]]
      validate(need(!is.null(resp), "Choose a response variable."))
      if (length(fx) == 0 && length(rnd) == 0) {
        rv$message <- list(type = "warning", text = "Choose at least one fixed or random effect."); return()
      }
      if (length(rnd) == 0) {
        rv$message <- list(type = "warning",
                           text = "No random effect selected \u2014 glmmTMB will fit an ordinary GLM (no grouping structure).")
      }
      if (binary) df[[resp]] <- if (is.factor(df[[resp]])) as.integer(df[[resp]]) - 1L else df[[resp]]

      main_f <- build_formula_string(resp, fx, rnd, input[[P("interactions")]],
                                     input[[P("ranslope")]], input[[P("slope_group")]])
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
      rv$fit <- list(mod = mod, response = resp, fixed = fx, random = rnd,
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
      cat("family:      ", if (binary) sprintf("binomial(link='%s')", rv$fit$family_key)
                          else rv$fit$family_key, "\n")
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
      sprintf(paste0(
        "library(glmmTMB); library(DHARMa); library(emmeans)\n\n",
        "mod <- glmmTMB(\n  %s,\n  ziformula = %s,\n  dispformula = %s,\n  family = %s,\n  data = your_data)\n\n",
        "summary(mod)\n",
        "sim <- simulateResiduals(mod, n = 250, plot = TRUE)\n",
        "testDispersion(sim)\n%s",
        "car::Anova(mod, type = 3)\n"),
        rv$fit$main_f, rv$fit$zi_f, rv$fit$disp_f,
        if (binary) sprintf("binomial(link = '%s')", rv$fit$family_key) else rv$fit$family_key,
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
    output[[P("disp_test")]] <- renderPrint({
      req(rv$fit)
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
