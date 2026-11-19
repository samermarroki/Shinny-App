library(shiny)
library(bslib)
library(tidyverse)
library(tsibble)
library(fable)
library(feasts)
library(fabletools)
library(lubridate)

# ---- Load data (project-relative path) --------------------------------
raw <- readr::read_csv("AustralianWines.csv",
                       col_types = cols(.default = col_character()),
                       na = c("", "NA", "*", "-"),
                       show_col_types = FALSE)

if ("Varietal" %in% names(raw)) {
  long <- raw
} else {
  month_col <- names(raw)[1]
  long <- raw |>
    pivot_longer(-all_of(month_col), names_to = "Varietal", values_to = "Sales") |>
    rename(Month = all_of(month_col))
}

parse_month <- function(x) {
  d <- suppressWarnings(lubridate::parse_date_time(
    as.character(x), orders = c("my", "bY", "b-y", "Y-m", "Ymd")))
  yearmonth(as.Date(d))
}

wines <- long |>
  mutate(Month = parse_month(Month),
         Sales = suppressWarnings(as.numeric(Sales))) |>
  filter(!is.na(Month)) |>
  as_tsibble(index = Month, key = Varietal) |>
  fill_gaps() |>
  group_by_key() |>
  mutate(Sales = zoo::na.approx(Sales, na.rm = FALSE)) |>
  ungroup() |>
  filter(!is.na(Sales))

all_varietals <- sort(unique(wines$Varietal))
date_min <- as.Date(min(wines$Month))
date_max <- as.Date(max(wines$Month))

# ---- UI ---------------------------------------------------------------
ui <- page_sidebar(
  title = "Australian Wine Sales: Model Comparison",
  sidebar = sidebar(
    width = 320,
    selectizeInput("varietals", "Varietal(s)", choices = all_varietals,
                   selected = all_varietals[1:2], multiple = TRUE),
    dateRangeInput("dates", "Date range", start = date_min, end = date_max,
                   min = date_min, max = date_max, format = "yyyy-mm"),
    sliderInput("val_months", "Validation window (months)", 6, 36, 12, step = 1),
    sliderInput("horizon", "Forecast horizon beyond validation (months)", 6, 36, 12, step = 1),
    checkboxGroupInput("models", "Models",
                       choices = c("TSLM", "ETS", "ARIMA", "Ensemble (mean of all)"),
                       selected = c("TSLM", "ETS", "ARIMA")),
    checkboxInput("log_tx", "Log-transform sales (stabilizes seasonal variance)", FALSE),
    checkboxInput("free_y", "Free y-axis per varietal", TRUE),
    sliderInput("level", "Prediction interval (%)", 50, 99, 80),
    helpText("Models are fit on data before the validation window, evaluated on the validation window, then projected forward.")
  ),
  navset_card_tab(
    nav_panel("Overview", plotOutput("overview", height = 450)),
    nav_panel("Decomposition", plotOutput("decomp", height = 500)),
    nav_panel("Forecasts", plotOutput("fc_plot", height = 550)),
    nav_panel("Model specs", tableOutput("specs")),
    nav_panel("Accuracy", tableOutput("acc")),
    nav_panel("About", markdown(
      "**Data:** AustralianWines.csv (monthly sales by varietal, 1980 to 1995). Missing months are linearly interpolated.

**Method:** TSLM, ETS, and ARIMA models are auto-specified with the fable package and fit on a training window. They are evaluated on a held-out validation window, then forecast forward.

**Packages:** shiny, bslib, tidyverse, tsibble, fable, feasts, fabletools, lubridate, zoo.

**Acknowledgments:** App code was drafted with assistance from Claude (Anthropic) and reviewed, tested, and deployed by the author.")),
    nav_panel("Debug", verbatimTextOutput("debug"))
  )
)

# ---- Server -----------------------------------------------------------
server <- function(input, output, session) {
  
  facet_scales <- reactive(if (isTRUE(input$free_y)) "free_y" else "fixed")
  
  dat <- reactive({
    req(length(input$varietals) >= 1, input$dates)
    wines |>
      filter(Varietal %in% input$varietals,
             Month >= yearmonth(input$dates[1]),
             Month <= yearmonth(input$dates[2]))
  })
  
  split_info <- reactive({
    d <- dat()
    validate(need(nrow(d) > 0, "No data in this range."))
    cutoff <- max(d$Month) - input$val_months
    train <- filter(d, Month <= cutoff)
    n_train <- train |> as_tibble() |> count(Varietal) |> pull(n) |> min()
    validate(need(n_train >= 36,
                  "Training window too short. Widen the date range or shrink the validation window."))
    list(train = train, valid = filter(d, Month > cutoff), cutoff = cutoff)
  })
  
  fits <- reactive({
    req(length(input$models) > 0)
    tr <- split_info()$train
    
    withProgress(message = "Fitting models...", value = 0.5, {
      f <- if (isTRUE(input$log_tx)) {
        tr |> model(
          TSLM  = TSLM(log(Sales) ~ trend() + season()),
          ETS   = ETS(log(Sales)),
          ARIMA = ARIMA(log(Sales))
        )
      } else {
        tr |> model(
          TSLM  = TSLM(Sales ~ trend() + season()),
          ETS   = ETS(Sales),
          ARIMA = ARIMA(Sales)
        )
      }
      f <- f |> mutate(Ensemble = (TSLM + ETS + ARIMA) / 3)
      keep <- intersect(c("TSLM", "ETS", "ARIMA", "Ensemble"),
                        sub(" \\(mean of all\\)", "", input$models))
      f |> select(Varietal, all_of(keep))
    })
  })
  
  fcst <- reactive({
    fits() |> forecast(h = input$val_months + input$horizon)
  })
  
  output$overview <- renderPlot({
    dat() |>
      as_tibble() |>
      mutate(Date = as.Date(Month)) |>
      ggplot(aes(Date, Sales)) +
      geom_line(colour = "steelblue") +
      facet_wrap(~Varietal, scales = facet_scales(), ncol = 1) +
      labs(y = "Sales (thousand litres)", x = NULL) +
      theme_minimal(base_size = 14)
  })
  
  output$decomp <- renderPlot({
    req(length(input$varietals) >= 1)
    v <- input$varietals[1]
    
    comp <- dat() |>
      filter(Varietal == v) |>
      model(STL(Sales ~ season(window = "periodic"))) |>
      components() |>
      as_tibble() |>
      transmute(Date = as.Date(Month),
                Data = Sales,
                Trend = trend,
                Seasonal = season_year,
                Remainder = remainder) |>
      pivot_longer(-Date, names_to = "Component", values_to = "Value") |>
      mutate(Component = factor(Component,
                                levels = c("Data", "Trend", "Seasonal", "Remainder")))
    
    ggplot(comp, aes(Date, Value)) +
      geom_line(colour = "steelblue") +
      facet_wrap(~Component, ncol = 1, scales = "free_y") +
      labs(title = paste("STL decomposition:", v), x = NULL, y = NULL) +
      theme_minimal(base_size = 13)
  })
  
  output$fc_plot <- renderPlot({
    si <- split_info()
    L <- input$level / 100
    
    fc <- fcst() |>
      mutate(lo = quantile(Sales, (1 - L) / 2),
             hi = quantile(Sales, 1 - (1 - L) / 2),
             pt = .mean) |>
      as_tibble() |>
      mutate(Date = as.Date(Month))
    
    hist <- dat() |> as_tibble() |> mutate(Date = as.Date(Month))
    
    ggplot() +
      geom_line(data = hist, aes(Date, Sales), colour = "grey30") +
      geom_ribbon(data = fc, aes(Date, ymin = lo, ymax = hi, fill = .model), alpha = 0.2) +
      geom_line(data = fc, aes(Date, pt, colour = .model), linewidth = 0.8) +
      geom_vline(xintercept = as.Date(si$cutoff + 1), linetype = "dashed") +
      facet_wrap(~Varietal, scales = facet_scales(), ncol = 1) +
      labs(y = "Sales (thousand litres)", x = NULL, colour = "Model", fill = "Model",
           caption = paste0("Dashed line = start of validation window. Shaded = ",
                            input$level, "% prediction interval.")) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
  
  output$specs <- renderTable({
    fits() |>
      as_tibble() |>
      pivot_longer(-Varietal, names_to = "Model", values_to = "fit") |>
      mutate(Specification = map_chr(fit, function(m) {
        txt <- tryCatch(format(m), error = function(e) NA_character_)
        gsub("^<|>$", "", txt)
      })) |>
      select(Varietal, Model, Specification)
  })
  
  output$acc <- renderTable({
    si <- split_info()
    train_acc <- accuracy(fits()) |> mutate(Window = "Train")
    valid_acc <- fcst() |>
      filter(Month <= si$cutoff + input$val_months) |>
      accuracy(si$valid) |>
      mutate(Window = "Validation")
    bind_rows(train_acc, valid_acc) |>
      select(Varietal, Model = .model, Window, RMSE, MAE, MAPE) |>
      arrange(Varietal, Window, RMSE) |>
      mutate(across(c(RMSE, MAE, MAPE), ~ round(.x, 2)))
  })
  
  # ---- Debug: shows why ARIMA fails in a given environment ------------
  output$debug <- renderPrint({
    tr <- split_info()$train |> filter(Varietal == input$varietals[1])
    
    cat("R:", R.version.string, "\n\nPackage versions:\n")
    print(sapply(c("fable", "fabletools", "feasts", "tsibble", "distributional"),
                 function(p) as.character(packageVersion(p))))
    
    cat("\nVarietal:", input$varietals[1],
        "| Training rows:", nrow(tr),
        "| NAs in Sales:", sum(is.na(tr$Sales)), "\n")
    
    try_fit <- function(label, expr) {
      msgs <- character()
      res <- withCallingHandlers(
        tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
        warning = function(w) {
          msgs <<- c(msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      out <- if (is.character(res)) res else tryCatch(format(res$A[[1]]),
                                                      error = function(e) "could not format")
      cat("\n", label, "\n  result: ", out, "\n", sep = "")
      if (length(msgs)) cat("  warnings: ", paste(msgs, collapse = " | "), "\n", sep = "")
    }
    
    try_fit("1. Auto ARIMA (what the app uses)",
            tr |> model(A = ARIMA(Sales)))
    try_fit("2. Auto ARIMA, stepwise = FALSE",
            tr |> model(A = ARIMA(Sales, stepwise = FALSE)))
    try_fit("3. Auto ARIMA, approximation = FALSE",
            tr |> model(A = ARIMA(Sales, approximation = FALSE)))
    try_fit("4. Fixed order (0,1,1)(0,1,1)[12]",
            tr |> model(A = ARIMA(Sales ~ pdq(0,1,1) + PDQ(0,1,1))))
  })
}

shinyApp(ui, server)
