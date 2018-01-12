#R

# AUTHOR CHristian Trachsel 2017


#data shaping

.calc.master.scan <- function(x){
  if("MasterScanNumber" %in% names(x)){
    x <- dplyr::mutate(x, MasterScanNumber = replace(MasterScanNumber, MasterScanNumber == 0, NA))
    return(x)
  } else {
    set1 <- x %>% 
      dplyr::filter(MSOrder == "Ms") %>% 
      dplyr::select(MasterScanNumber = scanNumber, CycleNumber)
    set2 <- dplyr::select(x, scanNumber, CycleNumber)
    res <- dplyr::left_join(set2, set1, by = "CycleNumber") %>% 
      dplyr::mutate(type = x$ScanType) %>% 
      dplyr::mutate(MasterScanNumber = replace(MasterScanNumber, scanNumber == MasterScanNumber, NA)) %>% 
      dplyr::select(MasterScanNumber) %>% 
      dplyr::bind_cols(x, .)
    return(res)
  }
}

.calc.transient <- function(x){
  if("OrbitrapResolution" %in% names(x)){
    x <- x %>% 
      dplyr::rename(FTResolution = OrbitrapResolution)
  } else {
    x <- x
  }
  res <- x %>%
    dplyr::mutate(FTResolution = replace(FTResolution, MassAnalyzer == "MassAnalyzerITMS", NA)) %>% 
    dplyr::mutate(transient = dplyr::case_when(FTResolution == 15000 ~ 32, 
                                               FTResolution == 17500 ~ 64,
                                               FTResolution == 30000 ~ 64,
                                               FTResolution == 50000 ~ 96,
                                               FTResolution == 35000 ~ 128,
                                               FTResolution == 60000 ~ 128,
                                               FTResolution == 70000 ~ 256,
                                               FTResolution == 120000 ~ 256,
                                               FTResolution == 140000 ~ 512,
                                               FTResolution == 240000 ~ 512
                                              )
                 )
  return(res)
}  

#plot functions
.TIC.BasePeak <- function(x){
  df <- x %>% 
    dplyr::filter(grepl("ms ", ScanType)) %>% 
    dplyr::select(StartTime, TIC, BasePeakIntensity) %>% 
    dplyr::rename(Base_Peak = BasePeakIntensity) %>% 
    tidyr::gather(key = "Type", value = "Intensity", TIC, Base_Peak)
  df$Type <- factor(df$Type, levels = c("TIC", "Base_Peak"))
 
   figure <- ggplot(df,aes(x= StartTime, y = Intensity)) +
    geom_line(size = 0.3) +
    facet_wrap(~Type, scales = "free", nrow = 2, ncol = 1) +
    labs(title = "TIC and Base-Peak plot", subtitle = "Plotting TIC intensity and base peak intensity against retention time") +
    labs(x = " Retention Time [min]", y = "Intensity Counts [arb. unit]") +
    scale_x_continuous(breaks = scales::pretty_breaks(15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    theme_light()
  return(figure)
}

.cycle.time <- function(x){
  df <- x %>% 
    dplyr::filter(MSOrder == "Ms") %>% 
    dplyr::select(StartTime) %>% 
    dplyr::mutate(CycleTime = (StartTime - lag(StartTime))*60) %>% 
    dplyr::select(CycleTime, StartTime) %>% 
    na.omit()
  
  figure<- ggplot(df, aes(x = StartTime, y = CycleTime)) + 
    geom_point(shape = ".") +
    geom_line(stat = "smooth", method = "gam", formula = y ~ s(x, bs= "cs"), colour = "deepskyblue3", se = FALSE) +
    labs(title = "Cycle time plot", subtitle = "Caclulated cycle time vs retention time") +
    labs(x = "Retention Time [min]", y = "Cycle Time [sec]") +
    geom_hline(yintercept = quantile(df$CycleTime, 0.95), colour = "red3", linetype = "longdash") +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    theme_light()
  return(figure)
}


.mz.dist <- function(x){
  df <- x %>% 
    dplyr::filter(MSOrder == "Ms2")
  
  figure <- ggplot(df, aes(x = StartTime, y = PrecursorMass)) + 
    geom_point(shape = ".") +
    geom_line(stat = "smooth", method = "gam", formula = y ~ s(x, bs= "cs"), size = 1.1, alpha = 0.6, colour = "deepskyblue3", se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    labs(title = "Retention Time to m/z correlation plot") +
    labs(subtitle = "Plotting retention time against m/z value of all selected precursors") +
    labs(x = "Retention Time", y = "Presursor m/z value") +
    theme_light()
  return(figure)
}

.mz.frequency <- function(x){
  res <-  x %>% 
    dplyr::filter(MSOrder == "Ms2") %>% 
    dplyr::select(ChargeState, PrecursorMass) %>% 
    dplyr::mutate(deconv = round((PrecursorMass -1.00782)*ChargeState, 0))
  
  figure <- ggplot(res, aes(x=deconv, fill = factor(ChargeState), colour = factor(ChargeState))) +
    geom_histogram(binwidth= 100, alpha=.3, position="identity") +
    labs(title = "Precursor mass to charge frequency plot ") +
    labs(subtitle = "Plotting frequency of precursor masses for each charge state") +
    labs(x = "Precursor mass [neutral mass]", y = "Frequency [counts]") +
    labs(fill="Charge State", colour = "Charge State") +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    theme_light()
  return(figure)
}

.charge.states <- function(x){
  res <- x %>% 
    dplyr::filter(MSOrder == "Ms2") %>% 
    dplyr::count(ChargeState) %>% 
    dplyr::rename(Counts = n) %>% 
    dplyr::mutate(percentage = (100/sum(Counts)*Counts))
  xbreaks <- unique(res$ChargeState)
  
  figure <- ggplot(res, aes(x = ChargeState, y = percentage)) +
    geom_bar(stat = "identity", fill = "deepskyblue2") +
    geom_text(aes(label = Counts), vjust=-0.3, size=3.5) +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(breaks = scales::pretty_breaks(15), expand = c(0, 0), limits = c(0, (max(res$percentage))+3)) +
    labs(title = "Charge state plot", subtitle = "Plotting the number of occurances of all selected precursor charge states") +
    labs(x = "Charge States", y = "Percent [%]") +
    theme_light()
  return(figure)
}

.charge.states.quantiles <- function(x){
  res <- x %>% 
    dplyr::filter(MSOrder == "Ms2") %>% 
    dplyr::select(ChargeState) %>% 
    dplyr::mutate(steps = cume_dist(.$ChargeState)) %>% 
    dplyr::group_by(ChargeState) %>% 
    dplyr::do(unique(.))
  zero <- data.frame(ChargeState = 0, steps = 0)
  res <- dplyr::bind_rows(zero, res)
  xbreaks <- res$ChargeState
  
  figure <- ggplot(res, aes(x = ChargeState, y = steps)) +
    geom_step(colour = "deepskyblue2", size = 1.3, alpha = 0.3) +
    geom_point(shape = 18, colour = "black") +
    geom_hline(yintercept = 0.95, colour = "red3", linetype = "longdash") +
    geom_text(aes(label = signif(steps, 3)), vjust = -0.5, hjust = 1) +
    scale_y_continuous(breaks = scales::pretty_breaks(10)) +
    scale_x_continuous(breaks = xbreaks) +
    labs(title = "Cumulative charge state percentage plot") +
    labs(x = "Charge States", y = "Percent [%]") +
    theme_light()
  return(figure)
}

.scan.times <- function(x){
    if("ElapsedScanTimesec" %in% names(x)){
      res <- x %>% 
        dplyr::mutate(ElapsedScanTimesec = ElapsedScanTimesec * 1000)
    } else {
      res <- x %>% 
        dplyr::mutate(ElapsedScanTimesec = (lead(x$StartTime)-x$StartTime)*60000) %>% 
        dplyr::select(StartTime, scanNumber, MSOrder, ElapsedScanTimesec, transient) %>% 
        dplyr::filter(!is.na(.$ElapsedScanTimesec))
    } 
    figure <- ggplot(res, aes(x = StartTime, y = ElapsedScanTimesec))+
      geom_point(shape = ".")+
      facet_grid(MSOrder~., scales = "free")+
      geom_line(stat = "smooth", method = "gam", formula = y~s(x), colour = "deepskyblue3", se = FALSE)+
      labs(title = "Scan time plot", subtitle = "Plotting the elapsed scan time for each individual scan")+
      labs(x = "Retentione Time [min]", y = "Scan Time [ms]")+
      scale_x_continuous(breaks = scales::pretty_breaks((n = 8)))+
      scale_y_continuous(breaks = scales::pretty_breaks((n = 8)))+
      theme_light()+
      geom_hline(data = res, aes(yintercept = transient), colour = "red3")
    return(figure)
}

.injection.times <- function(x){
  if("Max.IonTimems" %in% names(x)){
    maxtimes <- x %>% 
      dplyr::group_by(MSOrder) %>% 
      dplyr::summarise(maxima = max(Max.IonTimems))
  } else {
    maxtimes <- x %>% 
      dplyr::group_by(MSOrder) %>% 
      dplyr::summarise(maxima = max(IonInjectionTimems))
  }
  
  figure <- ggplot(x, aes(x = StartTime, y = IonInjectionTimems)) +
    geom_hline(data = maxtimes, aes(yintercept = maxima), colour = "red3", linetype = "longdash") +
    geom_point(shape = ".") +
    geom_line(stat = "smooth", method = "gam", formula = y~s(x), colour = "deepskyblue3", se = FALSE) +
    facet_grid(MSOrder~., scales = "free") +
    scale_y_continuous(breaks = scales::pretty_breaks((n = 8))) +
    scale_x_continuous(breaks = scales::pretty_breaks((n = 8))) +
    labs(title = "Injection time plot", subtitle = "Plotting injection time against retention time for MS and MSn level") +
    labs(x = "Retentione Time [min]", y = "Injection Time [ms]") +
    theme_light()
  return(figure)
}

.lm.correction <- function(x){
  if("LMCorrectionppm" %in% names(x)){
    res <- x %>% 
      dplyr::filter(MSOrder == "Ms")
  } else {
    res <- x %>% 
      dplyr::filter(MSOrder == "Ms") %>% 
      dplyr::rename(LMCorrectionppm = LMmZCorrectionppm)
  }  
  
  figure <- ggplot(res, aes(x = StartTime , y = LMCorrectionppm)) +
    geom_hline(yintercept = c(-5,5), colour = "red3", linetype = "longdash") +
    geom_line(size = 0.3) +
    geom_line(stat = "smooth", method= "gam", formula = y ~ s(x, bs ="cs"), colour = "deepskyblue3", se = FALSE) +
    labs(title = "Lock mass correction plot", subtitle = "Plotting lock mass correction value over time") +
    labs(x = "Retention Time [min]", y = "Lock Mass Correction [ppm]") +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(8), limits = c(-10,10)) +
    theme_light()
  return(figure)
}



#' Title
#'
#' @param x 
#'
#' @return
#' @import dplyr
#'
#' @examples
.ms2.frequency <- function(x){
  NoMS2 <- x %>% 
    dplyr::filter(MSOrder == "Ms") %>% 
    dplyr::count() %>% 
    dplyr::rename(Counts = n) %>% 
    dplyr::pull()
  res <- x %>% 
    dplyr::filter(MSOrder == "Ms2") %>% 
    dplyr::count(MasterScanNumber) %>% 
    dplyr::count(n) %>% 
    dplyr::rename(NumberOfMS2Scans = n, Counts = nn) %>% 
    rbind(c(0, NoMS2 -sum(.$Counts))) %>% 
    dplyr::arrange(NumberOfMS2Scans) %>% 
    dplyr::mutate(percentage = signif((100/sum(Counts) * Counts), 2))
  xbreaks <- res$NumberOfMS2Scans
  if(max(res$NumberOfMS2Scans >=25)){
    res <- res %>% 
      dplyr::bind_cols(parts = unlist(cut(.$NumberOfMS2Scans, breaks = c(-1,25,50,75,100))))
    levels(res$parts) <- list("1-25" = levels(res$parts)[1], 
                              "26-50" = levels(res$parts)[2], 
                              "51-75" = levels(res$parts)[3], 
                              "76-100" = levels(res$parts)[4]) 
    res <- res %>% 
      dplyr::select(NumberOfMS2Scans, Counts, percentage) %>% 
      dplyr::mutate(parts = as.factor("ALL")) %>% 
      dplyr::bind_rows(res) %>% 
      dplyr::mutate(x_min = case_when(parts == "ALL" ~ 1,
                                      parts == "1-25" ~ 1,
                                      parts == "26-50" ~ 26,
                                      parts == "51-75" ~ 51,
                                      parts == "76-100" ~ 76)) %>% 
      dplyr::mutate(x_max = case_when(parts == "ALL" ~ max(res$NumberOfMS2Scans),
                                      parts == "1-25" ~ 25,
                                      parts == "26-50" ~ 50,
                                      parts == "51-75" ~ 75,
                                      parts == "76-100" ~ 100))
    res$parts <- factor(res$parts, levels = c("ALL", "1-25", "26-50", "51-75", "76-100"))  
    
    figure <- ggplot(res, aes(x=NumberOfMS2Scans, y = percentage)) + 
      geom_bar(stat = "identity", fill = "deepskyblue2") + 
      geom_text(aes(label = Counts), vjust=-0.3, size=3.5) + 
      facet_wrap(~parts, scale = "free", nrow = 5, ncol = 1) +
      ylim(0, max(res$percentage+5)) + 
      geom_blank(aes(x = x_min)) +
      geom_blank(aes(x = x_max))
  } else {
    figure <- ggplot(res, aes(x = NumberOfMS2Scans, y = percentage)) +
      geom_bar(stat = "identity", fill = "deepskyblue2") +
      geom_text(aes(label = Counts), vjust=-0.3, size=3.5)
  }
  figure +
    scale_x_continuous(breaks = xbreaks) +
    labs(title = "Cycle load plot") +
    labs(subtitle = "Plotting the number of MS2 scans associated with each MS1 scan") +
    labs(x = "Number of MS2 associated with an MS1 scan", y = "Percentage [%]") +
    theme_light()
}

.ms.data.points <- function(x){
  binSize <- 15
  binNumber <- ceiling(nrow(x)/binSize)
  binVector <- rep(1:binNumber, each = binSize)
  res <- x %>% 
    dplyr::filter(MSOrder == "Ms") %>% 
    dplyr::select(StartTime) %>% 
    dplyr::mutate(CycleTime = (dplyr::lead(StartTime) - StartTime)*60) %>% 
    dplyr::filter(!is.na(CycleTime)) %>% 
    dplyr::mutate("10sec" = floor(10/CycleTime)) %>% 
    dplyr::mutate("20sec" = floor(20/CycleTime)) %>% 
    dplyr::mutate("30sec" = floor(30/CycleTime)) %>% 
    dplyr::select(StartTime, "10sec", "20sec", "30sec") %>% 
    dplyr::mutate(Bins = binVector[1:nrow(.)]) %>% 
    dplyr::group_by(Bins) %>% 
    dplyr::summarise_all(funs(mean)) %>% 
    dplyr::mutate_at(c("10sec","20sec","30sec"), funs(floor)) %>% 
    tidyr::gather("PeakWidthAtBaseline", "Points", 3:5)
  
  figure <- ggplot(res, aes(x = StartTime, y = Points, colour = PeakWidthAtBaseline)) +
    geom_point(size = 0.3) +
    scale_colour_manual(values = c("red3", "darkorchid3", "deepskyblue3")) +
    scale_y_continuous(breaks = scales::pretty_breaks((n = 20))) + 
    scale_x_continuous(breaks = scales::pretty_breaks((n = 8))) +
    labs(title ="Point Over Chromatographic Peak") + 
    labs(subtitle = "Plotting the number of Ms data points over different preselected chromatographic peak widths") +
    labs(x = "Retention Time", y = "Points over Peak") +
    labs(colour = "Peak width") +
    theme_light()
  return(figure)
}

.ms2.vs.RT <- function(x){
  MS2 <- x %>% 
    dplyr::filter(MSOrder == "Ms2") %>% 
    dplyr::count(MasterScanNumber) %>% 
    dplyr::rename(scanNumber = MasterScanNumber)
  MS <- x %>% 
    dplyr::select(StartTime, scanNumber)
  res <- dplyr::inner_join(MS, MS2, by = "scanNumber")
  
  figure <- ggplot(res, aes(x = StartTime, y = n)) +
    geom_point(shape = ".") +
    geom_line(stat = "smooth", method = "loess", span = 0.2, colour = "deepskyblue3", se = FALSE) +
    scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    coord_cartesian(ylim = c(0, max(res$n)+1)) +
    labs(title = "Time resolved number of Ms2 scans") +
    labs(subtitle = "Plotting the number of Ms2 per Ms1 scan versus retention time") +
    labs(x = "Retention Time [min]", y = "Number of Ms2 per Ms1 [counts]") +
    theme_light()
  return(figure)
}