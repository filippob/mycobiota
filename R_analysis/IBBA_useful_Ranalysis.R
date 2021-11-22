#' title: "Data analysis on mycobiota data of USEFUL project"
#' author: Francesco Vitali
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' ```{r global_options, echo = FALSE, include = FALSE}
#' options(width = 9999)
#' knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
#'                      cache = FALSE, tidy = FALSE, size = "small",
#'                      fig.height = 10, fig.width = 15,fig.align = "center")
#' ```

#' # INTRODUCTION
#' 
#' In this script we will perform classical alpha and beta diversity analysis on the results of mycobiota sequencing 
#' of Taleggio cheese. We will perform analysis on the R1 and R2 separately (i.e. using picking results from non merged reads)
#' then comparing results (i.e. Mantel test or procsutes rotation). Finally, we will explore sthe effect of some experimental metadata
#' 
#' 

#+ echo=FALSE, message = FALSE
