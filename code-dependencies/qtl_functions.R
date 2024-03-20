### DOCKER VERSION
### Functions for QTL analysis (SARS and peanut allergy)
library(lme4) # for modeling random effects
library(MESS) # for AUC calculation 
library(plyr)
library(shades)
library(scales)
library(extRemes)
library(tidyverse)
source("code-dependencies/lmmultiresponse.R")

#-------------------------------Miscellaneous----------------------------------#
which_chr <- function(cross, marker) {
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o] 
  return(chr)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# has_xchr_perms: indicates whether permutatino object has x-chr-specific permutations  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
has_xchr_perms <- function(perm){
  return('xchr' %in% names(attributes(perm)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# add_attributes: add an attribute to an object  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_attribute <- function(obj, attr_name, attr_value){
  attr(obj, attr_name) <- attr_value
  return(obj)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# logit transformation (for proportions/percentages)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
logit <- function(p){log(p/(1-p))}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Will's truncated-RINT transformation 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
trint <- function(y, theta=0.01){
  p <- theta + (rank(y)-1)*(1-2*theta)/(length(y)-1)
  qnorm(p)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to ensure directory exists before creating files in it
#+++++++++++++++++++++++++++++++++++++++++++++++++++
ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create logger function
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_logger <- function(filename, sep="\n"){
  if(file.exists(filename)){
    file.remove(filename);
  }
  function(...){
    text <- sprintf(...);
    cat(text, file=filename, sep=sep, append=T);
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to calculate number of markers in a cross object 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
calc.num.markers <- function(cross){
  p = 0
  for (c in 1:length(cross$geno)){
    p = p + ncol(cross$geno[[c]]$data)
  }
  return(p)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to add phenotype data to R/qtl data table 
# Input:
#   pheno.name = name of phenotype column in phenotype spreadsheet 
#   new.col = empty column 
#   col.pos = position of new column in R/qtl table 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_pheno <- function(pheno.name, new.col, col.pos){
  rqtl <- rqtl %>% add_column(placeholder.name = new.col, .before = col.pos)
  names(rqtl)[names(rqtl) == "placeholder.name"] <- pheno.name
  for (i in 3:nrow(rqtl)){
    rqtl[[pheno.name]][i] <- ifelse(rqtl$mouse_ID[i] %in% pheno$Geno_ID, 
                                    pheno[[pheno.name]][which(pheno$Geno_ID == rqtl$mouse_ID[i])], NA)
  }
  return(rqtl)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# readGRMBin: R script to read the GRM binary file
# Output:
#     cross: r/qtl cross object updated with new aggregate phenotype 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
readGRMBin <- function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N, grm=grm))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# fill_in_res:  when a phenotype is missing values (e.g. flow data in SARS 
#               experiments), residuals(fit) will be shorter than the original 
#               phenotype vector/column. This function expans the residuals 
#               vector to the original length, filling in the missing spaces 
#               with NAs. 
# Input:
#     pheno: vector of phenotype values (original length)
#     fit: fit with residuals to use in new data vector 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
fill_in_res <- function(pheno, fit){
  has.data <- as.numeric(names(residuals(fit)))
  missing.data <- !seq(1,length(pheno)) %in% has.data # boolean vector (TRUE if missing)  
  
  tmp <- vector(length = length(pheno))
  tmp[missing.data] <- NA
  tmp[has.data] <- residuals(fit)
  
  return(tmp)
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# load_themes: loads ggplot themes 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
load_themes <- function(){
  # To find default colors
  #show_col(hue_pal()(2))
  
  poster_theme <<- theme(axis.title.x = element_text(size = 20), 
                         axis.text.x = element_text(size = 18), 
                         axis.title.y = element_text(size = 20),
                         axis.text.y = element_text(size = 18),
                         legend.title = element_text(size = 20),
                         legend.text = element_text(size = 18),
                         plot.title = element_text(size = 22))
  
  poster_theme2 <<- theme(axis.title.x = element_text(size = 22), 
                          axis.text.x = element_text(size = 20), 
                          axis.title.y = element_text(size = 22),
                          axis.text.y = element_text(size = 20),
                          legend.title = element_text(size = 22),
                          legend.text = element_text(size = 20),
                          plot.title = element_text(size = 24))
  
  # For things needing big text  
  big_theme <<- theme(axis.title = element_text(size = 24), 
                      axis.text = element_text(size = 22),
                      plot.title = element_text(size = 24),
                      legend.title = element_text(size = 22),
                      legend.text = element_text(size = 18))
  
  # For things needing big text but smaller x-axis labels   
  # (like HS and titer plots with x-axis = strain)
  big_theme2 <<- theme(axis.title = element_text(size = 24), 
                       axis.text.y = element_text(size = 22),
                       axis.text.x = element_text(size = 18),
                       plot.title = element_text(size = 24),
                       legend.title = element_text(size = 22),
                       legend.text = element_text(size = 20))
  
  # For PxG plots needing big text  
  big_pxg_theme <<- theme(axis.title = element_text(size = 28), #26
                          axis.text = element_text(size = 26),
                          legend.title = element_text(size = 26),
                          legend.text = element_text(size = 24))
  
  # themes for Rmarkdown
  rmd_theme <<- theme(axis.title.x = element_text(size = 16), 
                      axis.text.x = element_text(size = 14), 
                      axis.title.y = element_text(size = 16),
                      axis.text.y = element_text(size = 14),
                      legend.title = element_text(size = 16),
                      legend.text = element_text(size = 14),
                      plot.title = element_text(size = 18))
  
  # theme for publication (multiple plots per figure)
  pub_theme <<- theme(legend.key.size = unit(0.03, 'npc'), # legend key is 0.03 of the plot size (basically)
                      axis.title = element_text(size = 8), 
                      axis.text = element_text(size = 6), 
                      legend.title = element_text(size = 8),
                      legend.text = element_text(size = 6),
                      legend.box.spacing = unit(0, 'pt')) # no space between plot and legend
  
  # theme for publication (multiple plots per figure)
  pub_theme_wt <<- theme_bw() + theme(legend.key.size = unit(0.03, 'npc'), # legend key is 0.03 of the plot size (basically)
                      axis.title = element_text(size = 8), 
                      axis.text = element_text(size = 6), 
                      legend.title = element_text(size = 8),
                      legend.text = element_text(size = 6),
                      legend.box.spacing = unit(0, 'pt')) # no space between plot and legend
  
  # theme for publication (one plot)
  pub_theme2 <<- theme(axis.title = element_text(size = 22), 
                       axis.text = element_text(size = 20), 
                       plot.title = element_text(size = 24))
  
  # theme for publication PxG plots 
  pub_pxg_theme <<- theme(axis.title = element_text(size = 22), 
                          axis.text = element_text(size = 18))
  # add axis.title.y = element_text(hjust = 0) or something to move y-axis away from plot
}

#--------------------------Calculate derived measures--------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# calc_auc: this function calculates an "area under the curve" metric, which is 
# really the area above the curve and below the horizontal line at the first data 
# point. 
# Input:
#     cross: r/qtl cross object 
#     col.name: column name to be used for the new aggregate phenotype 
#     steps: x-axis data 
#     phenos: y-axis data - phenotypes to be used to calculate the new aggregate 
#             phenotype
# Output:
#     cross: r/qtl cross object updated with new aggregate phenotype 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
calc_auc <- function(cross, steps, col.name, phenos){
  cross$pheno[,col.name] <- rep(NA,nrow(cross$pheno))
  
  #avgs <- colMeans(cross$pheno[,phenos], na.rm = T)
  
  for (i in 1:nrow(cross$pheno)){
    line <- cross$pheno[i,phenos]
    auc <- auc(x = steps, y = line)
    aac <- line[1]*tail(steps, n=1) - auc
    cross$pheno[i,col.name] <- aac 
  }
  
  return(cross)
}


#----------------------------Working with QTL objects--------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# create_models: this function creates all single-QTL models specified in the 
#                models table
# Input:
#     cross.obj: r/qtl cross object 
#     models: table containing information about phenotypes to create single-QTL 
#             models for 
#     covar: dataframe containing covariates 
#     method: whether to use existing R/qtl functions ("rqtl", default);
#             linear mixed modeling code ("lmm"); or
#             gene-by-treatment effect modeling code ("gxt")
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_models <- function(cross.obj, models, covar = NULL, method = "rqtl", 
                          trt = NULL, mod.dir){
  ensure_directory(mod.dir)
  
  for (i in 1:nrow(models)){
    filepath = paste(c(mod.dir, models[i,'obj'], ".Rdata"), collapse = "")
    
    if (file.exists(filepath)){ # model exists already, load .Rdata object 
      print(paste("Loading scanone object from", filepath))
      load(filepath, envir = .GlobalEnv)
    } else { # file doesn't exist, create model 
      print(paste("Creating scanone object for", models[i,'name'], "using", models[i,'type'], "model."))
      
      pheno.col <- models$colname[i]
      mod.type <- models$type[i]
      
      if ((method == "lmm") & (mod.type == "normal")){ 
        assign(models[i,'obj'], hk(cross.obj, pheno.col = pheno.col, envir = .GlobalEnv))
      } else if ((method == "gxt") & (mod.type == "normal")){  
        assign(models[i,'obj'], gxt(cross.obj, pheno.col = pheno.col, 
                                    covar = covar, trt = trt), envir = .GlobalEnv)
      } else { # use rqtl::scanone
        # define covariates, if specified
        if (!is.na(models$cov[i])){  
          addcovar <- covar[,models$cov[i]]
        } else { addcovar <- NULL }
        # create model 
        assign(models[i,'obj'], scanone(cross.obj, pheno.col = pheno.col, model = mod.type, 
                                        method = "hk", addcovar = addcovar), 
               envir = .GlobalEnv) 
      }
      
      print(paste("Saving object to", filepath))
      save(list = models[i,'obj'], file = filepath)
    }
  }
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# summary.permgev: summary() function for permutation object
#+++++++++++++++++++++++++++++++++++++++++++++++++++
summary.permgev <- function(object, ...){
  # Calculate 5% and 10% significance thresholds 
  thresholds <- matrix(qevd(p = c(0.95, 0.90), 
                            loc = object$results$par[1], 
                            scale = object$results$par[2], 
                            shape = object$results$par[3]))
  colnames(thresholds) <- c('lod')
  rownames(thresholds) <- c('5%', '10%')
  
  return(thresholds)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# create_perms: this function creates permutation objects for all single-QTL models 
#               specified in the models table 
# Input:
#     cross.obj: r/qtl cross object 
#     models: table containing information about phenotypes with single-QTL models 
#     covar: covariates 
#     perm.dir: directory with existing perm objects / to save perm objects to 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_perms <- function(cross.obj, models, covar = NULL, perm.dir, n.perm = 1000,
                         perm.strata = NULL, perm.Xsp = FALSE){
  
  ensure_directory(perm.dir)
  
  #--------------------Create or load permutation object-----------------------#
  for (i in 1:nrow(models)){
    perm <- models$perm.obj[i]
    perm.fname <- gsub("\\.", "", perm) # remove period from permutation object name, if present
    permfile <- paste(perm.dir, perm.fname, ".Rdata", sep = "")
    
    # create or load permutation object
    if(file.exists(permfile)){ # permutation object already exists, load from Rdata file 
      load(permfile, envir = .GlobalEnv)
      print(paste("Loaded", perm, "from file."))
    } else {
      print(paste("Creating", perm, "..."))
      
      addcovar <- NULL
      if(models$cov[i] != ""){
        addcovar <- covar[,models$cov[i]]
      }
      
      assign(perm, 
             scanone(cross.obj, pheno.col = models[i, 'colname'], 
                     model = models[i,'type'], addcovar = addcovar, n.perm = n.perm,
                     perm.strata = perm.strata, perm.Xsp = perm.Xsp), 
             envir = .GlobalEnv)
      save(list = perm, file = permfile)
    }
    
    #------------------------Fit GEV distribution(s)---------------------------#  
    # no X-chr-specific thresholds, 1 GEV distribution
    gevname <- paste0(perm, '.gev')
    print(paste("Creating GEV fit", gevname, "..."))
    fitgev <- fevd(as.data.frame(get(perm))$lod, type = "GEV")
    assign(gevname, fitgev, envir = .GlobalEnv)
    # add custom 'permgev' class so that we can use custom summary function to get GEV-defined thresholds  
    assign(gevname,
           add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
           envir = .GlobalEnv)
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# doc_peaks: this function finds significant LOD peaks, calculates 95% Bayes credible 
#           intervals for those peaks, and creates a summary table with data for 
#           each one, including marker, chromosome, position, LOD and positions of 
#           Bayes credible intervals.
# Input: 
#     models: table containing information about phenotypes with single-QTL models 
#     sig.level: significance level (can be 0.05 or 0.10)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
doc_peaks <- function(models, sig.level = 0.10){
  
  #--------------------Make sure all permutations exist------------------------#
  perms.exist = T
  for (m in 1:nrow(models)){
    if (!exists(models$perm.obj[m])){
      print(paste(models$perm.obj[m], "does not exist. All permutation objects must exist before attempting to identify significant peaks."))
      perms.exist = F
    }
  }
  if (perms.exist == F){stop("Missing permutation objects, execution halted.")}
  
  
  #-------------------------Count significant peaks----------------------------# 
  num.peaks = 0
  num.peaks.per.mod <- rep(NA, nrow(models))
  
  ### Delete?
  gevnames <- vector(mode = 'list', length = nrow(models))
  thresholds <- vector(mode = 'list', length = nrow(models))
  
  sig.ix <- ifelse(sig.level == 0.10, 2, ifelse(sig.level == 0.05, 1, NA)) 
  for (m in 1:nrow(models)){
    mod <- models$obj[m]
    perm <- models$perm.obj[m]
    
    # no X-chr specific thresholds, 1 GEV dist 
    num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
    
    gevnames[[m]] <- paste0(perm, '.gev') # GEV dist fit from perm object
    thresholds[[m]] <- summary(get(gevnames[[m]]))[sig.ix]
    num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
  }
  num.peaks <- sum(num.peaks.per.mod)
  
  
  #---------------------Get data for significant peaks-------------------------# 
  peak.data <- c('model', 'marker', 'chr', 'pos', 'lod', 'Bayes CI')
  peaks <- matrix(ncol = length(peak.data), nrow = num.peaks)
  colnames(peaks) <- peak.data
  
  row.ix = 0 # row in peaks table
  for (m in 1:nrow(models)){ # for each phenotype / model
    mod <- models$obj[m]
    perm <- models$perm.obj[m]
    gev <- gevnames[[m]]
    
    if (num.peaks.per.mod[m] == 0){
      print(paste("For model", models[m,'obj'], ", no LOD peaks above threshold"))
    } else {
      for (p in 1:num.peaks.per.mod[m]){ # for each peak 
        row.ix = row.ix + 1 # increment row index
        chrom = summary(get(mod), threshold = thresholds[[m]])[p,1] # get QTL chromosome
        # bayesint() produces a table with 3 rows (lower interval, peak, upper interval), columns = chr, pos, lod 
        bayesCI <- bayesint(get(models[m,'obj']), chr = chrom, prob = 0.95)
        # Populate peaks dataframe 
        peaks[row.ix,] <- c(mod, rownames(bayesCI)[2], chrom, bayesCI[2,2], bayesCI[2,3], 
                            paste(bayesCI[1,2],"-",bayesCI[3,2]))
      } 
    }
  }
  
  peaks <- as.data.frame(peaks)
  
  # Remove leading/trailing white space 
  for (i in 1:ncol(peaks)){
    peaks[,i] <- trimws(peaks[,i], which = c("both"))
  }
  
  return(peaks)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_unadj_pval: this function gets the unadjusted p-value from the linear regression
#           between a particular phenotype and the genotype at a particular marker.
# Input: 
#     cross: 
#     pheno.col: 
#     marker:
#     covariates 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_unadj_pval <- function(cross, pheno.col, marker, covar = NULL){
  pheno <- pull.pheno(cross, pheno.col) ### COMMENT THIS OUT??
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar), "", " + "),
                    paste(covar, collapse = " + "))
  
  if (is.null(covar)){
    lmdata <- data.frame(WB02$pheno[,pheno.col], marker.geno)
    names(lmdata)[1] <- pheno.col # have to rename column 
  } else {
    lmdata <- data.frame(WB02$pheno[,c(pheno.col,covar)], marker.geno)
  }
  
  lmod <- lm(formula, lmdata)
  unadj_pval <- summary(lmod)$coefficients['marker.geno', 'Pr(>|t|)']
  return(unadj_pval)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_unadj_pval: this function gets the unadjusted p-value from the linear regression
#           between a particular phenotype and the genotype at a particular marker.
# Input: 
#     cross: 
#     pheno.col: 
#     marker:
#     covariates 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_var_expl <- function(cross, pheno.col, marker, covar = NULL){
  pheno <- pull.pheno(cross, pheno.col)
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar), "", " + "),
                    paste(covar, collapse = " + "))
  
  if (is.null(covar)){
    lmdata <- data.frame(WB02$pheno[,pheno.col], marker.geno)
    names(lmdata)[1] <- pheno.col # have to rename column 
  } else {
    lmdata <- data.frame(WB02$pheno[,c(pheno.col,covar)], marker.geno)
  }
  
  lmod <- lm(formula, lmdata)
  sigma2 <- summary(lmod)$sigma^2
  
  uAA <- mean(pheno[which(marker.geno == 1)]) # 1 = AA
  uAB <- mean(pheno[which(marker.geno == 2)]) # 2 = AB 
  a <- uAB - uAA 
  
  h2 <- a^2 / (a^2 + 4*sigma2)
  
  return(h2)
}


#-----------------------------------Plotting-----------------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans: this function plots genome scans for all models specified in the 
#             models table. Adds a solid line for 5 % significance threshold and 
#             a dotted line for 10 % significance threshold, if the permutation
#             object exists. 
# Input:
#     models: table containing information about phenotypes with single-QTL models 
#     ylim: y limits, if they need to be consistent across plots. If not set, will 
#           use the default. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
plot_scans <- function(models, save = FALSE, save.dir = NULL, ...){
  for (m in 1:nrow(models)){
    mod_name <- models$obj[m] # model object name
    mod <- get(mod_name) # model object 
    perm <- models$perm.obj[m] # permutation object
    
    if (save == TRUE){ # saving, so open connection to png and make axis labels bigger 
      fname <- paste0(save.dir, mod_name, '.png')
      png(fname, width = 750)
      plot(mod, ylab = "", xlab = "", main = models$name[m], 
           alternate.chrid = T, bandcol = "gray90", cex.main = 2, cex.axis = 2, ...) 
      title(ylab = "LOD", line = 2.5, cex.lab = 2)
      title(xlab = "Chromosome", cex.lab = 2.3, line = 3) # cex.lab was 2
    } else { # not saving, just printing
      plot(mod, ylab = "LOD", main = models$name[m], alternate.chrid = T, 
           bandcol = "gray90", ...)
    }
    
    # threshold lines are the same regardless of whether we are saving or printing
    if (exists(perm)){
      permgev <- paste0(perm, '.gev')
      perm.sum <- summary(get(permgev))
      abline(h = c(perm.sum[1], perm.sum[2]), lty=1:2)
    }
    
    if (save == TRUE){ dev.off() } # close connection 
  }
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots for significant peaks
# Input:
#     cross.obj: r/qtl cross object 
#     cross.type: bc or f2
#     raw.data: raw phenotypes (cross.obj$pheno before transformations)
#     geno.map: list mapping A and B allele to more informative labels
#     peaks: dataframe of significant peaks 
#     qtl_map: dataframe mapping chromosome to QTL name (for plot titles)
#     theme: ggplot theme to use for PxG
#+++++++++++++++++++++++++++++++++++++++++++++++++++
plot_pxg <- function(cross.obj, cross.type, raw.data, peaks, qtl_map, 
                     theme, type, save = FALSE, save.dir = NULL, ...){
  
  for (k in 1:nrow(peaks)){
    
    if (save == TRUE){ # saving, so open connection to png 
      plotname <- paste0(peaks[k,'model'], '-chr', peaks[k,'chr'])
      fname <- paste0(save.dir, plotname, '.png')
      png(fname, width = 550)
    }
  
    pheno.col = models[which(models$obj == peaks[k,'model']), 'colname']
    pheno = raw.data[,pheno.col]
    marker.name = peaks[k,'marker']
        
    p <- pxg(cross.obj, pheno = pheno, 
               marker = marker.name,  
               ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
               theme = theme, ...)
    print(p)
    
    if (save == TRUE){ dev.off() }
  }
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots with ggplot. 
# Alternative to qtl::plotPXG, uses ggplot and more intuitive error bars  
# Input:
#     cross: r/qtl cross object 
#     pheno: vector of phenotype values for all mice 
#     marker: marker whose genotype to plot
#     geno.map: list mapping A and B allele to more informative labels 
#     qtl.map: list mapping chromosomes to QTL names for x-axis label
#     xlab: x axis label
#     ylab: y axis label 
#     title: title
#     theme: ggplot theme
#     type: scatter (default), boxplot or violin 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg <- function(cross, pheno, marker, geno.map, qtl.map = NULL, xlab = NULL, ylab, ylim = NULL,
                title = NULL, bestfit = TRUE, theme = rmd_theme, type = 'scatter', ...){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if ((chr == 'X')&('M' %in% levels(cross$pheno$sex))){bestfit = FALSE} 
  
  if(is.null(xlab) & !is.null(qtl.map)){
    qtl.name <- qtl.map$qtl_name[qtl.map$chr == chr]
    xlab <- as.expression(bquote(italic(.(qtl.name))*" (chr"*.(chr)*") "*Genotype))
  } else if (is.null(xlab) & is.null(qtl.map)){
    xlab <- "Genotype"
  }
  
  cross.imp <- fill.geno(cross) 
  marker.genos <- cross.imp$geno[[chr]]$data[,marker]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos,
                   pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotypes to genotype labels                 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      X.data <- data.frame(sex = sex, geno = df$geno, pgm = pgm)
      X.data$X.geno <- with(df, 
                            ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                            ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                            ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                            ifelse(sex==1 & geno==1, 'AY', 'BY'))))))
      
      df$geno <- mapvalues(X.data$X.geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
      geom_jitter(width = 0.1, height = 0) + 
      scale_y_continuous(limits = ylim) + 
      geom_errorbar(data = dfsum, mapping = aes(x = geno, y = mean, ymin = mean, ymax = mean), 
                    width = 0.2, color = brewer.pal(n=5,name='Dark2')[5], size = 2) + 
      geom_errorbar(data = dfsum, mapping = aes(x = geno, y = mean, ymin = mean-sd, ymax = mean+sd), 
                    width = 0.3, color = brewer.pal(n=5,name='Dark2')[5], size = 2) +
      labs(title = title, y = ylab, x = xlab) + 
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = brewer.pal(n=5,name='Dark2')[4])} + 
      theme
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_boxplot(fill = "#1B9E77", notch = TRUE) + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  } else if (type == 'violin'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_violin(fill = "#1B9E77") + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((type == 'scatter') & (bestfit == TRUE)){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}





#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for coloring dots by 2nd marker genotype
# Remove color-by-imputed-or-not because we're coloring by the genotype at a 2nd marker 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_colbym <- function(cross, pheno, marker, geno.map, xlab = 'Genotype', ylab, 
                title = NULL, bestfit = FALSE, theme = rmd_theme, marker2,
                type = 'scatter'){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # what chromosome is marker2 on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker2)
  chr2 <- names(cross$geno)[o]
  
  cross <- fill.geno(cross)
  marker.genos <- cross$geno[[chr]]$data[,marker]
  marker2.genos <- cross$geno[[chr2]]$data[,marker2]
  
  df <- data.frame(geno = marker.genos,
                   geno2 = marker2.genos,
                   pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotypes to genotype labels                 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
    df$geno2 <- mapvalues(df$geno2, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    
    ### NEED TO VALIDATE THIS
    
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      geno <- df$geno
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  df$geno2 <- factor(df$geno2, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, col = geno2)) + 
      geom_jitter(width = 0.1, height = 0) + 
      labs(title = title, y = ylab, x = xlab) + 
      scale_color_brewer(palette = 'Dark2', name = bquote(italic(Qpa2)~Genotype)) +
      theme
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, fill = geno2)) + 
      geom_boxplot() +
      scale_fill_brewer(palette = 'Dark2', name = bquote(italic(Qpa2)~Genotype)) + 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((bestfit == TRUE) & (type == 'scatter')){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}






