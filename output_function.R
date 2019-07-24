

library(knitr)
library(ggplot2)
library (rmarkdown)
library(dplyr)     
library(tidyr)
library(clinfun)
library(tidyselect)
library(MASS)
library(boot)
library(pROC)
library(rms)
library(ncvreg)
library(biostatR)

# ------------------------------
# univariate logistic regression
# output wald type confidence intervals
# ------------------------------
# a function to output OR(95%CI) and Wald p value for a single variable
# input the logistic regression model (glm)
sumfit.fun <- function(fit1){
  smfit1 <- summary(fit1)
  mp1 <- round(smfit1$coefficients[-1,4], 3)
  e1 <- round(exp(smfit1$coefficients[-1,1]), 
              2)
  lb1 <- round(exp(smfit1$coefficients[-1,1] - 
                     1.96*smfit1$coefficients[-1,2]), 
               2)
  up1 <- round(exp(smfit1$coefficients[-1,1]+ 
                     1.96*smfit1$coefficients[-1,2]), 
               2)
  mOR1 <- paste(e1," (",lb1,", ",up1,")", 
                sep="")
  cbind(mOR1,mp1)
}

# output uva table for a list of variables
# input the variable list and a single outcome
uvaout.logistic <- function(varlist,outcome){
  outab <- NULL
  for (i in 1:ncol(varlist)){
    fit1 <- glm(outcome ~ varlist[,i], 
                family="binomial")
    coef1 <- sumfit.fun(fit1)
    # for categorical var
    # level 1 is the referece level
    if (is.character(varlist[,i]) | is.factor(varlist[,i]))
      rownames(coef1) <- names(table(varlist[,i]))[-1] 
    else rownames(coef1) <- names(varlist)[i]
    outab <- rbind(outab, coef1)
  }
  # rownames(outab) <- names(varlist)
  colnames(outab) <- c("OR (95%CI)","p Value")
  return(outab)
}

# ------------------------------------
# function output mva logistic results
# ------------------------------------
mvaout.logistic = function(mv1){
  sum.mv1= summary(mv1)
  if (nrow(sum.mv1$coefficients)==2) {
    sum1 = c(sum.mv1$coefficients[-1,c(1,2)],round(sum.mv1$coefficients[-1,4],3))
    sumout = matrix(c(paste(round(sum1[1],2)," (",round(sum1[1]-1.96*sum1[2],2),", ",round(sum1[1]+1.96*sum1[2]),")",sep=""),sum1[4]),1,2)
  } else {
    sum1 = data.frame(sum.mv1$coefficients[-1,c(1,2)],round(sum.mv1$coefficients[-1,4],3))
    sumout = data.frame(paste(round(exp(sum1[,1]),2)," (",round(exp(sum1[,1]-1.96*sum1[,2]),2),", ",round(exp(sum1[,1]+1.96*sum1[,2]),2),")",sep=""),sum1[,3])
  }
  colnames(sumout) = c("OR (95%CI)","p Value")
  rownames(sumout) = rownames(sum.mv1$coefficients)[-1]
  return(data.frame(sumout))
}

# --------------------------------------------------
# univariate Cox regression for continuous variables
# --------------------------------------------------
uni.cox.table <- function(timefup, status, varlist, cindex=F){
  tab = NULL
  for (i in 1:ncol(varlist)){
    c1 = coxph(Surv(timefup, status)~varlist[,i])
    p1 = summary(c1)
    s1 = round(p1$conf.int[c(1,3,4)],2)
    pvalue = print.pv(p1$waldtest[3])
    
    if (cindex == T) {
      cin = coxphCPE(c1)
      cint <- paste(round(cin[1],3)," (",
                    round(cin[1]-1.96*cin[3],3),", ",
                    round(cin[1]+1.96*cin[3],3),")",
                    sep="")
      tab = rbind(tab, 
                  unname(c(colnames(varlist)[i],p1$n,
                           paste(s1[1],"(",s1[2],", ",s1[3],")",
                                 sep=""),
                           pvalue,cint)))
    } else
      tab = rbind(tab, 
                  unname(c(colnames(varlist)[i],p1$n,
                           paste(s1[1],"(",s1[2],", ",s1[3],")",sep=""),
                           pvalue)))
  }
  if (cindex==F) colnames(tab)=c("Variable","N","HR (95%CI)","p Value") else colnames(tab)=c("Variable","N","HR(95%CI)","p Value","c-index (95%CI)")
  return((tab))
}

# -----------------------------------------------------
# Univariate Cox regression for categorical variables
# Output HR and Wald test p value, 
# -----------------------------------------------------
uni.cox.cat.table <- 
  function(timefup, status, varlist, cindex = F ){
    tab = NULL
    for (i in 1:ncol(varlist)){
      tdat = data.frame(timefup, status, varlist[,i])
      tdat.nona = tdat[!is.na(varlist[,i]),]
      c1=coxph(Surv(timefup, status)~factor(varlist...i.),data=tdat.nona)
      p1 = summary(c1)
      s1 = unname(round(p1$conf.int[,c(1,3,4)],2))
      if (is.matrix(s1)) {s2 = s1} else {s2 = t(as.matrix(s1))}
      c0 = coxph(Surv(timefup, status)~1, data=tdat.nona)
      pvalue = print.pv(anova(c0,c1)$P[2])
      if (cindex == T) {
        cin <- coxphCPE(c1)
        cint <- paste(round(cin[1],3)," (",
                      round(cin[1]-1.96*cin[3],3),", ",
                      round(cin[1]+1.96*cin[3],3),")",
                      sep="")
        tab <- rbind(tab, 
                     c(colnames(varlist)[i],rep("",2),
                       pvalue,
                       cint),
                     cbind(names(table(varlist[,i])), table(varlist[,i]),
                           c("1",paste(s2[,1]," (",s2[,2],", ",s2[,3],")",sep="")),
                           rep("",length(table(varlist[,i]))),
                           rep("",length(table(varlist[,i])))
                     ))
        colnames(tab) <- c("Variable","N","HR (95%CI)","p Value", "c Index")
        
      } else {
        tab = rbind(tab, c(colnames(varlist)[i],rep("",2),pvalue),
                    cbind(names(table(varlist[,i])), table(varlist[,i]),
                          c("1",paste(s2[,1]," (",s2[,2],", ",s2[,3],")",sep="")),
                          rep("",length(table(varlist[,i])))))
        colnames(tab) <- c("Variable","N","HR (95%CI)","p Value")
      }
    }
    rownames(tab) <- NULL
    return(tab)
  }


# -------------------------------
# function output mva cox results
# -------------------------------
mvaout.cox = function(mv1){
  sum.mv1= summary(mv1)
  if (nrow(sum.mv1$coefficients)==1) {
    sum1 = c(round(sum.mv1$conf.int[,c(1,3,4)],2),round(sum.mv1$coefficients[,5],3))
    sumout = matrix(c(paste(sum1[1]," (",sum1[2],", ",sum1[3],")",sep=""),sum1[4]),1,2)
  } else {
    sum1 = data.frame(round(sum.mv1$conf.int[,c(1,3,4)],2),round(sum.mv1$coefficients[,5],3))
    sumout = data.frame(paste(sum1[,1]," (",sum1[,2],", ",sum1[,3],")",sep=""),sum1[,4])
  }
  colnames(sumout) = c("HR (95%CI)","p Value")
  rownames(sumout) = rownames(sum.mv1$coefficients)
  return(data.frame(sumout))
}

# ----------------------------------------------
# revise fmt_pvalue function in biostatR package
# ----------------------------------------------
# edit fmt_pvalue function in biostatR
fmt_pvalue <-
  function (x, digits = 1, prepend_p = FALSE) {
    if (digits == 2) {
      p_fmt <- dplyr::case_when(
        x > 1 ~ NA_character_, x < 
          0 ~ NA_character_, x > 0.99 ~ ">0.99", round(x, 2) >= 
          0.1 ~ sprintf("%.2f", x), x >= 0.001 ~ sprintf("%.3f", x), x < 0.001 ~ "<0.001")
    }
    if (digits == 1) {
      p_fmt <- dplyr::case_when(
        x > 1 ~ NA_character_, x < 
          0 ~ NA_character_, x > 0.9 ~ ">0.9", round(x, 1) >= 
          0.2 ~ sprintf("%.1f", x), round(x, 2) >= 0.1 ~ sprintf("%.1f", x), 
        x >= 0.001 ~ sprintf("%.3f", x), x < 0.001 ~ "<0.001")
    }
    if (digits == 3) {
      p_fmt <- dplyr::case_when(x > 1 ~ NA_character_, x < 
                                  0 ~ NA_character_, x > 0.999 ~ ">0.999", round(x, 1) >= 
                                  0.2 ~ sprintf("%.3f", x), round(x, 2) >= 0.1 ~ sprintf("%.3f", 
                                                                                         x), x >= 0.001 ~ sprintf("%.3f", x), x < 0.001 ~ 
                                  "<0.001")
    }
    # if (prepend_p == TRUE) {
    #     p_fmt <- dplyr::case_when(is.na(x) ~ NA_character_, 
    #         stringr::str_sub(x, end = 1L) %in% c("<", ">") ~ 
    #             paste0("p", x), TRUE ~ paste0("p=", x))
    # }
    return(p_fmt)
  }

# ----------------
# p value output
# 3 decimal points
# ----------------
# a function round and present a single p value with 3 decimal points in a table/txt
print.pv <- function(pvalue) {
  if (!is.numeric(pvalue)) stop("Input value is not numeric.")
  if (is.vector(pvalue)){
    if (is.na(pvalue)) pvalue3=pvalue
    else if (pvalue < 0.001) pvalue3 <- "<0.001"
    else if (pvalue > 0.999) pvalue3 <- 0.999
    else if (is.na(pvalue)) pvalue3 <- NA
    else pvalue3 <- sprintf("%.3f",round(pvalue,3))    
    return(pvalue3)
  } 
}
