# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Author: Jeremias N. Brand
# This is a collection of helper functions I use for comparative analysis.
# Use them at your own risk!


pgls_diag_gls <- function(name, mod_type, model, phy, formula, data,correlation, weight, method, control) {
  # http://blog.phytools.org/2013/02/a-comment-on-distribution-of-residuals.html
  m <- model
  par <- m$modelStruct[[1]]
  if(mod_type == "lambda"){
    pRes <- chol(solve(vcv(corPagel(par,phy=phy,fixed=TRUE)))) %*% residuals(m)
  } else if (mod_type == "OU"){
    pRes <- chol(solve(vcv(corMartins(par,phy,fixed=TRUE)))) %*% residuals(m)
  } else {
    pRes <- chol(solve(vcv(corBrownian(phy=phy)))) %*% residuals(m)
  }
  
  par(mfrow=c(2,2), oma=c(0.5,0.5,2,0.5), cex=0.8)
  # residual density plot
  plot(density(pRes), main=paste("Phylogenetic Residuals", mod_type, round(par,2)))
  # residual distribution
  plot(pRes ~ fitted(m), main="Phylogenetic Residuals ~ fitted values")
  # qq plot
  qqnorm(pRes)
  qqline(pRes)
  # parameter profile plot
  if(mod_type == "lambda"){
    optimize_lambda(formula=formula, data = data, phy=phy,
                    weight = weight, method = method, control=control)
    abline(v = par, col = "red")
  } else if (mod_type == "OU"){
    optimize_alpha(formula=formula, data = data, phy=phy,
                    weight = weight, method = method, control=control)
    abline(v = par, col = "red")
  } else {
    print("BM model")
  }
  
  mtext(paste0("PGLS Diagnostics of ",name), outer=TRUE,  cex=1, line=-0.5)
  }


optimize_lambda <- function(formula, data,phy, correlation,weight,method, form, control) {
  lam <- seq(0, 1, length.out = 500)
  lik <- sapply(lam,
                function(lambda) {
                  logLik(gls(formula, data = data, 
                             correlation = corPagel(lambda, phy=phy,form = form, fixed=TRUE), 
                             weight = weight, method = method,
                             control=control))
                })
  plot(lik ~ lam, type = "l", main = expression(paste("Likelihood Plot for ",lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
  return(lam[which.max(lik)])
  }

optimize_alpha <- function(formula, data,phy, correlation,weight,method,form, control) {
  lam <- seq(0, 10, length.out = 500)
  lik <- sapply(lam,
                function(lambda) {
                  logLik(gls(formula, data = data, 
                             correlation = corMartins(lambda, phy=phy,form=form,fixed=TRUE), 
                             weight = weight, method = method,
                             control=control))
                })
  plot(lik ~ lam, type = "l", main = expression(paste("Likelihood Plot for ",alpha)), ylab = "Log Likelihood", xlab = expression(alpha))
  return(lam[which.max(lik)])
}

AICc_pgls_nlme <- function(mod){
  k <- attr(logLik(mod),"df")
  lnL <- as.numeric(mod$logLik)
  n <- mod$dims$N
  aic <- round(2 * k - 2 * lnL,2)
  aicc <- round(2 * k * (n/(n - k - 1)) - (2 * lnL),2)
  return(c(AIC=aic, AICc=aicc))
}


AIC_ace <- function (v, n) 
{
  # AIC function from geiger adapted for an ace object
  k <- length(v[["rates"]])
  lnL <- v[["loglik"]]
  v$aic <- 2 * k - 2 * lnL
  v$aicc <- 2 * k * (n/(n - k - 1)) - (2 * lnL)
  v$k <- k
  v$n <- n
  return(v)
}

plot_MK <- function(file, trait, tree, cols) {
  library(phytools)
  pdf(file = file, width = 3*6, height = 20)
  margin = c(2.5, 2, 2, 2)
  piesize = 0.8
  title_size = 5
  label_size = 1.5
  par(mfrow=c(1,3), oma=c(0,0,1,0))
  
  # we remove tips that are NA
  to_drop <- names(trait)[is.na(trait)]
  
  tree <- drop.tip(tree, to_drop)
  trait <- trait[!is.na(trait)]

  fitER <- ace(trait,tree,model="ER",type="discrete")
  fitER
  # then these are the estimates for each node
  round(fitER$lik.anc,3)
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=fitER$lik.anc, piecol=cols, cex=piesize)
  tiplabels(pie=to.matrix(trait, sort(unique(trait))), piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main="MK - equal rates", outer = F, cex = title_size)
  
  
  # Here we need to used all numbers
  fitSYM <- ace(trait,tree,model="SYM",type="discrete")
  fitSYM
  # then these are the estimates for each node
  round(fitSYM$lik.anc,3)
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=fitSYM$lik.anc, piecol=cols, cex=piesize)
  tiplabels(pie=to.matrix(trait, sort(unique(trait))), piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main="MK - symetrical rates", outer = F, cex = title_size)
  
  
  fitARD <- ace(trait,tree,model="ARD",type="discrete")
  fitARD
  round(fitARD$lik.anc,3)
  plotTree(tree,type="phylogram", fsize = label_size,ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=fitARD$lik.anc,piecol=cols, cex=piesize)
  tiplabels(pie=to.matrix(trait, sort(unique(trait))), piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  # quite similar results
  title(main="MK - ARD", outer = F, cex = title_size)
  dev.off()
  
  # add AIC and AICc to the models
  n <- length(trait)
  fitER <- AIC_ace(fitER, n)
  fitSYM <- AIC_ace(fitSYM, n)
  fitARD <- AIC_ace(fitARD, n)
  
  return(list(fitER = fitER, fitSYM = fitSYM, fitARD = fitARD, trait = trait, tree = tree))
}


plot_MK_onemod <- function(file, trait, tree, model, cols) {
  library(phytools)
  pdf(file = file, width = 1*6, height = 20)
  margin = c(2.5, 2, 2, 2)
  piesize = 0.8
  title_size = 5
  label_size = 1.5
  par(mfrow=c(1,1), oma=c(0,0,1,0))
  
  # we remove tips that are NA
  to_drop <- names(trait)[is.na(trait)]
  tree <- drop.tip(tree, to_drop)
  trait <- trait[!is.na(trait)]

  fitmod <- ace(trait,tree,model=model,type="discrete")
  # then these are the estimates for each node
  round(fitmod$lik.anc,3)
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=fitmod$lik.anc, piecol=cols, cex=piesize)
  tiplabels(pie=to.matrix(trait, sort(unique(trait))), piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main=paste("Model:", deparse(substitute(model))), outer = F, cex = title_size)
  dev.off()
  
  n <- length(trait)
  fitmod <- AIC_ace(fitmod, n)
  
  return(list(mod = fitmod, trait = trait, tree = tree))
}



plot_simmap <- function(file, trait, tree, nsim=100, cols) {
  prob_matrix <- to.matrix(trait, sort(unique(trait)))
  No_states <- length(na.omit(unique(trait)))
  prob_matrix[names(trait[is.na(trait)]),] <- rep(1/No_states, No_states)
  
  
  ER_trees<-make.simmap(tree, prob_matrix, nsim=nsim, model="ER", Q = "mcmc")
  SYM_trees<-make.simmap(tree, prob_matrix, nsim=nsim, model="SYM", Q = "mcmc")
  ARD_trees<-make.simmap(tree, prob_matrix, nsim=nsim, model="ARD", Q = "mcmc")
  
  
  ER<-describe.simmap(ER_trees,plot=FALSE)
  SYM<-describe.simmap(SYM_trees,plot=FALSE)
  ARD<-describe.simmap(ARD_trees,plot=FALSE)
  
  
  pdf(file = file, width = 3*6, height = 20)
  par(mfrow=c(1,3), oma=c(0,0,1,0))
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=ER$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main="stochastic character mapping - ER", outer = F, cex = title_size)
  
  
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=SYM$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main="stochastic character mapping - SYM", outer = F, cex = title_size)
  
  
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=ARD$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main="stochastic character mapping - ARD", outer = F, cex = title_size)
  dev.off()
  
  return(list(ER_trees = ER_trees, SYM_trees = SYM_trees, ARD_trees = ARD_trees,
              ER = ER, SYM = SYM, ARD = ARD))
}

project.phylomorphospace_mod <- function(tree, X, nsteps = 200, sleep = 0, direction = c("to", 
                                                          "from", "both"), ...) 
{
  lab = "off"
  direction <- direction[1]
  #tree <- minRotate(reorder(tree, "cladewise"), X[, 2], print = FALSE)
  X <- X[tree$tip.label, ]
  A <- cbind(fastAnc(tree, X[, 1]), fastAnc(tree, X[, 2]))
  cladogram <- tree
  cladogram$edge.length <- NULL
  mar <- par()$mar
  plotTree(cladogram, type = "cladogram", nodes = "centered", 
           plot = FALSE)
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  obj$xx <- obj$xx * (max(c(X[, 1], A[, 1])) - min(c(X[, 1], 
                                                     A[, 1]))) + min(c(X[, 1], A[, 1]))
  xlim <- range(obj$xx)
  xlim[2] <- xlim[2] + abs(diff(xlim)) * (obj$x.lim[2] - 1)
  obj$yy <- (obj$yy - min(obj$yy))/(max(obj$yy) - min(obj$yy)) * 
    (max(c(X[, 2], A[, 2])) - min(c(X[, 2], A[, 2]))) + min(c(X[, 
                                                                2], A[, 2]))
  ylim <- range(obj$yy)
  X0 <- cbind(obj$xx[1:Ntip(tree)], obj$yy[1:Ntip(tree)])
  rownames(X0) <- tree$tip.label
  A0 <- cbind(obj$xx[1:tree$Nnode + Ntip(tree)], obj$yy[1:tree$Nnode + 
                                                          Ntip(tree)])
  rownames(A0) <- 1:tree$Nnode + Ntip(tree)
  par(mar = mar, new = TRUE)
  dev.hold()
  if (direction %in% c("to", "both")){ 
    phylomorphospace(tree, X0, A0, label = lab, 
                     xlim = xlim, ylim = ylim, ...)
    ani.record()
    }
  else if (direction == "from") {
    phylomorphospace(tree, X, A, label = lab, xlim = xlim, 
                     ylim = ylim, ...)
  ani.record()
  dev.flush()}
  if (direction == "both") 
    nsteps <- ceiling(nsteps/2)
  if (direction %in% c("to", "both")) {
    for (i in 2:nsteps) {
      Sys.sleep(sleep)
      dev.hold()
      phylomorphospace(tree, ((nsteps - i) * X0 + i * X)/nsteps, 
                       ((nsteps - i) * A0 + i * A)/nsteps, xlim = xlim, 
                       ylim = ylim, label = lab,...)
      ani.record()
      dev.flush()
      
    }
  }
  if (direction %in% c("from", "both")) {
    for (i in (nsteps - 1):1) {
      Sys.sleep(sleep)
      dev.hold()
      phylomorphospace(tree, ((nsteps - i) * X0 + i * X)/nsteps, 
                       ((nsteps - i) * A0 + i * A)/nsteps, xlim = xlim, 
                       ylim = ylim, label = lab, ...)
      ani.record()
      dev.flush()
    }
    dev.hold()
    phylomorphospace(tree, X0, A0, label = lab, 
                     xlim = xlim, ylim = ylim, label = lab, ...)
    ani.record()
    dev.flush()
  }
  invisible(NULL)
}


plot_simmap_onemod <- function(file, trait, tree, model,
                                   nsim=100, burnin=100, samplefreq=100,
                                   width=6, height=25, cols, label_offset = 0.1, alpha=1, beta = 1,
                                   label_size=1, title_size=1, margin= c(0.5,0.5,2,0.5), piesize=2, trans_tbl=NULL) {
  library(ape)
  Ntips_scale <- length(trait)/100
  prob_matrix <- to.matrix(trait, sort(unique(trait)))
  No_states <- length(na.omit(unique(trait)))
  prob_matrix[names(trait[is.na(trait)]),] <- rep(1/No_states, No_states)
  
  trees<-make.simmap(tree, prob_matrix, nsim=nsim, model= model, Q = "mcmc",
                     burnin=burnin, samplefreq = samplefreq,
                     prior = list(alpha=alpha, beta=beta, use.empirical=FALSE))
  mod<-describe.simmap(trees,plot=FALSE)
  
  if (!is.null(trans_tbl)) {
    tree$tip.label <-  trans_tbl$sci_name_short[match(tree$tip.label, trans_tbl$Short)]
    tree$tip.label <- gsub(" ", "_", tree$tip.label)
  }
  
  pdf(file = file, width = width, height = height*Ntips_scale)
  par(mfrow=c(1,1), mar=margin, oma=c(0,0,0,0))
  plot(tree, type="phylogram", show.tip.label = TRUE, label.offset = label_offset)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=mod$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols, cex=piesize)
    add.simmap.legend(leg = sort(unique(trait)),
                    colors=cols, prompt=FALSE, x=0.2, y=4, fsize=0.8,
                    vertical = TRUE, cex=10)
  title(main=paste(deparse(substitute(model)),"nsim:", nsim, "N:", length(trait)), outer = F, cex = title_size)
  dev.off()
  
  return(list(trees = trees, mod = mod))
}

REplot_simmap <- function(file, var, tree, simmap, df, treename, type = "phylogram",
                          width=6, height=17, cols, label_offset = 0.04,
                          label_size=1, title_size=1, margin= c(0,0,0,0), piesize=0.5, trans_tbl=NULL) {
  
  drop <- tree$tip.label[!tree$tip.label %in% df[[treename]]]
  t <-  drop.tip(tree, drop)
  Ntips_scale <- length(t$tip.label)/100
  # get trait
  trait <- get_orderd_trait(t=t, df=df, var=var, treename=treename)
  prob_matrix <- to.matrix(trait, sort(unique(trait)))
  # find the nodes in the new tree
  node_loc <- matchNodes(simmap$trees[[1]],t,"distances")
  # Sci names (needs to be after node matching because that uses tips!)
  if (!is.null(trans_tbl)) {
    t$tip.label <-  trans_tbl$sci_name_short[match(t$tip.label, trans_tbl[[treename]])]
    t$tip.label <- gsub(" ", "_", t$tip.label)
  }
  
  pdf(file = file, width = width, height = height*Ntips_scale, useDingbats = FALSE)
  par(mfrow=c(1,1), mar=margin, oma=c(0,0,0,0))
  plot(t, type=type, show.tip.label = TRUE, label.offset = label_offset)
  nodelabels(node=node_loc[,2],
             pie=simmap$mod$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols, cex=piesize)
  dev.off()
}


REplot_simmap_fan <- function(file, var, tree, simmap, df, treename, type = "phylogram",
                              width=7, height=7, cols, cols2, label_offset = 0.04, f_size,
                              label_size=1, title_size=1, margin= c(0,0,0,0), piesize=0.5, trans_tbl=NULL) {
  
  
  drop <- tree$tip.label[!tree$tip.label %in% df[[treename]]]
  t <-  drop.tip(tree, drop)
  Ntips_scale <- length(t$tip.label)/100
  # get trait
  trait <- get_orderd_trait(t=t, df=df, var=var, treename=treename)
  prob_matrix <- to.matrix(trait, sort(unique(trait)))
  # find the nodes in the new tree
  node_loc <- matchNodes(simmap$trees[[1]],t,"distances")
  # Sci names (needs to be after node matching because that uses tips!)
  if (!is.null(trans_tbl)) {
    t$tip.label <-  trans_tbl$sci_name_short[match(t$tip.label, trans_tbl[[treename]])]
    t$tip.label <- gsub(" ", "_", t$tip.label)
  }
  
  pdf(file = file, width = width, height = height*Ntips_scale, useDingbats = FALSE)
  par(mfrow=c(1,1), mar=margin, oma=c(0,0,0,0))
  plot(t, type=type, show.tip.label = TRUE, cex=f_size,  label.offset = label_offset,
       edge.width = 1.5,
       rotate.tree = 100,
       open.angle = 2.5)
  nodelabels(node=node_loc[,2],
             pie=simmap$mod$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols2, cex=piesize)
  dev.off()
}



extract_Q <- function(mod) {
  Q <- lapply(mod, function(x) x$Q)
  return(Reduce("+", Q)/length(Q))
}

plot_simmap_onemod_bku <- function(file, trait, tree, model,
                               nsim=100, burnin=100, samplefreq=100,
                               width=6, height=20, cols,
                               label_size=1, title_size=1, margin= c(1,1,1,1), piesize=2) {
  prob_matrix <- to.matrix(trait, sort(unique(trait)))
  No_states <- length(na.omit(unique(trait)))
  prob_matrix[names(trait[is.na(trait)]),] <- rep(1/No_states, No_states)
  
  trees<-make.simmap(tree, prob_matrix, nsim=nsim, model= model, Q = "mcmc", burnin=burnin, samplefreq = samplefreq)
  mod<-describe.simmap(trees,plot=FALSE)

  print(prob_matrix)
  pdf(file = file, width = width, height = height)
  par(mfrow=c(1,1), oma=c(0,0,1,0))
  plotTree(tree,type="phylogram", fsize = label_size,
           ftype="i", mar = margin, offset=0.5)
  nodelabels(node=1:tree$Nnode+Ntip(tree),
             pie=mod$ace, piecol=cols, cex=piesize)
  tiplabels(pie=prob_matrix, piecol=cols, cex=piesize)
  add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                    y=-max(nodeHeights(tree)),fsize=0.8, vertical = FALSE)
  title(main="stochastic character mapping", outer = F, cex = title_size)
  dev.off()
  
  return(list(trees = trees, mod = mod))
}

read_bt_discrete_mcmc <- function(path) {
  require(coda)
  require(dplyr)
  BT_file <- path
  FullOut <-  scan(file = BT_file, what="c", quiet=T, sep="\n")
  BT_log <- read.table(BT_file, skip = (grep("Iteration	Lh", FullOut) - 1), sep = "\t", header = TRUE, quote="\"")
  BT_log <- select(BT_log, -X, -Tree.No, -Iteration)
  BT_log_mcmc <- mcmc(BT_log)
  return(BT_log_mcmc)
}


read_bt_discrete_marginalLH <- function(path) {
  require(coda)
  require(dplyr)
  BT_file <- path
  FullOut <-  scan(file = BT_file, what="c", quiet=T, sep="\n")
  BT_log <- read.table(BT_file, skip = (grep("Log marginal likelihood:", FullOut) - 1), sep = "\t", header = FALSE, quote="\"")
  return(BT_log$V2)
}

check_for_mag  <- function(x, ref) {
  item <- x[2]
  mag_needed <- x[1]
  !is.na(ref[,x[2]]) & is.na(ref[,x[1]])
}

mean_NA_real <- function(x, na.rm = TRUE){
  m <- mean(x, na.rm = na.rm)
  res <- ifelse(is.na(m), NA_real_, m)
  return(res)
}

sum_NA_real <- function(x, na.rm = TRUE){
  if (sum(is.na(x)) == 0) {
    return(NA_real_)
  }
  return(sum(x, na.rm = na.rm))
}

ICC_with_subsetting <- function(x,df) {
  # we want to calculate Inter class correlation
  # but only if there are enough intra class samples
  s <- na.omit(unique(df$Short))
  d <- na.omit(df[,c("Short","ID", x)])
  no_icc_samples <- sum(table(factor(d$ID)) >= 2)
  if (no_icc_samples > 1){
    icc_res <- ICCest(factor(d[,"ID"]), d[,x])
    # add specific info
  } else {
    icc_res <- rep(list(NA), 7)
  }
  icc_res$N_ICC <- no_icc_samples
  #icc_res$Short <- s
  names(icc_res) <- c("ICC","LowerCI","UpperCI","N","k",
                      "varw","vara","N_ICC")
  return(unlist(icc_res))
}

get.landscape.FPK_df <- function (fit, Npts = 100, xlim = NULL, ylim = NULL) {
  # modification of the original function to just give points to plot independently
  if ("a" %in% names(fit$par)) {
    a = fit$par$a
  }
  else {
    a = fit$par_fixed$a
  }
  if ("b" %in% names(fit$par)) {
    b = fit$par$b
  }
  else {
    b = fit$par_fixed$b
  }
  if ("c" %in% names(fit$par)) {
    c = fit$par$c
  }
  else {
    c = fit$par_fixed$c
  }
  bounds = fit$par_fixed$bounds
  SEQ = seq(from = -1.5, to = 1.5, length.out = Npts)
  V = a * SEQ^4 + b * SEQ^2 + c * SEQ
  step = (bounds[2] - bounds[1])/(Npts - 1)
  if (is.null(xlim)) {
    xlim = c(bounds[1], bounds[2])
  }
  if (is.null(ylim)) {
    ylim = c(0, max(exp(-V)/sum(exp(-V) * step)))
  }
  return(data.frame(x = seq(from = bounds[1], 
                            to = bounds[2], length.out = Npts), y = (exp(-V)/sum(exp(-V) * step)))
  )
}


plot_tips <- function(df, tree, var, treenamecol="Short", palette) {
  trait <- as.character(df[,var])
  names(trait) <- df$sci_name_short
  trait <- ifelse(is.na(trait), "no observation", trait)
  plotTree(tree, type = "phylogram", fsize = 2, ftype = "i", offset = 0.4, mar = c(1,1,1,1))
  cols<-setNames(palette[1:length(unique(trait))], sort(unique(trait)))
  tiplabels(pie = to.matrix(trait, sort(unique(trait))), piecol = cols, cex = 0.6, vjust = 0.1)
  add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1],
                    y = 10 * max(nodeHeights(tree)), fsize = 2)
  title(main = var, outer = F, cex = 4, fsize = 20)
}


stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}


get_orderd_trait  <- function(t,df,var,treename) {
  order <- match(t$tip.label, df[,treename])
  trait <- as.character(df[order,var])
  names(trait) <- df[order,treename]
  trait <- ifelse(is.na(trait), "no observation", trait)
  return(trait)
}

get_pies_trait <- function(trait, palette){
  cols<-setNames(palette[1:length(unique(trait))], sort(unique(trait)))
  return(list(to.matrix(trait, sort(unique(trait))), cols))
}


prepare_phyl_pca <- function(df,tree, ignore=c(""), treename="Short") {
  library(dplyr)
  # produces a comparative dataset and gives reasons for omissions
  print("Producing comparative dataset an corresponding pruned tree.")
  print(tree)
  if(ignore != ""){
    print("Ignoring following columns for pruning:")
    print(ignore)
  }
  ppca <- df[df[,treename] %in% tree$tip.label, setdiff(names(df), ignore)]
  missing <- apply(is.na(ppca[setdiff(names(ppca), treename)]), 2, which)
  # remove empty lists from the results
  missing <- missing[lapply(missing, length)>0]
  print("The following columns have missing data:")
  taxa <- lapply(missing, function(x){as.character(ppca[x,treename])})
  print(taxa, row.names = FALSE)
  sum_taxa <- unique(do.call(c, taxa))
  print(paste(length(sum_taxa),"taxa will be removed:"))
  print(sum_taxa, row.names = FALSE)
  # remove missing data
  ppca <- na.omit(ppca)
  rownames(ppca) <- ppca[,treename]
  
  to_drop <- tree$tip.label[! tree$tip.label %in% rownames(ppca)]
  t <- drop.tip(tree, tip = to_drop)
  print(paste("pruned data has:",length(t$tip.label), "tips"))
  # add the ignored columns back in 
  if(ignore != ""){
    ppca <- left_join(ppca, df[,c(treename, ignore)])
    ppca <- select(ppca, one_of(c(treename, ignore)),
                   everything()) %>%
      as.data.frame()
    rownames(ppca) <- ppca[,treename]
  }
  return(list(ppca, t))
}


find_influ <- function(df) {
  SDval = sd(df$val)
  res <- df %>%
    mutate(raw_D = val - full_model,
           stand_D = raw_D/SDval,
           pcnt_D = round(abs(abs(raw_D)/full_model) * 100, 1))
  return(res)
}


pgls_for_sample_size_se <- function(formula, tree, df, SScol, min_species = 10, method = "ML", model = "BM", treename = "Short", se=NULL, startval=NULL) {
  require(dplyr)
  require(phytools)
  require(rr2)
  correlation = ifelse(model == "BM","corBrownian",
                       ifelse(model == "OU","corMartins",
                              ifelse(model == "lambda","corPagel","err")
                       ))
  
  correlationgls = ifelse(model == "BM","corBrownian(1, phy = tree)",
                          ifelse(model == "OU","corMartins(1, phy = tree, fixed = FALSE)",
                                 ifelse(model == "lambda","corPagel(1, phy = tree, fixed = FALSE)","err")
                          ))
  # Sometimes pgls.SEy fails so then we can use gls this is the switch variable for this
  USEGLS = FALSE
  
  if (correlation == "err") {stop("Don't recognize evolutionary model. must be BM,OU or lambda")}
  if (class(formula) != "formula") 
    stop("formula must be class 'formula'")
  
  # incooperate SE in y
  if (!is.null(se)){
    fullse = setNames(df[, se], df[,treename])
  } else {
    fullse = NULL
  }
  pgls <- tryCatch(
    {if(!is.null(startval)){
      pgls.SEy_startval(formula, data = df, tree = tree,
               corClass =eval(parse(text=correlation)), se=fullse,
               method = method, startval = startval, fixed = FALSE)
    } else {
      pgls.SEy(formula, data = df, tree = tree,
               corClass =eval(parse(text=correlation)), se=fullse, method = method)
    }
   },
    error=function(cond){
      message("\npgls.SEy throws an error we will use gls() in this analysis!")
      message(paste("ERROR:", cond, "\n"))
      # switch to gls
      return(
        tryCatch(
          {gls(model = formula,
               correlation = eval(parse(text=correlationgls)),
               data = df, method = method)},
          error=function(cond) {
            message("gls also did not work!")
            stop("aborting run!")
          }
        )
      )},
    warning=function(cond){
    })
  
  if(is.null(pgls$call$weights)) {USEGLS = TRUE} 
  
  print(summary(pgls))
  r2pred = round(rr2::R2.pred(pgls),3)
  
  Ns = sort(unique(df[,SScol]))
  
  s_base <- summary(pgls)
  print(s_base)
  nVars = nrow(s_base$tTable)
  
  # pgls_table
  res_df <- as.data.frame(matrix(NA, ncol = 2+(4*nVars), nrow = length(Ns)))
  s <- as.data.frame(summary(pgls)$tTable)
  s$par <- row.names(s); s$Min_sample_size <- min(Ns); s$N_species_in_model <- length(tree$tip.label)
  s_res <- reshape(s, idvar = c("Min_sample_size", "N_species_in_model"), timevar = "par", direction = "wide")
  nam <- gsub("\\(|\\)","",names(s_res))
  nam <- gsub("Std.Error","SE",nam)
  nam <- gsub("(.+)\\.(.+)","\\2:\\1",nam)
  names(res_df) <- nam
  
  original_values <- as.data.frame(t(s_res[2:length(s_res)]))
  original_values$par <- as.character(nam[2:length(nam)])
  names(original_values) <- c("full_model", "stat")
  original_values$stat <- as.factor(original_values$stat)
  prg <- txtProgressBar(0,length(Ns), initial = 0, style = 2)
  for (i in 1:length(Ns)){
    setTxtProgressBar(prg, i)
    # remove everything with lower sample size
    index = df[,SScol] >= Ns[[i]]
    r_df <- df[index,]
    if(nrow(r_df) < min_species){
      res <- c(Ns[[i]], nrow(r_df), rep(NA,(4*nVars)))
    } else {
      rownames(r_df) <- r_df[,treename]
      if(sum(index != FALSE) == 0) {
        r_tree <- drop.tip(tree, tip = df[!index,treename])#1
      }else{
        r_tree <- tree 
        }
      
      if (!is.null(se)){
        r_se = setNames(r_df[, se], r_df[,treename])
      } else {
        r_se = NULL
      }
      if(USEGLS == FALSE) {
        pgls_r <- tryCatch(
          {pgls.SEy(formula, data = r_df, tree = r_tree,
                    corClass =eval(parse(text=correlation)), se=r_se, method = method)},
          error=function(cond){
            message("\npgls.SEy did not converge, trying gls")
            print(r_df)
            message(paste("ERROR:", cond))
            return(
              tryCatch(
                {gls(model = formula,
                     correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
                     data = r_df, method = method)},
                error=function(cond) {
                  message("\ngls also did not work! Skipping this iteration.")
                  return(NA)
                }
              ))},
          warning=function(cond){
          })
      } else {
        pgls_r <- tryCatch(
          {gls(model = formula,
               correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
               data = r_df, method = method)},
          error=function(cond) {
            message("\ngls did not work! Skipping this iteration.")
            return(NA)
          })
        }

      if(is.na(pgls_r[1])){
        res <- c(Ns[[i]], nrow(r_df), rep("MODEL FAIL",(4*nVars)))
      } else {
        s <- as.data.frame(summary(pgls_r)$tTable)
        s$par <- row.names(s)
        s$Min_sample_size <- Ns[[i]]
        s$N_species_in_model <- nrow(r_df)
        res <- reshape(s, idvar = c("Min_sample_size", "N_species_in_model"), timevar = "par", direction = "wide")
        names(res) <- nam
      }
    }
    res_df[i,] <- res
  }
  close(prg)
  suppressWarnings(
  res_sum <- res_df[,2:length(res_df)] %>%
    select(-N_species_in_model) %>%
    gather(key = stat, value = val) %>%
    mutate(val = as.numeric(val),
           stat = factor(stat, levels = nam)) %>%
    group_by(stat)%>%
    summarise_all( funs(mean = mean(., na.rm =TRUE), 
                        sd = sd(., na.rm =TRUE),
                        n = n())) %>%
    mutate(se = sd / sqrt(n),
           CI_low = mean - qnorm(0.975) * se,
           CI_high = mean + qnorm(0.975) * se) %>%
    select(-n,-sd, -se) %>%
    left_join(original_values, by = "stat") %>%
    mutate_at(vars(-stat), round, 3)%>%
    as.data.frame()
  )
  res_df <- na.omit(res_df)
  return(list(pgls=pgls, r2pred=r2pred, stats=res_sum, raw=res_df))
}



pgls_for_sample_size_se2 <- function(formula, tree, df, SScol, min_species = 10, method = "ML", model = "BM",
                                     treename = "Short", tree_used = NULL, se=NULL, startval=NULL, UseFirstEstimate = FALSE) {
  ####
  ## DEBUG
  # formula = total_sprmLen_um_log_mean ~ Rec_binary
  # tree = comp_df[[2]]
  # df = comp_df[[1]]
  # SScol = "N_sperm_length"
  # min_species = 40
  # method = "ML"
  # model = "OU"
  # treename = "Short"
  # se = "total_sprmLen_um_log_stderr"
  # startval=NULL
  # UseFirstEstimate = TRUE
  ## END DEBUG
  ####
  
  require(dplyr)
  require(phytools)
  require(rr2)
  correlation = ifelse(model == "BM","corBrownian",
                       ifelse(model == "OU","corMartins",
                              ifelse(model == "lambda","corPagel","err")
                       ))
  correlationgls = ifelse(model == "BM","corBrownian(1, phy = tree)",
                          ifelse(model == "OU","corMartins(1, phy = tree, fixed = FALSE)",
                                 ifelse(model == "lambda","corPagel(1, phy = tree, fixed = FALSE)","err")
                          ))
  # Sometimes pgls.SEy fails so then we can use gls this is the switch variable for this
  USEGLS = FALSE
  
  ## TODO add to loo
  interval <- if (model == "lambda") c(0,1) else c(0,100)
  
  if (correlation == "err") {stop("Don't recognize evolutionary model. must be BM,OU or lambda")}
  if (class(formula) != "formula") 
    stop("formula must be class 'formula'")
  
  # incooperate SE in y
  if (!is.null(se)){
    fullse = setNames(df[, se], df[,treename])
  } else {
    fullse = NULL
  }
  pgls <- tryCatch(
    {if(!is.null(startval)){
      pgls.SEy_startval(formula, data = df, tree = tree,
                        corClass =eval(parse(text=correlation)), se=fullse,
                        method = method, startval = startval, fixed = UseFirstEstimate,interval = interval)
    } else {
      pgls.SEy(formula, data = df, tree = tree,
               corClass =eval(parse(text=correlation)), se=fullse, method = method, interval = interval,
               fixed = FALSE)
    }
    },
    error=function(cond){
      message("\npgls.SEy throws an error we will use gls() in this analysis!")
      message(paste("ERROR:", cond, "\n"))
      # switch to gls
      return(
        tryCatch(
          {gls(model = formula,
               correlation = eval(parse(text=correlationgls)),
               data = df, method = method, interval = interval)},
          error=function(cond) {
            message("gls also did not work!")
            stop("aborting run!")
          }
        )
      )},
    warning=function(cond){
    })
  # turned the switch off!
  if(sum(grepl("SEy", pgls$call)) == 0) {USEGLS = TRUE; print("USING GLS")} 
  if (model %in% c("OU", "lambda")) {
    "estimated parameter of the model, will us this in subsampling"
    STARTVAL = pgls$modelStruct$corStruct
    print(STARTVAL)
    #print(paste("half-life:", round(log(2)/STARTVAL,2)))
  }
  
  r2pred = round(rr2::R2.pred(pgls),3)
  # consistent subsetting
  # Ns = sort(unique(df[,SScol]))
  Ns = 1:10
  s_base <- summary(pgls)
  #print(s_base)
  nVars = nrow(s_base$tTable)
  # prepare a pgls_table
  res_df <- as.data.frame(matrix(NA, ncol = 2+(4*nVars), nrow = length(Ns)))
  s <- as.data.frame(summary(pgls)$tTable)
  s$par <- row.names(s); s$Min_sample_size <- min(Ns); s$N_species_in_model <- length(tree$tip.label)
  s_res <- reshape(s, idvar = c("Min_sample_size", "N_species_in_model"), timevar = "par", direction = "wide")
  nam <- gsub("\\(|\\)","",names(s_res))
  nam <- gsub("Std.Error","SE",nam)
  nam <- gsub("(.+)\\.(.+)","\\2:\\1",nam)
  names(res_df) <- nam
  
  original_values <- as.data.frame(t(s_res[2:length(s_res)]))
  original_values$par <- as.character(nam[2:length(nam)])
  names(original_values) <- c("full_model", "stat")
  original_values$stat <- as.factor(original_values$stat)
  ### Sample size resampling
  prg <- txtProgressBar(0,length(Ns), initial = 0, style = 2)
  for (i in 1:length(Ns)){
    setTxtProgressBar(prg, i)
    # remove everything with lower sample size
    index = df[,SScol] >= Ns[[i]]
    r_df <- df[index,]
    # check that with this filtering we have enough sp and still all factors
    if(nrow(r_df) < min_species){
      res <- c(Ns[[i]], nrow(r_df), rep("SMALL N",(4*nVars)))
    } else {
      rownames(r_df) <- r_df[,treename]
      if(sum(index == FALSE) != 0) {
        r_tree <- drop.tip(tree, tip = df[!index,treename])#1
      }else{
        r_tree <- tree 
      }
      
      if (!is.null(se)){
        r_se = setNames(r_df[, se], r_df[,treename])
      } else {
        r_se = NULL
      }
      if(USEGLS == FALSE) {
        pgls_r <- tryCatch(
          {if(UseFirstEstimate & model %in% c("OU", "lambda")) {
            # print(paste0("r_df:",nrow(r_df)))
            # print(paste0("length tree:",length(r_tree$tip.label)))
            # print(paste0("startval:", STARTVAL))
            # print(paste0("r_se:",length(r_se)))
            pgls.SEy_startval(formula, data = r_df, tree = r_tree,
                              corClass = eval(parse(text=correlation)),
                              se=r_se, method = method, fixed = UseFirstEstimate, startval = STARTVAL, interval = interval)
          } else{
            pgls.SEy(formula, data = r_df, tree = r_tree,
                     corClass =eval(parse(text=correlation)), se=r_se, method = method, interval = interval)
          }
          },
          error=function(cond){
            message("\npgls.SEy did not converge, trying gls")
            message(paste("ERROR:", cond))
            return(
              tryCatch(
                {gls(model = formula,
                     correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
                     data = r_df, method = method, interval = interval)},
                error=function(cond) {
                  message("\ngls also did not work! Skipping this iteration.")
                  return(NA)
                }
              ))},
          warning=function(cond){
          })
      } else {
        pgls_r <- tryCatch(
          {gls(model = formula,
               correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
               data = r_df, method = method, interval = interval)},
          error=function(cond) {
            message("\ngls did not work! Skipping this iteration.")
            return(NA)
          })
      }
      
      if(is.na(pgls_r[1])){
        res <- c(Ns[[i]], nrow(r_df), rep("MODEL FAIL",(4*nVars)))
      } else {
        s <- as.data.frame(summary(pgls_r)$tTable)
        s$par <- row.names(s)
        s$Min_sample_size <- Ns[[i]]
        s$N_species_in_model <- nrow(r_df)
        res <- reshape(s, idvar = c("Min_sample_size", "N_species_in_model"), timevar = "par", direction = "wide")
        # reduction can lead to loss of factors. we catch this here and print a warning
        if(ncol(res) != length(nam)){
          res <- c(Ns[[i]], nrow(r_df), rep("FACTOR MISSING",(4*nVars)))
        }
        names(res) <- nam
      }
    }
    res_df[i,] <- res
  }
  close(prg)
  suppressWarnings(
    res_sum <- res_df[,2:length(res_df)] %>%
      select(-N_species_in_model) %>%
      gather(key = stat, value = val) %>%
      mutate(val = as.numeric(val),
             stat = factor(stat, levels = nam)) %>%
      group_by(stat)%>%
      summarise_all( funs(mean = mean(., na.rm =TRUE), 
                          sd = sd(., na.rm =TRUE),
                          n = n())) %>%
      mutate(se = sd / sqrt(n),
             CI_low = mean - qnorm(0.975) * se,
             CI_high = mean + qnorm(0.975) * se) %>%
      select(-n,-sd, -se) %>%
      left_join(original_values, by = "stat") %>%
      mutate_at(vars(-stat), round, 3)%>%
      as.data.frame()
  )
  res_df <- na.omit(res_df)
  
  # add some input information and tree size to the output tables
  
  if(is.null(tree_used)) {
    tree_used = "missing"
  } 
  
  res_sum$tree <- tree_used
  res_sum$r2pred <- r2pred
  res_sum$N <- length(tree$tip.label)
  res_sum$model <- model

  res_df$tree <- tree_used
  res_df$model <- model
  
  return(list(pgls=pgls, r2pred=r2pred, stats=res_sum, raw=res_df))
}



leave_one_out_se2 <- function(formula, tree, df, method = "ML", model = "BM", treename = "Short", se=NULL,
                              startval = NULL,  UseFirstEstimate=FALSE, tree_used = NULL) {
  require(phytools)
  require(dplyr)
  require(rr2)
  ## DEBUG
  # comp_df <- species_df %>%
  #   select(Short, Rec_binary, bristle_binary, total_sprmLen_um_log_mean,
  #          total_sprmLen_um_log_stderr, N_sperm_length,
  #          bdSize_um_logsqrt_mean
  #   ) %>%
  #   mutate(Rec_binary = factor(Rec_binary, levels=c(1,0)) )%>%
  #   prepare_phyl_pca(tree = trees[[i]], treename = "Short")
  # 
  # UseFirstEstimate = TRUE
  # startval = 0.73
  # formula = total_sprmLen_um_log_mean ~ Rec_binary
  # se = "total_sprmLen_um_log_stderr"
  # tree = comp_df[[2]]
  # df = comp_df[[1]]
  # method = "ML"
  # model = "lambda"
  ## /DEBUG
  
  influ_threshold = 5
  correlation = ifelse(model == "BM","corBrownian",
                       ifelse(model == "OU","corMartins",
                              ifelse(model == "lambda","corPagel","err")
                       ))
  
  correlationgls = ifelse(model == "BM","corBrownian(1, phy = tree)",
                          ifelse(model == "OU","corMartins(1, phy = tree, fixed = FALSE)",
                                 ifelse(model == "lambda","corPagel(1, phy = tree, fixed = FALSE)","err")
                          ))
  
  
  if (correlation == "err") {stop("Don't recognize evolutionary model. must be BM,OU or lambda")}
  if (class(formula) != "formula") 
    stop("formula must be class 'formula'")
  # Sometimes pgls.SEy fails so then we can use gls this is the switch variable for this
  USEGLS = FALSE
  
  
  # incooperate SE in y
  if (!is.null(se)){
    fullse = setNames(df[, se], df[,treename])
  } else {
    fullse = NULL
  }
  pgls <- tryCatch(
    {if(!is.null(startval)){
      pgls.SEy_startval(formula, data = df, tree = tree,
                        corClass =eval(parse(text=correlation)), se=fullse,
                        method = method, startval = startval, fixed = UseFirstEstimate)
    } else {
      pgls.SEy(formula, data = df, tree = tree,
               corClass =eval(parse(text=correlation)), se=fullse, method = method)
    }},
    error=function(cond){
      message("\npgls.SEy throws an error we will use gls() in this analysis!")
      message(paste("ERROR:", cond))
      return(
        tryCatch(
          {gls(model = formula,
               correlation = eval(parse(text=correlationgls)),
               data = df, method = method)},
          error=function(cond) {
            message("gls also did not work!")
            stop("aborting run!")
          }
        )
      )},
    warning=function(cond){
    })
  # the gls call does not have a weights variable
  # liam changed phytools and the old check does not work anymore
  if(sum(grepl("SEy", pgls$call)) == 0) {USEGLS = FALSE} 
  print(pgls)
  print(sum(grepl("SEy", pgls$call)))
  if (model %in% c("OU", "lambda")) {
    "estimated parameter of the model, will us this in subsampling"
    STARTVAL = pgls$modelStruct$corStruct
    print(STARTVAL)
  }
  s_base <- summary(pgls)
  print(s_base)
  nVars = nrow(s_base$tTable)
  print(s_base$tTable)
  r2pred = round(rr2::R2.pred(pgls),3)
  # pgls_table
  res_df <- as.data.frame(matrix(NA, ncol = 1+(4*nVars), nrow = nrow(df)))
  s <- as.data.frame(summary(pgls)$tTable)
  s$par <- row.names(s); s$sp_removed <- "None"
  s_res <- reshape(s, idvar = "sp_removed", timevar = "par", direction = "wide")
  nam <- gsub("\\(|\\)","",names(s_res))
  nam <- gsub("Std.Error","SE",nam)
  nam <- gsub("(.+)\\.(.+)","\\2:\\1",nam)
  names(res_df) <- nam
  original_values <- as.data.frame(t(s_res[2:length(s_res)]))
  original_values$par <- as.character(nam[2:length(nam)])
  names(original_values) <- c("full_model", "par")
  prg <- txtProgressBar(0,length(df[,treename]), initial = 0, style = 2)
  # Downsample data frame
  for (i in 1:length(df[,treename])){
    setTxtProgressBar(prg, i)
    index = !df[,treename] == df[i,treename]
    r_df <- df[index,]
    rownames(r_df) <- r_df[,treename]
    r_tree <- drop.tip(tree, tip = as.character(df[i,treename]))
    if (!is.null(se)){
      r_se = setNames(r_df[, se], r_df[,treename])
    } else {
      r_se = NULL
    }
    # print(paste0("r_df:",nrow(r_df)))
    # print(paste0("length tree:",length(r_tree$tip.label)))
    # print(paste0("startval:", STARTVAL))
    # print(paste0("r_se:",length(r_se)))
    if(USEGLS == FALSE) {
      pgls_r <- tryCatch(
        {if(UseFirstEstimate & model %in% c("OU", "lambda")) {
          pgls.SEy_startval(formula, data = r_df, tree = r_tree,
                            corClass = eval(parse(text=correlation)), se=r_se, method = method,
                            fixed = UseFirstEstimate, startval = STARTVAL)
        } else{
          pgls.SEy(formula, data = r_df, tree = r_tree,
                   corClass =eval(parse(text=correlation)), se=r_se, method = method)
        }
        },
        error=function(cond){
          message("\npgls.SEy did not converge, trying gls")
          message(paste("ERROR:", cond))
          return(
            tryCatch(
              {gls(model = formula,
                   correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
                   data = r_df, method = method)},
              error=function(cond) {
                message("\ngls also did not work! Skipping this iteration.")
                return(NA)
              }
            ))},
        warning=function(cond){
        })
    } else {
      pgls_r <- tryCatch(
        {gls(model = formula,
             correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
             data = r_df, method = method)},
        error=function(cond) {
          message("\ngls did not work! Skipping this iteration.")
          return(NA)
        })
    }
    if(is.na(pgls_r[1])){
      res <- c(df[i,treename], rep("MODEL FAIL",(4*nVars)))
    } else {
      s <- as.data.frame(summary(pgls_r)$tTable)
      s$par <- row.names(s)
      s$sample <- df[i,treename]
      res <- reshape(s, idvar = "sample", timevar = "par", direction = "wide")
      names(res) <- nam
    }
    res_df[i,] <- res
  }
  close(prg)
  # Now calculate statistics
  res_sum <- suppressWarnings(res_df[,2:ncol(res_df)] %>%
    gather(key = par, value = val) %>%
    mutate(val = as.numeric(val),
           par = factor(par, levels = nam)) %>%
    group_by(par)%>%
    summarise_all(funs(mean = mean(., na.rm =TRUE), 
                        sd = sd(., na.rm =TRUE),
                        n = n())) %>%
    mutate(se = sd / sqrt(n),
           CI_low = mean - qnorm(0.975) * se,
           CI_high = mean + qnorm(0.975) * se,
           tree = ifelse(is.null(tree_used),"missing",tree_used),
           model = model) %>%
    as.data.frame() %>%
    left_join(original_values, by="par") %>%
    select(par, tree, model, full_model, mean, CI_low, CI_high) %>%
    mutate_if(is.numeric, round,3))
  
  # finding outliers and how strong their influence is
  res_long <- res_df %>%
    gather(key = par, value = val, -sp_removed) %>%
    mutate(val = as.numeric(val)) %>%
    left_join(original_values, by="par") %>% 
    # only do the influence measure for the Values
    filter(grepl("Value", par)) %>% as.data.frame()
  
  res_ls <- split(res_long, res_long$par)
  res_ls <- lapply(res_ls, find_influ)
  outliers <- do.call(rbind, res_ls) %>%
    select(sp_removed, par, full_model, everything()) %>%
    filter(pcnt_D >= influ_threshold) %>%
    arrange(par, desc(pcnt_D)) %>%
    mutate_if(is.numeric, round, 3) 
  
  if(nrow(outliers)>0) {
    outliers$tree <- ifelse(is.null(tree_used),"missing",tree_used)
    outliers$model <- model
    }
  res_df$tree <- ifelse(is.null(tree_used),"missing",tree_used)
  res_df$model <- model

  return(list(pgls=pgls, r2pred=r2pred, stats=res_sum, raw=mutate_if(res_df, is.numeric, round, 3), outliers=outliers))
}


leave_one_out_se <- function(formula, tree, df, method = "ML", model = "BM", treename = "Short", se=NULL, UseFirstEstimate=FALSE) {
  require(phytools)
  require(dplyr)
  require(rr2)
  influ_threshold = 50
  correlation = ifelse(model == "BM","corBrownian",
                       ifelse(model == "OU","corMartins",
                              ifelse(model == "lambda","corPagel","err")
                       ))
  
  correlationgls = ifelse(model == "BM","corBrownian(1, phy = tree)",
                   ifelse(model == "OU","corMartins(1, phy = tree, fixed = FALSE)",
                   ifelse(model == "lambda","corPagel(1, phy = tree, fixed = FALSE)","err")
                                       ))
  
  
  if (correlation == "err") {stop("Don't recognize evolutionary model. must be BM,OU or lambda")}
  if (class(formula) != "formula") 
    stop("formula must be class 'formula'")
  # Sometimes pgls.SEy fails so then we can use gls this is the switch variable for this
  USEGLS = FALSE
  
  
  # incooperate SE in y
  if (!is.null(se)){
    fullse = setNames(df[, se], df[,treename])
  } else {
    fullse = NULL
  }
  pgls <- tryCatch(
    {pgls.SEy(formula, data = df, tree = tree,
              corClass =eval(parse(text=correlation)), se=fullse, method = method)},
    error=function(cond){
      message("\npgls.SEy throws an error we will use gls() in this analysis!")
      message(paste("ERROR:", cond))
      return(
        tryCatch(
          {gls(model = formula,
                      correlation = eval(parse(text=correlationgls)),
                      data = df, method = method)},
          error=function(cond) {
          message("gls also did not work!")
            stop("aborting run!")
                   }
          )
        )},
    warning=function(cond){
    })
  # the gls call does not have a weights variable
  if(is.null(pgls$call$weights)) {USEGLS = TRUE} 
  if (model %in% c("OU", "lambda")) {
    "estimated parameter of the model, will us this in subsampling"
    STARTVAL = pgls$pgls$modelStruct$corStruct
  }
  s_base <- summary(pgls)
  print(s_base)
  nVars = nrow(s_base$tTable)
  r2pred = round(rr2::R2.pred(pgls),3)
  # pgls_table
  res_df <- as.data.frame(matrix(NA, ncol = 1+(4*nVars), nrow = nrow(df)))
  s <- as.data.frame(summary(pgls)$tTable)
  s$par <- row.names(s); s$sp_removed <- "None"
  s_res <- reshape(s, idvar = "sp_removed", timevar = "par", direction = "wide")
  nam <- gsub("\\(|\\)","",names(s_res))
  nam <- gsub("Std.Error","SE",nam)
  nam <- gsub("(.+)\\.(.+)","\\2:\\1",nam)
  names(res_df) <- nam
  original_values <- as.data.frame(t(s_res[2:length(s_res)]))
  original_values$par <- as.character(nam[2:length(nam)])
  names(original_values) <- c("full_model", "par")
  prg <- txtProgressBar(0,length(df[,treename]), initial = 0, style = 2)
  for (i in 1:length(df[,treename])){
    setTxtProgressBar(prg, i)
    index = !df[,treename] == df[i,treename]
    r_df <- df[index,]
    rownames(r_df) <- r_df[,treename]
    r_tree <- drop.tip(tree, tip = as.character(df[i,treename]))
    if (!is.null(se)){
      r_se = setNames(r_df[, se], r_df[,treename])
    } else {
      r_se = NULL
    }
    
    if(USEGLS == FALSE) {
      pgls_r <- tryCatch(
        {pgls.SEy(formula, data = r_df, tree = r_tree,
                  corClass =eval(parse(text=correlation)), se=r_se, method = method)},
        error=function(cond){
          message("\npgls.SEy did not converge, trying gls")
          message(paste("ERROR:", cond))
          return(
            tryCatch(
              {gls(model = formula,
                   correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
                   data = r_df, method = method)},
              error=function(cond) {
                message("\ngls also did not work! Skipping this iteration.")
                return(NA)
              }
            ))},
        warning=function(cond){
        })
    } else {
      pgls_r <- tryCatch(
        {gls(model = formula,
             correlation = eval(parse(text=gsub(pattern = "tree", replacement = "r_tree", correlationgls))),
             data = r_df, method = method)},
        error=function(cond) {
          message("\ngls did not work! Skipping this iteration.")
          return(NA)
        })
    }
    if(is.na(pgls_r[1])){
      res <- c(df[i,treename], rep("MODEL FAIL",(4*nVars)))
    } else {
      s <- as.data.frame(summary(pgls_r)$tTable)
      s$par <- row.names(s)
      s$sample <- df[i,treename]
      res <- reshape(s, idvar = "sample", timevar = "par", direction = "wide")
      names(res) <- nam
    }
    res_df[i,] <- res
  }
  close(prg)
  # Now calculate statistics
  res_sum <- res_df[,2:ncol(res_df)] %>%
    gather(key = par, value = val) %>%
    mutate(val = as.numeric(val),
           par = factor(par, levels = nam)) %>%
    group_by(par)%>%
    summarise_all( funs(mean = mean(., na.rm =TRUE), 
                        sd = sd(., na.rm =TRUE),
                        n = n())) %>%
    mutate(se = sd / sqrt(n),
           CI_low = mean - qnorm(0.975) * se,
           CI_high = mean + qnorm(0.975) * se 
    ) %>%
    select(-n,-sd, -se) %>%
    mutate_at(vars(-par), round, 3)%>%
    as.data.frame() %>%
    left_join(original_values) %>%
    mutate_if(is.numeric, round,3)
  
  # finding outliers and how strong their influence is
  res_long <- res_df %>%
    gather(key = par, value = val, -sp_removed) %>%
    mutate(val = as.numeric(val)) %>%
    left_join(original_values) %>% 
    # only do the influence measure for the Values
    filter(grepl("Value", par)) %>%
    as.data.frame()
  
  res_ls <- split(res_long, res_long$par)
  res_ls <- lapply(res_ls, find_influ)
  outliers <- do.call(rbind, res_ls) %>%
    select(-full_model) %>%
    filter(pcnt_D >= influ_threshold) %>%
    arrange(desc(pcnt_D)) %>%
    mutate_if(is.numeric, round, 3)
  
  return(list(pgls=pgls, r2pred=r2pred, stats=res_sum, raw=mutate_if(res_df, is.numeric, round, 3), outliers=outliers))
}





robust_pgls <- function(trees, df, min_species = 40, outnumber, response, response_err, response_N, Lh_model, evo_model, response_plot_text, diagdir, predictor, predictor_plot_text, predictor_tick_labels = NA) {
  ltre <- length(trees)
  simple_pgls_res_ls <- list()
  main_pgls<- list()
  lool <- list()
  lool_stat = list(); lool_raw = list(); lool_outlier = list()
  pgls_stat = list(); pgls_raw = list()
  
  pdf(paste0("fig/",outnumber,"_", response, "_", predictor,"_", evo_model, "_pgls.pdf"))
  for (i in 1:ltre){
    comp_df <- df[, c("Short", predictor, response, response_err, response_N)]  %>%
      prepare_phyl_pca(tree = trees[[i]], treename = "Short")
    
    pgls <- pgls_for_sample_size_se(as.formula(paste0(response," ~ ", predictor)),
                                    SScol = response_N, se = response_err,
                                    tree = comp_df[[2]], df = comp_df[[1]],
                                    min_species = min_species, method = Lh_model, model = evo_model)
    
    loo <- leave_one_out_se(as.formula(paste0(response," ~ ", predictor)),
                            tree = comp_df[[2]], df = comp_df[[1]],
                            method = Lh_model, model = evo_model)
    treename <- names(trees)[[i]]
    main_pgls <- append(main_pgls, setNames(list(pgls), names(trees)[[i]]))
    lool <- append(lool, setNames(list(loo), names(trees)[[i]]))
    main_pgls[[i]]$stats$tree <- treename
    main_pgls[[i]]$stats$r2pred <- main_pgls[[i]]$r2pred
    main_pgls[[i]]$stats$N <- main_pgls[[i]]$pgls$dims$N
    main_pgls[[i]]$raw$tree <- treename
    lool[[i]]$stats$tree <- treename
    if(nrow(lool[[i]]$outliers)>0) {
      lool[[i]]$outliers$tree <- treename
    }
    lool[[i]]$raw$tree <- treename
    lool_stat <- append(lool_stat, list(lool[[c(i,3)]]))
    lool_raw <- append(lool_raw, list(lool[[c(i,4)]]))
    lool_outlier <- append(lool_outlier, list(lool[[c(i,5)]]))
    pgls_stat <- append(pgls_stat, list(main_pgls[[c(i,3)]]))
    pgls_raw <- append(pgls_raw, list(main_pgls[[c(i,4)]]))
    
    if(is.na(predictor_tick_labels)) {
      p <- ggplot(comp_df[[1]], aes(x = eval(parse(text=predictor)), y = eval(parse(text=response)))) +
        xlab(predictor_plot_text)+
        ylab(response_plot_text) +
        geom_point(size = 2) +
        geom_abline(intercept = main_pgls[[i]]$pgls$coefficients[1], slope = main_pgls[[i]]$pgls$coefficients[2],
                    color = "black", size= 1) +
        ggtitle(treename)
    } else {
      p <- ggplot(comp_df[[1]], aes(x = eval(parse(text=predictor)), y = eval(parse(text=response)))) +
        xlab(predictor_plot_text)+
        ylab(response_plot_text) +
        geom_point(size = 2) +
        geom_abline(intercept = main_pgls[[i]]$pgls$coefficients[1], slope = main_pgls[[i]]$pgls$coefficients[2],
                    color = "black", size= 1) +
        ggtitle(treename)  +
        scale_x_discrete(labels = predictor_tick_labels)  
    }
    print(p)
  }
  dev.off()
  
  lool_stat <- do.call(rbind, lool_stat)
  lool_raw <- do.call(rbind, lool_raw)
  lool_outlier <- do.call(rbind, lool_outlier)
  pgls_stat <- do.call(rbind, pgls_stat)
  pgls_raw <- do.call(rbind, pgls_raw)
  
  print(pgls_stat)
 
  
  outtables <- list(lool_stat, lool_raw, lool_outlier, pgls_stat, pgls_raw)
  outnames <- c("leave_one_out_stats", "leave_one_out_raw", "leave_one_out_outlier", 
                "pgls_stat", "pgls_sample_size")
  
  for (i in 1:length(outtables)) {
    write.csv(outtables[[i]],
              file = paste0(diagdir,outnumber, "_", response, "_", predictor, "_", evo_model, "_",
                            outnames[[i]], ".csv"),row.names = FALSE)
  }
  return(format_pgls_publication(pgls_stat, predictor = predictor))
}

  

prepare_trait_mode <- function(df,tree, varcol, SEcol, treename="Short") {
  library(dplyr)
  # produces a comparative dataset and gives reasons for omissions
  print("Producing comparative dataset and corresponding pruned tree.")
  df_c <- df[df[,treename] %in% tree$tip.label,c(treename, varcol, SEcol)]
  sum_taxa <- as.character(df_c[,treename][is.na(df_c[,varcol])])
  
  print(paste(length(sum_taxa),"taxa will be removed:"))
  print(sum_taxa)
  # replace SE if it is NA
  df_c <- df_c[!is.na(df_c[,varcol]),]
  meanSE <- mean(df_c[,SEcol],na.rm=TRUE)
  df_c[,SEcol] <- ifelse(is.na(df_c[,SEcol]), meanSE, df_c[,SEcol])
  to_drop <- tree$tip.label[! tree$tip.label %in% df_c[,treename]]
  t <- drop.tip(tree, tip = to_drop)
  print(paste("pruned data has:",length(t$tip.label), "tips"))
  # create named vectors
  x <- setNames(df_c[,varcol], df_c[,treename])
  SEx <- setNames(df_c[,SEcol], df_c[,treename])
  return(list(x, SEx, t))
}



format_pgls_publication <- function(pgls_stat_object, predictor) {
  require(dplyr)
  summary_pgls_stat <- pgls_stat_object %>%
    select(full_model, mean, CI_low, CI_high, everything()) %>%
    separate(col = stat, into = c("parameter","statistic"), sep = ":") %>%
    reshape(idvar = c("parameter", "tree","r2pred", "N", "model"), timevar = "statistic", direction = "wide") %>%
    mutate(predictor = predictor,
           parameter = ifelse(parameter == "Intercept", "intercept", "slope")) %>%
    select(predictor, tree, model, N, r2pred, everything())

  eps=0.0001 # precision of the p-values
  format_stat <- summary_pgls_stat %>%
    mutate(Value = paste0(full_model.Value,
                          " (", mean.Value, ", ",round(abs(mean.Value - CI_low.Value),3),")"),
           SE = paste0(full_model.SE,
                       " (", mean.SE, ", ",round(abs(mean.SE - CI_low.SE),3),")"),
           `t-value` = paste0(`full_model.t-value`,
                              " (", `mean.t-value`, ", ",round(abs(`mean.t-value` - `CI_low.t-value`),3),")"),
           `p-value` = paste0(format.pval(`full_model.p-value`,eps=eps, scientific =1),
                              " (", format.pval(`mean.p-value`,eps=eps, scientific =1), ", ",round(abs(`mean.p-value` - `CI_low.p-value`),3),")")
    ) %>%
    select(predictor,tree, model, N,r2pred,parameter, Value,SE,`t-value`, `p-value`)
  return(list(summary_pgls_stat, format_stat))
}


pgls.SEy_startval <- function (model, fixed=FALSE, startval=1, data, corClass = corBrownian, tree = tree, se = NULL, 
          method = c("REML", "ML"), interval = c(0, 1000), ...) 
{
  Call <- match.call()
  corfunc <- corClass
  data <- data[tree$tip.label, ]
  if (is.null(se)) 
    se <- setNames(rep(0, Ntip(tree)), tree$tip.label)
  lk <- function(sig2e, data, tree, model, ve, corfunc) {
    tree$edge.length <- tree$edge.length * sig2e
    ii <- sapply(1:Ntip(tree), function(x, e) which(e == 
                                                      x), e = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + ve[tree$tip.label]
    vf <- diag(vcv(tree))
    w <- varFixed(~vf)
    if (fixed) {
      COR <- corfunc(startval, tree, fixed = TRUE, ...)
    } else {
      COR <- corfunc(startval, tree, ...)
    }
    fit <- gls(model, data = cbind(data, vf), correlation = COR, 
               method = method, weights = w)
    -logLik(fit)
  }
  fit <- optimize(lk, interval = interval, data = data, tree = tree, 
                  model = model, ve = se^2, corfunc = corfunc)
  tree$edge.length <- tree$edge.length * fit$minimum
  ii <- sapply(1:Ntip(tree), function(x, e) which(e == x), 
               e = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + se[tree$tip.label]^2
  vf <- diag(vcv(tree))
  w <- varFixed(~vf)
  if (fixed) {
    COR <- corfunc(startval, tree, fixed = TRUE, ...)
    obj <- gls(model, data = cbind(data, vf), correlation = COR, 
        method = method, weights = w)
  } else {
    obj <- gls(model, data = cbind(data, vf),
               correlation = corfunc(startval, tree), weights = w, method = method)
  }
  obj$call <- Call
  obj
}


mkBound <- function(x) {
  # the Uncertainty_FPK function has problems with linear algebra if 
  # the bounds are too extreme. This returns bounds that do include 0
  # but don't get too extreme.
  if (x < 0){
    high = 1
    low = x + x
  } else {
    low = -1
    high = x + x
  }
  return(c(low, high))
}
