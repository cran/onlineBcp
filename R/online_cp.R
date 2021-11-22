#'   Online change point detection algorithm
#'       for normally distributed data.
#'

#' @param x the normalized data
#' @param theta the probability of occurrence of a change point, default 0.9
#' @param alpha the hyperparameter of posterior distribution, default 1.0
#' @param beta the hyperparameter of posterior distribution, default 1.0
#' @param th_cp threshold level for the posterior distribution of change point, default 0.5
#   alpha: 1                                                      #
#   beta: (0.5, 1)                                                #
#   th_cp: threshold level for the posterior distribution of      #
#          change point                                           #
#   recommended value 0.5                                         #
#                                                #
###################################################################

#' @return An object of the BayesCP class
#' @export
online_cp <- function(x, theta = 0.9, alpha = 1, beta = 1, th_cp = 0.5) {
  n=length(x)

  # theta=0.90; 		    		        # theta=(1-H)-> H=0.10 hazard function is taken to be geometric or exponential with theta
  # alpha=1; beta=1				          # these are the values of pior distribution's parameters
  p_ct=array(0,dim=c(n));     		# i=ct=0,1,2,...,t-1; p_ct: the probability of ct
  # t=1,2,...,n
  p_cs=array(0,dim=c(n));     		# s=t+1; j=0,1,2,...,t; p_cs: the probability of ct+1
  max_p=array(0, dim=c(n,3)); 		# keeps max probability and change point locations
  t=1;
  p_ct[t]=1; j=0; tau=0

  while (t<n){   								  # 1st loop

    j=tau; 								        # it is required for 2nd loop "while (j<=t)"

    if(t %% 100 == 0) cat("processing data point ",t,"\n")

    w1=array(0,dim=c(t))          #

    while (j<=t){   							# 2nd loop

      #cat("j=",j,"\n")
      y1=sum(x[(j+1):(t+1)])
      y2=sum(x[(j+1):t])
      y3=sum(x[(t+1):(t+1)])
      ys1=sum(x[(j+1):(t+1)]^2)
      ys2=sum(x[(j+1):t]^2)
      ys3=sum(x[(t+1):(t+1)]^2)
      a1=ys1-((y1^2)/(t-j+2))
      a2=ys2-((y2^2)/(t-j+1))
      a3=ys3/2

      us1=(t-j+2*alpha)/2
      us2=(t-j+2*alpha+1)/2
      us3=(2*alpha+1)/2
      if (j<t){
        w1[j+1]=(1/sqrt(pi))*(gamma(us2)/gamma(us1))*(sqrt(t-j+1)/sqrt(t-j+2))*((a2+2*beta)^us1)/((a1+2*beta)^us2);

        if (is.nan(w1[j+1])==TRUE | is.infinite(w1[j+1])==TRUE) {w1[j+1]=0}
        p_cs[j+1]=w1[j+1]*(theta)*p_ct[j+1]
      }
      else {
        w2=((2^alpha)*gamma(us3))/(sqrt(2*pi)*(a3+beta)^us3)
        sum_p=sum(p_ct[1:t])
        p_cs[j+1]=w2*(1-theta)*sum_p
        # cat("p_cs", p_cs[j+1], "\n")
      }
      j=j+1              # it is required for 2nd loop
    } 		               # end of 2nd loop

    p_ct=p_cs/sum(p_cs)

    max_p[t,1]=max(p_ct)

    # max_p[t,2]=which.max(p_ct)-1       # orginal code

    if ( max(p_ct) > th_cp ) {

      max_p[t,2]=which.max(p_ct)-1       # new code

      tau=max_p[t,2];

      p_cs=array(0,dim=c(n));  # it is required 2nd loop and to decrease comutation cost we give a threshold level and it is required 2nd loop

      # cat("tau=",tau,"\n");

    }

    ## run length is added

    if ( t==1) { max_p[t,3]=0;
    } else if (t>1 & max_p[t-1,2]==max_p[t,2]) { max_p[t,3]=max_p[t-1,3]+1
    } else { max_p[t,3]=0}

    t=t+1

  } 	# end of 1st loop

  bcp <- list(x=x, max_p = max_p, parameters = c(theta = theta, alpha = alpha, beta = beta, th_cp = th_cp),
              series_length = length(x), result = NULL)
  class(bcp) <-  "BayesCP"

  return(bcp)
}




#' Plot BayesCP object
#'
#' @param x the BayesCP class object to be plotted
#' @param xlab the default x-axis label, default "Index"
#' @param ylab the default y-axis label, default "x"
#' @param ... the plotting parameters passed to plot()
#' @return No return value, called for side effects
#' @export
plot.BayesCP <-  function(x, xlab = "Index", ylab = "x", ...) {
  if(is.null(x$result)) {
    print("No result to plot yet. Please call summary() first")
  }
  else {
    bcp <- x ## change x to bcp
    plot(seq_along(bcp$x), bcp$x, type="p", cex = .5, pch=16, col="blue", xlab = xlab, ylab = ylab, ...)
    #plot(bcp$x ~ xdata, type="p", cex = .5, xlim=c(1, n), pch=16, col="blue", xlab = xlab, ylab = ylab, ...)
    res <- bcp$result
    segment <- res$segment

    for (i in 1:dim(segment)[1] ) {
      lines(c(segment[i, 1], segment[i, 2]), c(segment[i, 3], segment[i, 3]), col="red", lty=1, lwd=2)
    }
  }



}



#' Summarize BayesCP object
#'
#' @param object the BayesCP class object to be summarized
#' @param norm.test logical value for normality test, default is false
#' @param ... parameters passed to summary()
#' @return An object of BayesCP class with updated summary result
#' @export
#' @examples
#' x <- c(rnorm(10, 0, 1), rnorm(10, 5, 1))
#' bcp <- online_cp(x)
#' summary(bcp)

summary.BayesCP <- function(object, norm.test = FALSE, ...){
  bcp <- object
  x <- bcp$x
  max_p <- bcp$max_p
  th_cp <- bcp$parameters["th_cp"]

  if(is.null(bcp$result) || (length(bcp$series_length) == 1 && dim(bcp$result$segment)[2] == 6 && norm.test) || (length(bcp$series_length) == 1 && dim(bcp$result$segment)[2] == 7 && !norm.test)) {  ## no result yet, do the summary, or change of normality test status

    ## Calculate mean, SD, CI for each segments

    ###################################################################
    #             calculate the means of segments                     #
    ###################################################################
    n <- length(max_p[,1])
    max_ps <- array(0, c(n,4))
    max_ps[,1] <- 1:n
    max_ps[,2:4] <- max_p
    if ( max_ps[1, 3]!=0) { max_ps=rbind(max_ps[1,], max_ps, deparse.level=1); max_ps[1, 2]=th_cp+0.01 }
    if ( max_ps[1, 3]==0 && max_ps[1, 2]<th_cp) { max_ps[1, 2]=th_cp+0.01 }



    #ta_=array(0, c(sum(max_ps[,2]>th_cp),3))
    #ta_=max_ps[max_ps[,2]>th_cp, ]
    #ta_1 <- max_ps[((max_ps[,2] > th_cp) & (max_ps[,4] == 0)),]
    ta_1 <- subset(max_ps, (max_ps[,2] > th_cp) & (max_ps[,4] == 0))

    if ( dim(ta_1)[1] < 2 ) {
      #n_ <- length(ta_1)
      mean_ <- array(0,c(n,3))
      k=1
      i=1

      for(i in 2:(length(ta_1)[1:1])){
        j=ta_1[3]+1
        if ((ta_1[3]!=0) && (ta_1[4]==0) ){ mean_[j:ta_1[3],1]=mean(x[j:ta_1[3]])
        mean_[j:ta_1[3],2]=k ; k=k+1 }
      }
      j=ta_1[3]+1
      mean_[j:n,1]=mean(x[j:n])
      # k=k+1
      mean_[j:n,2]=k
      mean_[1:n,3]=1:n

    } else  {
      n_=length(ta_1[ ,1])
      mean_=array(0,c(n,3))
      k=1
      i=1 ;
      for(i in 2:(dim(ta_1)[1:1])){
        j=ta_1[i-1,3]+1
        if ((ta_1[i,3]!=0) && (ta_1[i,4]==0) ){ mean_[j:ta_1[i,3],1]=mean(x[j:ta_1[i,3]])
        mean_[j:ta_1[i,3],2]=k ; k=k+1 }
      }
      j=ta_1[i,3]+1
      mean_[j:n,1]=mean(x[j:n])
      # k=k+1
      mean_[j:n,2]=k
      mean_[1:n,3]=1:n

    }



    ########################################################################
    # Segments:  begin , end,  Post. Prob., Mean, SD, LL of CI, UL of CI   #
    ########################################################################

    if ( dim(ta_1)[1] < 2 ) {

      # If there is no cp in data

      n_=length(ta_1)

      if(norm.test) {
        norm.res <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x))
        se_result=data.frame(1, n, "no cp", mean(x), sqrt(var(x)),
                             (mean(x)-sqrt(var(x)/n)*qnorm(0.05, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)),
                             (mean(x)+sqrt(var(x)/n)*qnorm(0.05, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)), norm.res$p.value)
        names(se_result)=c("begin", "end","no cp" ,"mean", "SD", "LL of CI", "UL of CI", "normality p-value")
      }
      else {
        se_result=data.frame(1, n, "no cp", mean(x), sqrt(var(x)),
                             (mean(x)-sqrt(var(x)/n)*qnorm(0.05, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)),
                             (mean(x)+sqrt(var(x)/n)*qnorm(0.05, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)))
        names(se_result)=c("begin", "end","no cp" ,"mean", "SD", "LL of CI", "UL of CI")

      }
      changepoint <- "This is no change points."
      res <- list(changepoint = changepoint, segments = se_result)

    } else {

      n_k=table(mean_[,2]); k_=dim(n_k)

      xc=split(mean_[,3], mean_[,2])      # 3th column is indices and splited by segment in 2nd column
      yc=split(x, mean_[,2])              # data is splited by segment in 2nd column
      #zc=split(mean_[,1], mean_[,2])      # mean is splited by segment in 2nd column

      ## merge segments with one observation
      l <- sapply(yc, length)
      index_long_seg = which(l > 1)[1]
      yc_temp <- NULL
      xc_temp <- NULL
      if(index_long_seg > 1) {
        for(i in 1:(index_long_seg-1)) {
          yc_temp <- c(yc_temp, yc[[i]])
          xc_temp <- c(xc_temp, xc[[i]])
        }
        yc[[index_long_seg]] <- c(yc_temp, yc[[index_long_seg]])
        xc[[index_long_seg]] <- c(xc_temp, xc[[index_long_seg]])
      }

      # index_long_seg = 1  ## index for the long segment to be merged with
      start_index <- index_long_seg
      for(i in start_index:k_) {
        if(length(yc[[i]]) > 1)  index_long_seg <-  i
        else {
          yc[[index_long_seg]] <- c(yc[[index_long_seg]], yc[[i]])
          xc[[index_long_seg]] <- c(xc[[index_long_seg]], xc[[i]])
        }

      }
      keep <- sapply(yc, length) > 1
      yc <- yc[keep]
      xc <- xc[keep]

      k_ <- length(yc)

      if(norm.test){
        se_result=array(0, dim=c(k_, 8))
        colnames(se_result)=c("begin", "end","post. prob." ,"mean", "SD", "LL of CI", "UL of CI", "normality p-value")
      }
      else {
        se_result=array(0, dim=c(k_, 7))
        colnames(se_result)=c("begin", "end","post. prob." ,"mean", "SD", "LL of CI", "UL of CI")

      }

      for (i in 1:k_) {

        begin_=xc[[i]][1];
        k_n=length(xc[[i]]); end_=xc[[i]][k_n];
        post_prob=ta_1[i,2];
        seg_mean=mean(yc[[i]]); seg_SD=sd(yc[[i]], na.rm = TRUE);
        LL_of_CI=(mean(yc[[i]])-sqrt(var(yc[[i]])/length(yc[[i]]))*qnorm(0.05, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
        UL_of_CI= (mean(yc[[i]])+sqrt(var(yc[[i]])/length(yc[[i]]))*qnorm(0.05, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
        if(norm.test) {
          norm.res <- ks.test(yc[[i]], "pnorm", mean = seg_mean, sd = seg_SD)
          se_result[i,]=c(begin_, end_, post_prob, seg_mean, seg_SD, LL_of_CI, UL_of_CI, norm.res$p.value)
        }
        else {
          se_result[i,]=c(begin_, end_, post_prob, seg_mean, seg_SD, LL_of_CI, UL_of_CI)

        }

      }
      # take the results for change points out
      changepoint <- data.frame(location = se_result[1:(k_-1), 2], post.prob = se_result[2:k_, 3])
      row.names(changepoint) <- NULL
      #print(changepoint)
      se_result[1, 3] <- NA
      segment <- se_result[, -3]
      res <- list(changepoint = changepoint, segment = segment)

      #print(segment)
    }
  }
  else {
    res <- bcp$result
  }




  ## Output results
  #changePoints <- max_p[(max_p[, 1] > th_cp &  max_p[, 3] == 0 & max_p[, 2]!= 0), , drop = FALSE]

  print("Change points")
  print(res$changepoint)
  print("Segments")
  print(res$segment)

  bcp$result <- res
  return(invisible(bcp))

}



#' Impute missing data
#'
#' @param x the normalized data with missing
#' @param method the imputation method
#'
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats ks.test
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom stats qnorm
#' @importFrom graphics lines
#'
#' @return The vector of imputed data with no missing values
#' @export
imputation <- function(x, method= c("Median", "kNN")) {

  method <- match.arg(method)
  n <- length(x)

  if (anyNA(x) == FALSE) { message( "there is no missing value")

  } else {


    if ( method=="Median") {

      # 4 observations are taken to compute the median of sequence.

      i_na=which(is.na(x)==TRUE, arr.ind = TRUE)

      for (i in length(i_na)) {

        if (i_na[i]<=3) {  x[i_na]=median(na.omit(x[1:5]) )

        } else if (i_na[i]>=(n-3)) { x[i_na]=median(na.omit(x[(n-4):n]))

        } else { x[i_na]=median(na.omit(x[(i_na[i]-2):(i_na[i]+2)])) }

      }

    } else if (method=='kNN') {

      x_=data.frame(x, 1)
      x_na=VIM::kNN(x_, dist_var = c("x", "X1"), k = 2)
      x=as.array(x_na$x)

    }
  }

  return(x)
}


#' Combine two BayesCP objects
#'
#' @param bcp1 the first BayesCP object to be combined
#' @param bcp2 the second BayesCP opbject to be combined
#'
#' @return The combined BayesCP object
#' @export
combine <- function(bcp1, bcp2) {
  # BayesCP has x, max_p, parameters, series_length, and result
  x <- c(bcp1$x, bcp2$x)
  max_p <- rbind(bcp1$max_p, bcp2$max_p)
  parameters <- rbind(bcp1$parameters, bcp2$parameters)
  series_length <- c(bcp1$series_length, bcp2$series_length)
  res1 <- bcp1$result
  res2 <- bcp2$result

  res2$changepoint$location <- res2$changepoint$location + length(bcp1$x)

  changepoint <- rbind(res1$changepoint, res2$changepoint)
  print("here!")

  print(res2$segment[,1])
  res2$segment[,1] <- res2$segment[,1] + length(bcp1$x) # 1st column is the begin index of the segment
  res2$segment[,2] <- res2$segment[,2] + length(bcp1$x) # 1st column is the end index of the segment


  segment <- rbind(res1$segment, res2$segment)
  result <- list(changepoint = changepoint, segment = segment)


  bcp <- list(x = x, max_p = max_p, parameters = parameters, series_length = series_length, result = result)

  class(bcp) <- "BayesCP"
  return(bcp)

}


