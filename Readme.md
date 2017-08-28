# Writing Package

In this tutorial I am trying to show an exampple of writing a package with R. 
This file is part of my project for statistical computing with Prof. Harner.

``` r
#' Implementing One-Way S3 update
#' @seealso  \code{\link[aov]{stat}}
#' @param z is the input in form of list OR data.frame
#' @param ... is the addistion arguments here
#' @seealso \code{\link{print.oneway}}
#' @seealso \code{\link{print.summary.oneway}}
#' @export
oneway <- function(z, ...) UseMethod("oneway")
  #' Calculating the anova parameters
  #' @param z is the input in form of list or data.frame
  #' @param ... is the addition arguments here
  #' @seealso \code{\link{print.oneway}}
  #' @seealso \code{\link{print.summary.oneway}}
  #' @examples
  #' a <- onewayanova:::oneway.default(coagulation)
  #' onewayanova:::print.oneway(a)
  #' @export
  oneway.default <- function(z, ...){
    library(faraway)
    data(coagulation)
    if (is.data.frame(z)==T) z <- with(z, split(z[,1],z[,2])) # Check and Change the input from data.frame to list

  groupNo <- length(z)
  totalObs <- length(unlist(z))
  groupLength <- sapply(z, length)

  groupMean <- sapply(z,mean)
  groupVar <- sapply(z, var)
  grandMean <- mean(unlist(z))

  dfModel <- groupNo - 1
  dfError <- totalObs - groupNo
  dfTotal <- totalObs - 1


  SSB <- sum(groupLength * ((groupMean - grandMean) ^ 2))
  SSW <- sum((groupLength - 1) * groupVar)
  SST <- SSB + SSW

  MSB <- SSB / dfModel
  MSW <- SSW / dfError

  F_stat <- MSB/MSW
  P_value <- 1 - pf(F_stat , dfModel, dfError)

  list(dfb = dfModel, dfw= dfError, dfTotal, ssb = SSB, ssw = SSW, SST, msb = MSB, msw = MSW, F_stat,      P_value,data = z,grpMean = groupMean,grplength= groupLength)
      }



#' oneway.factor combines two columns of factor and values togethere and then call the oneway.default function to do the calculations for anova
#' @param z as a factor input here
#' @param y as values input here
#' @seealso \code{\link{print.oneway}}
#' @seealso \code{\link{print.summary.oneway}}
#' @examples
#' library(faraway)
#' z <- coagulation$diet
#' y <- coagulation$coa
#' aa<- onewayanova:::oneway.factor(z,y)
#' onewayanova:::print.oneway(aa)
#' @export
oneway.factor <- function(z, y, ...) {
## Your code here
b <- data.frame(y,z) # I didn't use cbind as the book said we should avoid that one as much as we can. I also saw that by using cbind the factor are changing to the numbers.
nlist <- with(b, split(y,z))

out <- oneway.default(nlist)
out
}



#' Extracting factor and data from two input columns of data
#' @param formula is regular formula in r which is defining two variables by ~ sign
#' @param data is a list that can have more two columns so its two col will be extracted by formula
#' @seealso \code{\link{print.oneway}}
#' @seealso \code{\link{print.summary.oneway}}
#' @examples
#' zz <- list (factor = c(rep("A",4),rep("B",6), rep("C",6), rep("D",8)) ,
#'   values = c (62, 60 ,63 ,59 ,63 ,67, 71 ,64 ,65 ,66 ,68 ,66, 71 ,67 ,68 ,68, 56,62 ,60 ,61, 63, 64 ,63 ,59 ),
#'   factor2= c(rep("A",4),rep("B",2)) , factor3= c(rep("F",3)) , values2 = c( 63 ,67 ,71, 64, 65, 66),
#'   values2 =c(23,45,67) )
#' onewayformula2 <- onewayanova:::oneway.formula(factor ~ values , zz)
#' onewayanova:::print.oneway(onewayformula2)
#' @export
oneway.formula <- function(formula, data=list(), ...) {

  out <- model.frame(formula , data)
  fac <- out[,1]
  values <- out[,2]
  o2 <- oneway.factor(fac,values)
  return(o2)
}



#' Printing DF and SS
#' @param x x is the input from oneway.default which has the values of anova
#' @param ... ... is the addition argument here
#' @seealso \code{\link{print.oneway}}
#' @seealso \code{\link{print.summary.oneway}}
#' @examples
#' library(faraway)
#' data(coagulation)
#' attach(coagulation)
#' k <- onewayanova:::oneway.default(coagulation)
#' onewayanova:::print.oneway(k)
#' @export
print.oneway <- function(x, ...) {
  tab <- with(x, cbind ( Df = c(x$dfb , x$dfw),
                         SumSq = c(x$ssb , x$ssw) )
              )
  rownames(tab) <- c("diet" , "Residuals")
  printCoefmat(tab )
}


#' Providing summary information for the output of oneway.defual
#' @param object object is the ouput of oneway.default which has the values of anova table
#' @param ... is an addition argument
#' @seealso \code{\link{print.summary.oneway}}
#' @examples
#' library(faraway)
#' data(coagulation)
#' attach(coagulation)
#' k <- onewayanova:::oneway.default(coagulation)
#' kk <- onewayanova:::summary.oneway(k)
#' onewayanova:::print.oneway(k)
#' @export
summary.oneway <- function(object, ...) {
  msb <- object$ssb / object$dfb
  msw <- object$ssw / object$dfw
  F <- msb / msw
  p <- pf( F , object$dfb , object$dfw , lower.tail = FALSE)

  tab <- with(object, cbind ( Df = c(object$dfb , object$dfw),
                              SumSq = c(object$ssb , object$ssw) ))
  rownames(tab) <- c("diet" , "Residuals")
  return(tab)

}


#' print.summary.oneway gets the input from the summary.oneway to print it on the screen
#' @param x x  is the ouput of summary.oneway which has the values of anova table
#' @param ... ... is an addition argument
#' @seealso \code{\link{print.oneway}}
#' @examples
#' library(faraway)
#' data(coagulation)
#' attach(coagulation)
#' oo <- onewayanova:::oneway.default(coagulation)
#' pp <- onewayanova:::summary.oneway(oo)
#' onewayanova:::print.summary.oneway(pp)
#' @export
print.summary.oneway <- function(x, ...) {

 dfb <-  x[1]
 dfw <- x[2,1]
 ssb <- x[1,2]
 ssw <- x[2,2]

  # oneway.table  <- function(x){
  msb <- ssb / dfb
  msw <- ssw / dfw
  F <- msb / msw
  p <- pf( F , dfb , dfw , lower.tail = FALSE)

  tab <-  cbind ( Df = c(dfb , dfw),
                         SumSq = c(ssb , ssw),
                         MeanSq = c(msb,msw),
                         Fvalue= c(F,NA),
                         "Pr(>F)" = c(p,NA))
  rownames(tab) <- c("diet" , "Residuals")
  printCoefmat(tab, P.values = TRUE, has.Pvalue=TRUE, signif.stars = TRUE, na.print="")
}


#' lsmeans.oneway gets the input from the oneway.defual to calculate the LSD for it
#' @param object object is the ouput of oneway.default which has the values of anova table
#' @param ... gets an additional argument
#' @export
lsmeans.oneway <- function(object, ...) {


   LSD <- numeric(); sig_dif <- numeric();  LSD <- matrix(LSD,nrow= (object$dfb+1) ,ncol= (object$dfb+1) ) ;    sig_dif <- matrix(sig_dif,nrow= (object$dfb+1) ,ncol= (object$dfb+1) ) ;    a <- numeric()

   for ( i in seq(1,(object$dfb)+1))
     for ( j in seq(i+1,(object$dfb)+1)) {
      if (i!=object$dfb+1) {
      LSD[i,j] <- ( qt(0.975,object$dfw) * sqrt( object$msw * ( ( 1/ object$grplength[i]) + (1/ object$grplength[j]) )) )

      a <- (abs(object$grpMean[i]-object$grpMean[j]))

      if (a < LSD[i,j])
      {sig_dif[i,j] <- T ;
       cat('LSD between group =\"',i,'\" and group =\"',j , '\"' , 'is' , '\"', LSD[i,j], '\"') ;
       cat(';Therefore No Significant difference is found on Mean difference between Group' ,i, '& Group', j , 'by knowing that  the difference between their mean is', a, '\n\n');
      }
      else
      {sig_dif[i,j] <- F ;
      cat('LSD between group =\"',i,'\" and group =\"',j , '\"' , 'is' , '\"', LSD[i,j], '\"') ;
      cat('; Therefore Significant difference is found on Mean difference between Group' ,i, '& Group', j , 'by knowing that  the difference between their mean is', a, '\n\n');
      }
     }

     }
}


#' plot.oneway gets input from the oneway.defualt to plot the boxes side by side to give the information about each group
#' @param x x is the ouput of oneway.default which has the values of input list
#' @param ... ... is and additional argument
#' @examples
#' library(faraway)
#' data(coagulation)
#' attach(coagulation)
#' ll <- onewayanova:::oneway.default(coagulation)
#' onewayanova:::plot.oneway(ll)
#' @export
plot.oneway <- function(x, ...) {
## Your code here
  boxplot(x$data ,col="lawngreen")
}
```


