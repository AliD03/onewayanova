#' oneway for implementing One-Way S3 update
#' @param z z is the input in form of list/data.frame
#' @param ... ... is the addistion arguments here
#' @export
oneway <- function(z, ...) UseMethod("oneway")
  #' oneway.defualt for calculating the anova parameters
  #' @param z z is the input in form of list/data.frame
  #' @param ... ... is the addistion arguments here
  #' @examples
  #' onewayanova:::oneway.default(coagulation)
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
#' @param z z as a factor input here
#' @param y y as values input here
#' @examples
#' library(faraway)
#' data(coagulation)
#' attach(coagulation)
#' z <- coagulation$diet
#' y <- coagulation$coa
#' onewayanova:::oneway.factor(z,y)
#' @export
oneway.factor <- function(z, y, ...) {
## Your code here
b <- data.frame(y,z) # I didn't use cbind as the book said we should avoid that one as much as we can. I also saw that by using cbind the factor are changing to the numbers.
nlist <- with(b, split(y,z))

out <- oneway.default(nlist)
out
}



#' oneway.formula gets a data list and extracts two specific colmumns from the formula in order to make it ready to find the anova variables by passing it to oneway.default
#' @param formula formula is regular formula in r which is defining two variables by ~ sign
#' @param data data is a list that can have more two columns so its two col will be extracted by formula
#' @export
oneway.formula <- function(formula, data=list(), ...) {

  out <- model.frame(formula , data)
  fac <- out[,1]
  values <- out[,2]
  o2 <- oneway.factor(fac,values)
  return(o2)
}



#' print.oneway gets input from the oneway.defual and print the very basic information about anova
#' @param x x is the input from oneway.default which has the values of anova
#' @param ... ... is the addition argument here
#' @examples
#' library(faraway)
#' data(coagulation)
#'
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


#' Summary.oneway gets input from the oneway.defual to for a summary out of it by making a table
#' @param object object is the ouput of oneway.default which has the values of anova table
#' @param ... is an addition argument
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



#### 9. Your S3 class implementation should be illustrated with the *coagulation* data set. The data consists of blood coagulation times for 24 animals randomly assigned to four different diets.

library(faraway)
data(coagulation)

# **VERY IMPORTAN COMMENT** : VERY IMPORTAN COMMENT: VERY IMPORTAN COMMENT. I do not need to change the coagualtion from data.frame to list to get it as an input to my codes since I wrote the following code in the one.way.default to do it for me automatically. So in the following code if the input is data.frame instead of list it will change the input to the list and make our job easy and our code more generic.
#coagulation <- with(coagulation, split(coag,diet)).

onewaydefault <- oneway.default(coagulation)
print.oneway(onewaydefault)

onewayfactor <- oneway.factor(coagulation$diet,coagulation$coag)
print.oneway(onewayfactor)

onewayformula <- oneway.formula(diet ~ coag , coagulation)
print.oneway(onewayformula)

# By uncommenting this line this is also possible to check the oneway fromula with the following list too that has 4 columns. Two of its lists are the same as coag and diet of coagualtion which can be shown that if I choose those two columns still I get the same answer.

#z <- list (factor = c(rep("A",4),rep("B",6), rep("C",6), rep("D",8)) , values = c (62, 60 ,63 ,59 ,63 ,67, 71 ,64 ,65 ,66 ,68 ,66, 71 ,67 ,68 ,68, 56,62 ,60 ,61, 63, 64 ,63 ,59 ),     factor2= c(rep("A",4),rep("B",2)) , factor3= c(rep("F",3)) , values2 = c( 63 ,67 ,71, 64, 65, 66), values2 =c(23,45,67) )
#onewayformula <- oneway.formula(factor ~ values , z) ; print(onewayformula)


summaryoneway <- summary.oneway(onewaydefault)

print.summary.oneway(summaryoneway)

#For checkig to show that the answer is the same as the anova table I use the built-in function here
fit <- aov(coag ~ diet, data=coagulation); summary(fit)

lsmeans.oneway(onewaydefault)
plot.oneway(onewaydefault)




