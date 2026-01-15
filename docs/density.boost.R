density.boost <- function(x,x.test)
{
   max.components <- 40
   init.var.factor <- 4
   var.factor <- 1.0
   lambda <- 1.0
   min.acceptable.improvement <- 0.001
   max.iter.with.noimprovement <- 5

   d <- dim(x)[2]
   N <- dim(x)[1]
   N.test <- dim(x.test)[1]

   f.x <- f.x.test <- rep(0,N)

   mu.boost <- matrix(0,nrow=max.components,ncol=d)
   sigma.boost <- array(0,dim=c(d,d,max.components))
   alpha.boost <- c(1,rep(0,max.components-1))

   summary.x <- cov.wt(x)
   mu.boost[1,] <- summary.x$center
   sigma.boost[,,1] <- init.var.factor*summary.x$cov
   f.x <- dmvnorm(x,mu.boost[1,],sigma.boost[,,1])
   f.x.test <- dmvnorm(x.test,mu.boost[1,],sigma.boost[,,1])

   i.component <- 1
   c.noimprovement <- 0
   while((i.component < max.components) &&
         (c.noimprovement < max.iter.with.noimprovement))
   {
      mu <- x[sample(1:N,size=1),]
      sigma <- 4*var(x)
      alpha <- 0.2

      i.bag <- sample((1:N),replace=F,size=N/2)
      i.outofbag <- (1:N)[-i.bag]
      J.old <- -Inf
      min.eigen <- min(eigen(sigma)$values)
      if(min.eigen > 10*.Machine$double.eps)
      {
         g.x <- dmvnorm(x[i.bag,],mu,sigma)
         J.new <- mean(log((1-alpha)*f.x[i.bag] + alpha*g.x))
      }
      while(((J.new == -Inf) || (abs(J.old-J.new) > 1E-6)) 
            && (min.eigen > 10*.Machine$double.eps))
      {
         J.old <- J.new

         # E-step
         p <- alpha*g.x/((1-alpha)*f.x[i.bag] + alpha*g.x)
         if(sum(is.na(p)) > 0)
         {
            p[is.na(p)] <- rep(1,sum(is.na(p)))
            cat("NAs ")
         }

         # M-step
         summary.x <- cov.wt(x[i.bag,],wt=p)
         mu <- summary.x$center
         sigma <- summary.x$cov
         alpha <- mean(p)
         min.eigen <- min(eigen(sigma)$values)
         if(min.eigen > 10*.Machine$double.eps)
         {
            # E-step
            g.x <- dmvnorm(x[i.bag,],mu,sigma)
            J.new <- mean(log((1-alpha)*f.x[i.bag] + alpha*g.x))
         }
         print(c(J.old,J.new))
      }

      # update density estimate
      if(min.eigen > 10*.Machine$double.eps)
      {
         alpha <- lambda*alpha

         g.x <- dmvnorm(x,mu,sigma)
         J.new <- mean(log((1-alpha)*f.x[i.outofbag]+alpha*g.x[i.outofbag]))
         J.old <- mean(log(f.x[i.outofbag]))
         delta.J <- J.new-J.old
         if((J.new != -Inf) && (J.old != -Inf) &&
            (delta.J/abs(J.old) < min.acceptable.improvement))
         {
            c.noimprovement <- c.noimprovement + 1
            cat("Reject component\n")
         }
         else
         {  # add new component
            f.x <- (1-alpha)*f.x + alpha*g.x
            g.x <- dmvnorm(x.test,mu,sigma)
            f.x.test <- (1-alpha)*f.x.test + alpha*g.x
            i.component <- i.component + 1
            mu.boost[i.component,] <- mu
            sigma.boost[,,i.component] <- var.factor*sigma
            alpha.boost <- alpha.boost*(1-alpha)
            alpha.boost[i.component] <- alpha
            cat("Adding component",round(mean(log(f.x.test)),3),"\n")
         }
         print(c(J.new,J.old,delta.J/abs(J.old),c.noimprovement))
      }
      else
      {
         cat("Reject component, sigma < 0\n")
         c.noimprovement <- c.noimprovement + 1
      }
   }

   J.boost <- mean(log(f.x.test))
   return(list(f.x=f.x,f.x.test=f.x.test,J.boost=J.boost))
}
