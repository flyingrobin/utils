library(doSNOW)
library(tcltk)

cl <- makeSOCKcluster(2)
registerDoSNOW(cl)

pb <- txtProgressBar(max=100, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
r <- foreach(i=1:100, .options.snow=opts) %dopar% {
  Sys.sleep(1)
  sqrt(i)
}
close(pb)
