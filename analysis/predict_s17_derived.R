# Input data must contain age (years), whtr, and binary male; V5nf also needs
# binary htn and inactive. These codes must follow S17's cycle-specific definitions.
predict_s17_derived <- function(model, data) {
  if(is.character(model)) model <- readRDS(model)
  data$a <- data$age/10-5
  data$w <- 10*(data$whtr-.55)
  x <- as.matrix(data[,model$variables,drop=FALSE])
  p <- plogis(as.numeric(model$beta[1]+x%*%model$beta[-1]))
  data.frame(probability=p,R0=as.integer(p>=model$thresholds["single"]),R1=as.integer(p>=ifelse(data$age>=60,model$thresholds["old"],model$thresholds["single"])))
}
