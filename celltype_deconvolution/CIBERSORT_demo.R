

library(CIBERSORT)

data(LM22)
data(mixed_expr)

View(LM22)
View(mixed_expr)

results <- cibersort(sig_matrix=LM22, mixture_file=mixed_expr[rownames(LM22), 
                                                              c(1,1)])


mixed_expr_subset <- mixed_expr[,1]


