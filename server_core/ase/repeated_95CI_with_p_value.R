# read data, change table name for test
args <- commandArgs(trailingOnly = TRUE)
data <- read.table(args[1], header=FALSE)

# do t-test
result <- t.test(data$V1, conf.level = 0.95, mu = 0.5)

mean_value <- result$estimate
p_value <- result$p.value
lower_ci <- result$conf.int[1]
upper_ci <- result$conf.int[2]

# output result
output_file <- args[2]
write.table(data.frame(mean_value, lower_ci, upper_ci, p_value), file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE)
