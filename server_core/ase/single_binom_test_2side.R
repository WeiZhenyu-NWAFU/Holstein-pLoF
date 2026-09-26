# read data, change table name for test
data <- read.table("uniq.snp.matched_pLoF.stat", sep="\t", header=TRUE)

# do binomial test
result <- lapply(1:nrow(data), function(i) {
  trials <- as.numeric(data[i, "total_count"])
  successes <- as.numeric(data[i, "ref_count"])
  test_result <- binom.test(successes, trials, p=0.5, alternative="two.sided")
  cbind(data[i,], P_Value = test_result$p.value)
})

# output result
outputfile <- "/storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/matched_pLoF/snp/reads_count_filter/uniq.snp.matched_pLoF.stat.binom_test_2side"
write.table(do.call(rbind, result), file=outputfile, sep="\t", row.names=FALSE, col.names=TRUE)

