depth = read.csv("OsHK1_mutants_12_2017_clean.csv", header=T)
depth_df=data.frame(depth)

model=aov(depth.mm.~genotype, data=depth_df)

#pairwise ttests
pairwise.t.test(depth_df$depth.mm., depth$genotype, p.adj = "bonf", alternative ="two.sided")

#confidence intervals
confint(model)
