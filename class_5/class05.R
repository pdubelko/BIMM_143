#'---
#'title: "Hands on Assessment: Graphs"
#'author: "Paige Dubelko"
#'date: "January 22, 2019"
#'output: word_document
#'---

#Whisker plot
x <- rnorm(1000,0)
boxplot(x)
# Age vs Weight - offset graph
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(weight, type = "o", pch = 15, cex = 1.5, lwd = 2, ylim = c(2,10),
     xlab = "Age (months)", ylab = "Weight (kg)",
     main = "Baby Weight with Age")
#Bar Graph
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE,
                    sep = "\t")
par(mar = c(2.5,11,2,2))
barplot(mouse$Count, main = "Mouse Genes", horiz = TRUE,
        xlab = "Count", names.arg = mouse$Feature, las = 1,
        col = rainbow(nrow(mouse)))

#Histrogram
x <- c(rnorm(10000), rnorm(10000) + 4)
hist(x, breaks = 10)

#male_female_counts.txt
gender <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE,
                     sep = "\t")
par(mar = c(5.1,5,5,2.1))
barplot(gender$Count, col = c("blue", "pink"), main = "Gender Count",
        names.arg = gender$Sample, ylim = c(0,20), las = 2)

#up_down_expression.txt
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
palette(c("blue", "grey", "magenta"))
par(mar = c(5.1,5,5,2.1))
plot(genes$Condition1, genes$Condition2, col = genes$State, xlab = "Expression Condition 1",
     ylab = "Expression Condition 2", main = "Expression Levels")

#meth
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
index <- meth$expression > 0
color <- densCols(meth$gene.meth[index], meth$expression[index],
                  colramp = colorRampPalette(c("blue", "green", "red", "yellow")))
plot(meth$gene.meth[index], meth$expression[index], xlab = "Gene", ylab = "Expression",
     col = color, pch = 20, main = "Gene Expression")



