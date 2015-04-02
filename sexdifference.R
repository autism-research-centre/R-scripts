
# This script is to calculate sex-difference p-values between males and females as mentioned in Randall et al., 2013, PLOS Genetics. 
# Read the article here: http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003500
# The basic formula is t = bmen-bwomen/sqrt(SEmen^2 + SEwomen^2 -2r*SEmen*SEwomen). 
# This script has been written for 23andme data. Please modify the column names accordingly.


#Step 1: File management
traitM = read.table("filename.txt", header = T) #read input file
traitF = read.table("filename.txt", header = T) #read input file
traitmerge = merge(traitF, traitM, by = "all.data.id") #merge two files by a column, in this case, all.data.id

#Step 2: Calculate Spearman’s rank correlation coefficient
traitmerge2 = traitmerge[,c("all.data.id", "SNP.x", "P.x", "P.y")] #prune the file to make it manageable
dim(traitmerge2) #dimensions of the file, returned as no of rows, no of columns
traitmerge3 = traitmerge2[!duplicated(traitmerge2[,c("all.data.id")]),] #remove duplicates by values in a column, in this case, all.data.id
dim(traitmerge3) #compare the dimension with the earlier file
correlation2 = cor(rank(traitmerge3[,"P.x"]),rank(traitmerge3[,"P.y"])) #calculate Spearman’s rank correlation coefficient
correlation2 #print rank correlation coefficient

#Step 3: Calculate t-statistic
traitmerge = traitmerge[!duplicated(traitmerge[,c("all.data.id")]),] #remove duplicates
traitmerge$numerator = traitmerge[,"effect.x"]-traitmerge[,"effect.y"] #create a new column (numerator) which is the difference between the effects. This is the numerator in the equation.
traitmerge$stderr.x2 = traitmerge[,"stderr.x"]^2 #create a new column with stderr.x squared
traitmerge$stderr.y2 = traitmerge[,"stderr.y"]^2 #create a new column with stderr.y squared
traitmerge$thirdvar = (-correlation2)*traitmerge[,"stderr.x"]*traitmerge[,"stderr.y"] #create a new column that multiplies standarderrors with the negative correlation coefficient. 
traitmerge$denosquare = traitmerge[,"stderr.x2"]+traitmerge[,"stderr.y2"]+traitmerge[,"thirdvar"] #create a new column that adds the products of the above three steps. 
traitmerge$denominator = sqrt(traitmerge[,"denosquare"]) #create a new column that square roots the product of the above step. This is the denominator in the equation
traitmerge$sexdifft = traitmerge[,"numerator"]/traitmerge[,"denominator"] #create a new column with the t-statistics.

#Step 4: Calculate the p-value of the t-statistic
traitmerge$pvalt = 2*pt(abs(traitmerge[,"sexdifft"]), degreesoffreedom, lower = FALSE) #calculate the p-value (two tailed). The p-values are stored in a new column titled pvalt. Insert the degrees of freedom. 

#Step 5: Final data management and writing the file.
traitmergefinal = traitmerge[,c("SNP.x", "freq.a.x","freq.b.x", "CHR.x", "BP.x", "effect_allele.x", "other_allele.x","sexdifft","pval","effect.y", "stderr.y","effect.x", "stderr.x")] #prune the data table. 
traitmergefinal = traitmergefinal[order(traitmergefinal$pvalt),]#order according to p-values in ascending order.
write.table(traitmergefinal, file = "traitsexdifference.txt", row.names = F, col.names = T, quote = F) #write the file as a text file. 
