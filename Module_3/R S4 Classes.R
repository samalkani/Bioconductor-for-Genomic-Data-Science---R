# R S4 Classes

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "ALL" and "Genomic Ranges Packages"
BiocManager::install(c("ALL", "GenomicRanges"), force = TRUE)

# 3. Load the Libraries
library(ALL)
library(GenomicRanges)

# 4. LM Object from a linear model
df <- data.frame(y = rnorm(10), x = rnorm(10))
lm.object <- lm(y ~ x, data = df)
lm.object

# 5. Asking for the objects class
class(lm.object)

# 6. Using names() function on the S4 object
names(lm.object)

# 7. Create a list of two elements letters and numbers
xx = list(a = letters[1:3], b = rnorm(4))
xx

# 8. Class XX = lm
class(xx) = "lm"

# 9. Print the lm object
xx

# 10.Load library
library(ALL)

# 11. Load dataset
data(ALL)

# 12. Print dataset
ALL

# 13. Calling Class of the object
class(ALL)

# 14. Helper function isS4() determines whether an object is of the S4 class
isS4(ALL)

# 15. Getting help with an S4 class object
class?ExpressionSet

# 16. Alternative for help with S4 class object
?"ExpressionSet-class"

# 17. Use constructors all the time in R e.g. list() function
xx = list(a = 1:3)

# 18. Constructor for the expression set
ExpressionSet()

# 19. Calling for help with an S4 object gives information 
# on the object and its constructor
?ExpressionSet

# 20. Classic way of defining an Expression set
new("ExpressionSet")

# 21. Defining a S4 class
getClass("ExpressionSet")

# 22. Using the slots with the @ symbol to access data in the expression set
ALL@annotation

# 23. Adjusting screen options
opar = options(width = 80)

# 24. Class of expression set
getClass("ExpressionSet")

# 25. Using the slots with the @ symbol to access data in the expression set
ALL@annotation

# 26. Alternative way of accessing data in a S4 object expression set
slot(ALL, "annotation")

# 27. Alternative way of accessing data in a S4 object expression set
# using an accessor function
annotation(ALL)

# 28. Accessing the help for the class
?ExpressionSet

# 29. Updating old objects with newer versions
# OLD_OBJECT = updateObject(OLD_OBJECT)

# 30. Testing that an object is valid
validObject(ALL)





















