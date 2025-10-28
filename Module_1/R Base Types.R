# R Base Types

# 1. Atomic Vectors (same data type) or Vectors
x = 1:10
x

# 2. Vectors can have names
names(x) = letters[1:10]
class(x)

# 3. Subset Vectors
x[1:3]

# 4. Index based on its names
x[c("a","b")]

# 5. Names do not have to be unique on vectors
x = 1:3
names(x) = c("a","a","b")
x
x["a"]

# 6. Check for duplicates
anyDuplicated(names(x))

# 7. Second element in the names vector is a duplicate
names(x)

# 8. No duplicates
names(x) = letters[1:3]
names(x)
anyDuplicated(names(x))

# 9. The difference between Numerics and Integers
x = 1
class(x)

x = 1:3
class(x)

# 10. Getting a single integer with command line
x = 1L
class(x)

# 11. Machine limit to big integers
.Machine$integer.max
2^31 - 1 == .Machine$integer.max

# 12. A big number BUT slightly smaller than the human genome!
round(.Machine$integer.max / 10^6)

# 13. Convert from integer to numeric data type to avoid the problem of integer limits
as.numeric(1L)

# 14. Matrices - Constructing a matrix
x = matrix(1:9, ncol = 3, nrow = 3)
x

# 15. Rownames
rownames(x) = letters[1:3]
x

# 16. Dimensions of a matrix
dim(x)

# 17. Number of rows and columns in a matrix
nrow(x)
ncol(x)

# 18. First two rows
x[1:2,]

# 19. First two columns
x[,1:2]

# 20. Both at the same time
x[1:2,1:2]

# 21. Sub-setting using characters
x["a",]

# 22. One dimensional Matrix
x["a",, drop=FALSE]

# 23. Sub-setting Matrices
x[x>5]

# 24. Matrices - Constructing a matrix - fills up columns first
x = matrix(1:9, ncol = 3, nrow = 3)
x

# 25. Constructing a matrix filling by row first
x = matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
x

# 26. Constructing a list of different data types
x = list(a = rnorm(3), b=letters[1:5], matrix)
x

# 27. Sub-setting the first two elements in a list
x[1:2]

# 28. Sub-setting the first element
x[1]

# 29. Extracting the first element itself (referencing 2D array)
x[[1]]

# 30. Subsetting the first element by name (classic)
x["a"]

# 31. Subsetting the first element by name (using dollar operator)
x$a

# 32. Partial Matching
names(x) = c("a", "letters", "c")
x$letters
x$let

# 33. Partial Matching can be dangerous
names(x) = c("a", "letters", "letters2")
x$let

# 34. Single bracket notation avoids partial matching and potential errors
names(x) = c("a", "letters", "letters2")
x["letters"]

# 35. Constructing a list where each element in the list is a single number
as.list(1:3)

# 36. Construct a list where each element of the list is of the same type
x = list(rnorm(3), 3:9)
x

# 37. Using lapply() function
lapply(x, mean)

# 38. Generating a vector instead of a list
unlist(lapply(x, mean))

# 39. simplifying lapply() output with sapply() function
sapply(x, mean)

# 40. Constructing a dataframe
x <- data.frame(sex = c("M", "M", "F"), age = c(32, 34, 29))
x

# 41. Selecting columns in a data frame
x$sex

# 42. Using single or double brackets to select a specific column
x[["sex"]]

# 43. Sub-setting the data frame the same way we do for matrices
x[1, "sex"]

# 44. Row names in a data frame have to be unique (1,2,3)
x

# 45. Using sapply and lapply on dataframes
sapply(x, class)

# 46. Converting and data frame into a matrix (character matrix)
as.matrix(x)

# 47. Converting a matrix into a list
as.list(x)

# 48. Conversions in Bioconductor use the as() function from the methods package
library(methods)
as(x, "matrix")









































