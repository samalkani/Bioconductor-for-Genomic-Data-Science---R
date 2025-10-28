# IRanges - Basic Usage

# 1. Installing the IRanges package

# BiocManager::install("IRanges")

# 2. load IRanges Package
library(IRanges)

# 3. Specify IRanges 
#    (giving 2 out of 3 arguments, start and end with the 
#    last one width being inferred from start and end values)
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir1

# 4. Specify IRanges by giving start and width arguments
ir2 <- IRanges(start = c(1,3,5), width = 3)
ir2

# 5. Accessing different components of IRanges
start(ir1)

# 6. Set elements of IRanges
width(ir2) <- 1
ir2

# 7. IRanges can have names like any other vector
names(ir1) <- paste("A", 1:3, sep = "")
ir1

# 8. IRanges don't have a dimension they have a length
length(ir1)

# 9. Sub-setting IRanges is like any other vector

  # Integer indexing
  ir1[1]

  # Name indexing
  ir1["A1"]

# 10. Concatenating IRanges
c(ir1, ir2)


# 11. User-defined function to plot Iranges
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) {
  height <- 1
  if (is(xlim, "Ranges")) {
    xlim <- c(min(start(xlim)), max(end(xlim)))
  }
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

  # Display 2 plots one under the other
  par(mfrow = c(2,1))

# 12. Example usage
  # Create an IRanges object
  ir<-IRanges(start = c(1,3,7,9),end = c(4,4,8,10))

# 13. Plot the IRanges object
plotRanges(ir)

# 14. Normal IRange created by the reduce function
reduce(ir)
plotRanges(reduce(ir))

# 15. User-defined function to plot IRanges (Alternative)
plotIRanges <- function(iranges_obj, title = "IRanges Plot") {
  # Validate input
  if (!inherits(iranges_obj, "IRanges")) {
    stop("Input must be an IRanges object.")
  }
  
  # Extract start, end, and width
  starts <- start(iranges_obj)
  ends <- end(iranges_obj)
  widths <- width(iranges_obj)
  
  # Create a plot
  plot(
    NULL, xlim = c(min(starts), max(ends)), ylim = c(0, length(iranges_obj) + 1),
    xlab = "Position", ylab = "Range Index", main = title, type = "n"
  )
  
  # Add horizontal bars for each range
  for (i in seq_along(iranges_obj)) {
    segments(
      x0 = starts[i], x1 = ends[i], y0 = i, y1 = i,
      col = "blue", lwd = 2
    )
  }
  
  # Add grid for better visualization
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
}

  # Display 2 plots one under the other
  par(mfrow = c(2,1))


# 16. Example usage
  # Create an IRanges object
  example_ranges <- IRanges(start = c(1,3,7,9),end = c(4,4,8,10))
  reduced_ranges <- reduce(example_ranges)

# 17. Plot the IRanges object
plotIRanges(example_ranges, title = "Example IRanges Plot")
plotIRanges(reduced_ranges, title = "Reduced IRanges Plot")

# 18. Plot Disjoin IRange (opposite of reduce function)
plotRanges(ir)
plotRanges(disjoin(ir))

# 19. IRange of length 4
ir

# 20. Manipulating IRanges using the resize function
resize(ir, width = 1, fix = "start")

# 21. More useful is to resize from the center of the intervals
resize(ir, width = 1, fix = "center")

# 22. Create two IRanges and look at them
ir1 <- IRanges(start = c(1, 3, 5), width = 1)
ir2 <- IRanges(start = c(4, 5, 6), width = 1)
ir1
ir2

# 23. Union of ir1 and ir2
union(ir1, ir2)

# 24. Intersection of ir1 and ir2
intersect(ir1, ir2)

# 25. IRanges example displayed below
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3)
ir1
ir2

# 26. FindOverlaps function
ov <- findOverlaps(ir1, ir2)
ov

# 27. Verify findOverlaps result by hand
ir1[1]
ir2[1]

# 28. Query Hits
queryHits(ov)

# 29. Using the unique function
unique(queryHits(ov))

# 30. FindOverlaps function arguments
args(findOverlaps)

# 31. Find the number of overlaps
countOverlaps(ir1, ir2)

# 32. Which of the ranges in ir2 and closest to the ranges in ir1?
nearest(ir1, ir2)










