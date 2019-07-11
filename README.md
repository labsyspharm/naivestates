# naivestates - Inference of cell states using Naive Bayes

## Installation

The package can be installed directly from GitHub with the following commands:
``` r
if( !require(devtools) ) install.packages("devtools")
devtools::install_github( "ArtemSokolov/naivestates" )
```

Once installed, the package can be loaded using the standard `library()` interface. In the remainder of this introduction, we will also make use of several `tidyverse` functions.
``` r
library( naivestates )
library( tidyverse )
```

## Loading the data

The package comes with a small example of 10,000 cells, which is distributed as a stand-alone .csv file. Let's load and examine it:

``` r
fnData <- system.file( "examples/example1_data.csv.gz", package="naivestates" )
X <- read_csv( fnData, col_types=cols() )
# # A tibble: 6 x 45
#   CellId Cell_25546ONHoe… Cell_25546ONA488 Cell_25546ONA555 Cell_25546ONA647
#    <dbl>            <dbl>            <dbl>            <dbl>            <dbl>
# 1    173             9.79             8.26             7.31             7.93
# 2    235             9.57             9.34             8.25             8.42
# 3    345             9.26             8.67             7.52             8.22
# 4    350            10.0              8.23             7.46             7.70
# 5    378            10.8              8.34             7.57             7.66
# 6    392             9.84             8.57             7.59             7.93
# # … with 40 more variables
```

The variable `fnData` will contain the filename of where the data lives on disk (for example, `/home/sokolov/R/x86_64-pc-linux-gnu-library/3.6/naivestates/examples/example1_data.csv.gz`). This file can always be examined directly by other tools and programming languages.

After loading the .csv file, `X` now contains a cell-by-feature data frame consisting of Cell IDs in the first column, followed by marker expression values in subsequent columns. The last few columns also contain spatial characteristics, such as the cell position in the image, as well as its area and perimeter.

