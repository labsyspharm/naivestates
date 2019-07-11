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

### Channel -> Cell Type/State mapping

The Naive Bayes framework requires predefined mapping of channels/markers to cell types/states. Currently, the framework assumes a binary association between the two (e.g., Immune cells express CD45, Stroma cells express SMA, etc.) As these relationships become more refined, future development will include support for probabilistic assignment of markers to cell types.

The package includes an example mapping to along with the dataset above. As with the data itself, the mapping can be loaded into R or inspected directly by other tools and programming languages.

``` r
fnMap <- system.file( "examples/example1_chnlmap.csv", package="naivestates" )
M <- read_csv( fnMap, col_types=cols() )
# # A tibble: 19 x 2
#    Channel              Class      
#    <chr>                <chr>      
#  1 Cell_25546ONMITF     Tumor      
#  2 Cell_25546ONS100     Tumor      
#  3 Cell_25546ONSMA      Stroma     
#  4 Cell_25546ONVIMENTIN Stroma     
#  5 Cell_25546ONVIMENTIN Mesenchymal
#  6 Cell_25546ONCD4      Immune     
#  7 Cell_25546ONCD4      T-hlp      
#  8 Cell_25546ONCD4      T-reg      
#  9 Cell_25546ONCD3      Immune     
# 10 Cell_25546ONCD8      Immune     
# 11 Cell_25546ONCD8      T-ctx      
# 12 Cell_25546ONCD45RO   Immune     
# 13 Cell_25546ONFOXP3    Immune     
# 14 Cell_25546ONFOXP3    T-reg      
# 15 Cell_25546ONPD1      Immune     
# 16 Cell_25546ONECAD     Epithelial 
# 17 Cell_25546ONCATENIN  Tumor      
# 18 Cell_25546ONKERATIN  Epithelial 
# 19 Cell_25546ONCD45     Immune     
```

There are several important points to highlight:
1. When composing your own marker-cell type relationships, ensure that the resulting data frame has column names `Channel` and `Class`.
2. Not every channel present in the dataset needs to have a mapping.
3. The same channel may map to multiple classes. For example, VIMENTIN is mapped to Stroma and Mesenchymal classes in the example above.
