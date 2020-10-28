FROM rocker/tidyverse:4.0.2

RUN apt-get update && apt-get install libxt6

RUN R -e "install.packages(c('optparse','mixtools','egg','uwot','ggthemes'))"

COPY . /app/
RUN R CMD INSTALL /app/
