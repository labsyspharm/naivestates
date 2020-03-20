FROM rocker/tidyverse:3.6.2

RUN R -e "install.packages(c('optparse','mixtools','egg'))"

COPY . /app/
RUN R CMD INSTALL /app/
