FROM rocker/verse:4.2.1

# Install R packages
RUN R -e "install.packages('readxl')"
RUN R -e "install.packages('stringi')"
RUN R -e "install.packages('factoextra')"
RUN R -e "install.packages('stringr')"
RUN R -e "install.packages('bestNormalize')"
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('corrplot')"