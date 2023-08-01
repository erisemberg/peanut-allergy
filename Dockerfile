FROM rocker/verse:4.2.1

# The noted packages fail to install, even though they appear successful during docker build, but then they install fine using install.packages() from within the docker container 

# Install R packages
RUN R -e "install.packages('readxl')"
RUN R -e "install.packages('stringi')"
RUN R -e "install.packages('lme4')"
RUN R -e "install.packages('MESS')" # failed
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('shades')" # failed
RUN R -e "install.packages('scales')"
RUN R -e "install.packages('extRemes')" # failed

RUN R -e "install.packages('factoextra')"
RUN R -e "install.packages('stringr')"
RUN R -e "install.packages('bestNormalize')" # failed
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('corrplot')"