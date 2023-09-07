FROM rocker/verse:4.2.1

# The noted packages fail to install, even though they appear successful during docker build, but then they install fine using install.packages() from within the docker container 

# Install R packages
RUN R -e "install.packages('lme4')"
RUN R -e "install.packages('MESS')" 
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('shades')" 
RUN R -e "install.packages('extRemes')" 

RUN R -e "install.packages('factoextra')"
RUN R -e "install.packages('bestNormalize')" 
RUN R -e "install.packages('corrplot')"