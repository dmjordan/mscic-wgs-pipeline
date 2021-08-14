ml openmpi gcc
#wget https://cran.r-project.org/src/contrib/Rmpi_0.6-9.1.tar.gz
#R CMD INSTALL Rmpi_0.6-9.1.tar.gz --configure-args="   \
#   --with-Rmpi-include=$OPENMPI_HOME/include           \
#   --with-Rmpi-libpath=$OPENMPI_HOME/lib               \
#   --with-Rmpi-type=OPENMPI"
#echo 'install.packages("doMPI")' | R --vanilla
R CMD INSTALL doMPI
