cd $PARFLOW_DIR
if [ -d "$PARFLOW_DIR/build" ]; then 
  rm -rf "$PARFLOW_DIR/build"
fi
ORIGIN_DIR=${pwd}
mkdir ${PARFLOW_DIR}/build ; cd ${PARFLOW_DIR}/build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$PARFLOW_DIR \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif90` \
  -DHYPRE_ROOT=$HYPRE_DIR \
  -DPARFLOW_AMPS_LAYER=mpi1 \
  -DPARFLOW_AMPS_SEQUENTIAL_IO=true \
  -DPARFLOW_ENABLE_TIMING=true \
  -DSILO_ROOT=$SILO_DIR \
  -DPARFLOW_HAVE_CLM=ON \
  -DALQUIMIA_ROOT=$ALQUIMIA_DIR \
  -DCRUNCH_ROOT=$CRUNCHFLOW_DIR \
  -DPFLOTRAN_ROOT=$PFLOTRAN_DIR \
  -DCMAKE_BUILD_TYPE=DEBUG \
  -DCMAKE_C_FLAGS="-W -Wall -Wextra" \
  -DCMAKE_CXX_FLAGS="-W -Wall -Wextra" \
  -DCMAKE_Fortran_FLAGS="-W -Wall -Wextra"
make -j 6 VERBOSE=1
make install
cd $ORIGIN_DIR