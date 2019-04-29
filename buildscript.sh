cd $PARFLOW_DIR
if [ -d "$PARFLOW_DIR/build" ]; then 
  rm -rf "$PARFLOW_DIR/build"
fi
mkdir build ; cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$PARFLOW_DIR \
  -DHYPRE_ROOT=$HYPRE_DIR \
  -DPARFLOW_AMPS_LAYER=mpi1 \
  -DPARFLOW_AMPS_SEQUENTIAL_IO=true \
  -DPARFLOW_ENABLE_TIMING=true \
  -DSILO_ROOT=$SILO_DIR \
  -DPARFLOW_HAVE_CLM=ON \
  -DALQUIMIA_ROOT=$ALQUIMIA_DIR
make
make install
