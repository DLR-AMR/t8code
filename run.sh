export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/t8code_build/local/lib 

cd ../t8code_build

make -j V=0
make install -j
cd ../t8code

cd example/remove
mpicxx t8_ring.cxx -I$HOME/t8code_build/local/include -I$HOME/t8code_build/local/include/src -L$HOME/t8code_build/local/lib -lt8 -lp4est -lsc -lm -lz
./a.out
rm a.out