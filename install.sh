
currentDirect="$PWD"
cd $currentDirect/nlopt-2.4.2
./configure --prefix=$currentDirect
make 
make install

cd ..
make
