rm RunWrapper

rm -R output
mkdir output

cp ita_input output/ita_input

cd build
make clean
make install
cd ..
