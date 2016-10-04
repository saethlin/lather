rm -r build
rm lather.cpython-35m-x86_64-linux-gnu.so
python setup.py build_ext --inplace
python script.py