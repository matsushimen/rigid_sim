if find resultfree/position.txt; then
echo "resultfree exist"
icc -parallel -static-intel -o outputfree outputfree.c
./outputfree resultfree/position.txt
else
mkdir resultfree
icc make_pos.c -lm -o make_pos
./make_pos > resultfree/position.txt
icc -parallel -o outputfree outputfree.c
./outputfree resultfree/position.txt
fi
