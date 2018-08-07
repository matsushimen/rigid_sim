if find result/position.txt; then
echo "result exist"
icc -parallel -o output output.c
./output result/position.txt
else
mkdir result
icc make_pos.c -lm -o make_pos
./make_pos > result/position.txt
icc -parallel -o output output.c
./output result/position.txt
fi
