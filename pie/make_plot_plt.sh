#!/bin/zsh

echo "set pm3d map" > plot.plt
echo "set size square" >> plot.plt

echo "set terminal qt -1" >> plot.plt
echo "set title 'original object'" >> plot.plt
echo "splot 'data_rpie/object.txt' using 1:2:3\n" >> plot.plt

for i in $(seq 20); do
    echo "set terminal qt ${i}" >> plot.plt
    echo "set title 'epoch ${i}'" >> plot.plt
    echo "splot 'data_rpie/epoch_${i}_object.txt' using 1:2:3\n" >> plot.plt
done
