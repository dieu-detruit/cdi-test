set pm3d map

set xrange[-1e-7:1e-7]
set yrange[-1e-7:1e-7]

set terminal qt 0
set title "original"
splot "exit.txt" using 1:2:4

set terminal qt 100
set title "Epoch 100"
splot "epoch_100_exit.txt" using 1:2:4

set terminal qt 500
set title "Epoch 500"
splot "epoch_500_exit.txt" using 1:2:4

set terminal qt 750
set title "Epoch 750"
splot "epoch_750_exit.txt" using 1:2:4

set terminal qt 1000
set title "Epoch 1000"
splot "epoch_1000_exit.txt" using 1:2:4

set terminal qt 2000
set title "Epoch 2000"
splot "epoch_2000_exit.txt" using 1:2:4

set terminal qt 5000
set title "Epoch 5000"
splot "epoch_5000_exit.txt" using 1:2:4

set terminal qt 10000
set title "Epoch 10000"
splot "epoch_10000_exit.txt" using 1:2:4
