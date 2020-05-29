set pm3d map

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

set terminal qt 3000
set title "Epoch 3000"
splot "epoch_3000_exit.txt" using 1:2:4
