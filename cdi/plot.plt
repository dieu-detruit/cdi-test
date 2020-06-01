set pm3d map

set size square

set terminal qt 0
set title "original"
splot "exit.txt" using 1:2:4

set terminal qt 10
set title "Epoch 10"
splot "epoch_10_exit.txt" using 1:2:4

set terminal qt 20
set title "Epoch 20"
splot "epoch_20_exit.txt" using 1:2:4

set terminal qt 30
set title "Epoch 30"
splot "epoch_30_exit.txt" using 1:2:4

set terminal qt 40
set title "Epoch 40"
splot "epoch_40_exit.txt" using 1:2:4

set terminal qt 50
set title "Epoch 50"
splot "epoch_50_exit.txt" using 1:2:4

set terminal qt 100
set title "Epoch 100"
splot "epoch_100_exit.txt" using 1:2:4
