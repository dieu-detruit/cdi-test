set pm3d map

set terminal qt 0
set title "original"
splot "exit.txt" using 1:2:4

set terminal qt 100
set title "Epoch 100"
splot "epoch_99_exit.txt" using 1:2:4

set terminal qt 500
set title "Epoch 500"
splot "epoch_499_exit.txt" using 1:2:4

set terminal qt 750
set title "Epoch 750"
splot "epoch_749_exit.txt" using 1:2:4

set terminal qt 1000
set title "Epoch 1000"
splot "epoch_999_exit.txt" using 1:2:4
