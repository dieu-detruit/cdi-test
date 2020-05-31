import numpy as np
from PIL import Image

image = np.array(Image.open('taro_kono.jpg').convert('L'))

image = image.transpose()
image = np.fliplr(image).copy()

print(image.shape)

np.savetxt('transobj.txt', image, fmt='%d')
