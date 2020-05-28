import numpy as np
from PIL import Image

image = np.array(Image.open('toilet.png').convert('L'))

image = image.transpose()
image = np.fliplr(image).copy()

np.savetxt('transobj.txt', image, fmt='%d')
