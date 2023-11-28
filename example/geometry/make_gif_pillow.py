import os
import sys
from PIL import Image
from natsort import natsorted

exclude = 0
in_filename = []
not_in_filename = []
for arg in sys.argv:
  if ".py" in arg:
    continue
  if not exclude:
    if "!" in arg:
      exclude = 1
    else:
      in_filename.append(arg)
  else:
    not_in_filename.append(arg)
    exclude = 0

files = os.listdir()
files = natsorted(files)
img = []
white_img = []
next = 0
for file in files:
  next = 0
  for check in in_filename:
    if check not in file:
      next = 1
      continue
  if next:
    continue
  if len(not_in_filename):
    for check in not_in_filename:
      if check in file:
        next = 1
        continue
  if next:
    continue
  print(file)
  img.append(Image.open(file))
  white_img.append(Image.new("RGBA", img[-1].size, "WHITE"))
  white_img[-1].paste(img[-1], (0, 0), img[-1])

#(left, upper, right, lower) = (2250, 630, 2250, 1350)
#for i in range(len(img)):
#  img[i] = img[i].crop((left, upper, right, lower))
white_img[0].save("out.gif", save_all=1, append_images=white_img[1:], duration=1000, loop=0, optimize=1)
