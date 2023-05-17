import sys
input_path = sys.argv[1]
output_path = sys.argv[2]

from pdf2image import convert_from_path
images = convert_from_path(input_path, 300)
for image in images:
    image.save(output_path)