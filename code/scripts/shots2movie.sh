#!/bin/sh

PREFIX="${1}"
OUTPUT=${2:-1.mp4}

# quanta-25316-2.bmp BMP3 800x600 800x600+0+0 8-bit sRGB 1.44MB 0.000u 0:00.000
SIZE=`identify "${PREFIX}"-1.bmp | cut -d " " -f 3`

echo "Creating an $SIZE movie..."

avconv -r 60 -f image2 -s "${SIZE}" -i "${PREFIX}"-%d.bmp -vcodec libx264 -crf 25 "${OUTPUT}"

echo "All done. ${OUTPUT} was created."
