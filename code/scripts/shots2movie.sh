#!/bin/sh

PREFIX="${1}"
OUTPUT=${2:-1.mp4}
FPS_IN=30
FPS_OUT=${3:-30}

# quanta-25316-2.bmp BMP3 800x600 800x600+0+0 8-bit sRGB 1.44MB 0.000u 0:00.000
SIZE=`identify "${PREFIX}"-1.bmp | cut -d " " -f 3`

echo "Creating an $SIZE movie..."

avconv -r "${FPS_IN}" -f image2 -s "${SIZE}" -i "${PREFIX}"-%d.bmp -vcodec libx264 -crf 25 -r "${FPS_OUT}" "${OUTPUT}"

echo "All done. ${OUTPUT} was created."
