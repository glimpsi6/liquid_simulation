#!/bin/bash
set -e   # остановиться при первой ошибке

# 1. Очистить старые кадры
rm -f frames/frame*.ppm out.mp4
mkdir -p frames

# 2. Запустить симуляцию
./sph_sim

# 3. Собрать видео
ffmpeg -framerate 30 -i frames/frame%04d.ppm \
       -c:v libx264 -pix_fmt yuv420p out.mp4

echo "Готово: out.mp4"