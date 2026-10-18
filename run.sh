#!/bin/bash
set -e

# Цвета для вывода
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Функция: прогресс-бар
# $1 - текущее значение, $2 - всего, $3 - текст
progress_bar() {
    local current=$1
    local total=$2
    local label=$3
    local width=40
    
    local percent=$(( current * 100 / total ))
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))
    
    printf "\r${YELLOW}%s${NC} [" "$label"
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' ' '
    printf "] %3d%% (%d/%d)" "$percent" "$current" "$total"
    
    if [ "$current" -eq "$total" ]; then
        printf "\n"
    fi
}

# 1. Очистка
echo -e "${GREEN}[1/3]${NC} Очистка старых кадров..."
rm -f frames/frame*.ppm out.mp4
mkdir -p frames

# 2. Симуляция с прогрессом
echo -e "${GREEN}[2/3]${NC} Запуск симуляции..."

# Запускаем симуляцию в фоне, а сами следим за количеством кадров
./sph_sim &
SIM_PID=$!

# Ожидаемое количество кадров (можно вынести в переменную или узнать заранее)
# Если неизвестно — покажем просто счётчик
EXPECTED_FRAMES=${EXPECTED_FRAMES:-0}

FRAME_COUNT=0
while kill -0 "$SIM_PID" 2>/dev/null; do
    FRAME_COUNT=$(ls frames/frame*.ppm 2>/dev/null | wc -l)
    if [ "$EXPECTED_FRAMES" -gt 0 ]; then
        progress_bar "$FRAME_COUNT" "$EXPECTED_FRAMES" "Симуляция"
    else
        printf "\r${YELLOW}Симуляция${NC} — сгенерировано кадров: %d" "$FRAME_COUNT"
    fi
    sleep 0.2
done

wait "$SIM_PID"
TOTAL_FRAMES=$(ls frames/frame*.ppm 2>/dev/null | wc -l)
printf "\r${GREEN}Симуляция завершена${NC} — всего кадров: %d\n" "$TOTAL_FRAMES"

# 3. Сборка видео с парсингом прогресса ffmpeg
echo -e "${GREEN}[3/3]${NC} Сборка видео..."

TOTAL_SECONDS=$TOTAL_FRAMES   # для -framerate 30 число секунд = кадров/30, но считаем точнее:

# Считаем длительность в секундах (целое, с округлением вверх)
DURATION_SEC=$(( (TOTAL_FRAMES + 29) / 30 ))
[ "$DURATION_SEC" -lt 1 ] && DURATION_SEC=1

# Запускаем ffmpeg и читаем его stderr
ffmpeg -y -framerate 30 -i frames/frame%04d.ppm \
       -c:v libx264 -pix_fmt yuv420p out.mp4 \
       -loglevel error -stats 2>&1 | \
while IFS= read -r line; do
    if [[ "$line" =~ time=([0-9]+):([0-9]+):([0-9]+)\.([0-9]+) ]]; then
        h=${BASH_REMATCH[1]}; m=${BASH_REMATCH[2]}; s=${BASH_REMATCH[3]}
        cur=$((10#$h * 3600 + 10#$m * 60 + 10#$s))
        # не даём превысить максимум
        [ "$cur" -gt "$DURATION_SEC" ] && cur=$DURATION_SEC
        progress_bar "$cur" "$DURATION_SEC" "Кодирование"
    fi
done

# ffmpeg мог упасть — проверяем код возврата через PIPESTATUS
# (нужно после пайпа, но т.к. мы в while-подshell — проверяем по факту файла)
if [ -s out.mp4 ]; then
    # Дорисовываем финальный 100%
    progress_bar "$DURATION_SEC" "$DURATION_SEC" "Кодирование"
    echo -e "${GREEN}Готово:${NC} out.mp4"
else
    echo -e "\033[0;31mОшибка:${NC} out.mp4 не создан"
    exit 1
fi