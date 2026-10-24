CC = gcc
CFLAGS = -Wall -O2 -std=c99
LDFLAGS = -lm

SRCS = main.c \
       model.c \
       soil_layers.c \
       soil_params.c \
       solver.c \
       utils.c \
       plot.c

HEADERS = constants.h \
          model.h \
          soil_params.h \
          soil_layers.h \
          solver.h \
          utils.h \
          plot.h

OBJS = $(SRCS:.c=.o)
TARGET = baresoil_evap.exe
SOIL   ?= data/soil.txt
PARAM  ?= data/param.txt
DATA   ?= data/flux_data.csv

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "链接目标文件..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "编译完成: $(TARGET)"

%.o: %.c $(HEADERS)
	@echo "编译 $<..."
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET) mkdir
	@echo "=========================================="
	@echo "  运行裸土蒸发模型"
	@echo "=========================================="
	@echo "土壤质地文件: $(SOIL)"
	@echo "参数文件    : $(PARAM)"
	@echo "数据文件    : $(DATA)"
	@echo ""
	"$(CURDIR)/$(TARGET)" -soil "$(SOIL)" -param "$(PARAM)" "$(DATA)"
	@echo ""
	@echo "运行完成!"

clean:
	@echo "清理目标文件..."
	rm -f $(OBJS) $(TARGET)
	@echo "清理完成"

clean-all: clean
	@echo "清理输出文件..."
	rm -f output/*.csv output/*.png output/*.gnu output/*.tmp
	@echo "清理完成"

mkdir:
	@mkdir -p output

help:
	@echo "=========================================="
	@echo "  裸土蒸发模型 - Makefile"
	@echo "=========================================="
	@echo ""
	@echo "  目标:"
	@echo "    make                - 编译模型"
	@echo "    make run            - 运行模型（自动创建 output/）"
	@echo "    make clean          - 清理目标文件"
	@echo "    make clean-all      - 清理所有（含输出）"
	@echo "    make mkdir          - 创建输出目录"
	@echo ""
	@echo "  可覆盖的运行参数:"
	@echo "    SOIL=file   土壤质地文件（默认: data/soil.txt）"
	@echo "    PARAM=file  参数文件（默认: data/param.txt）"
	@echo "    DATA=file   观测数据文件（默认: data/flux_data.csv）"
	@echo ""
	@echo "  param.txt 格式（3 列）:"
	@echo "    第1列  地下水位深度 (m, 正向下)"
	@echo "    第2列  通量测量高度 (m)"
	@echo "    第3列  地下水耦合系数 (0~1)"
	@echo ""
	@echo "  示例:"
	@echo "    make run SOIL=data/mysoil.txt PARAM=data/myp.txt DATA=data/myflux.csv"
	@echo "=========================================="

.PHONY: all run clean clean-all mkdir help