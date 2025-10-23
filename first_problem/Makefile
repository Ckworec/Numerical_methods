# Имя исполняемого файла
EXE = a.exe

# Компилятор и флаги
MAKEFLAGS += --no-print-directory --silent
CXX = g++
CXXFLAGS = -O3

# Список исходников
SRC = $(wildcard *.cpp)

# Цель по умолчанию
all: $(EXE)

# Компиляция
$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXE)

# Запуск
run:
	./$(EXE)

# Очистка
clean:
	del $(EXE)

# Полный цикл: форсированная сборка + запуск
full:
	@$(MAKE) clean
	@$(MAKE) all
	@$(MAKE) run
