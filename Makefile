CC = gcc

CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors

TARGET = symnmf

all: $(TARGET)

$(TARGET): symnmf.o
	$(CC) -o $(TARGET) symnmf.o $(CFLAGS) -lm

symnmf.o: symnmf.c symnmf.h
	$(CC) -c symnmf.c $(CFLAGS)

clean:
	@echo "Cleaning up"
	@rm -f *.o symnmf