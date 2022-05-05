TARGET = doolittle

CC = gcc

CFLAGS = -I. -Werror -g 

DEPS := doolittle.h



%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(TARGET).o
	$(CC) -o $(TARGET) $(TARGET).o  

run: $(TARGET)
	./$(TARGET)

# ignore the error message if clean file is not existed.
clean:
	-@rm *.o $(TARGET) 2>/dev/null || true
