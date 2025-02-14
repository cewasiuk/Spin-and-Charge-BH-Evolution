TARGET    = formatting.x
OBJECT_01 = formatting.o
SOURCE_01 = formatting.c

$(TARGET): $(OBJECT_01)
	gcc -g -o $(TARGET) $(OBJECT_01) -lm -std=c99

$(OBJECT_01) : $(SOURCE_01)
	gcc -g -c $(SOURCE_01) -o $(OBJECT_01)

all : $(TARGET)

clean :
	-rm $(TARGET) $(OBJECT_01)
