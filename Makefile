PROG_ROOT = $(HOME)/Desktop/phys404_hw6
PROGRAM_NAME = Potential_1D

CC := g++

INCLUDES := \
	-I$(PROG_ROOT)/include
CXXFLAGS := 
CPPFLAGS :=
CFLAGS :=
LDFLAGS := -std=c++11
LIBS := 

TARGET = bin/$(PROGRAM_NAME)
SOURCE = src/$(PROGRAM_NAME).cpp
OBJECT = $(SOURCE: .cpp=.o)

all : $(TARGET)

%.o : %.cpp
	$(CC) $(INCLUDES) $(CFLAGS) -c -o $@ $<

$(TARGET) : $(OBJECT)
	$(CC) -o $@ $+ $(INCLUDES) $(LDFLAGS) $(LIBS)

clean : $(OBJECT)
	rm -f *.o $(TARGET)
	@echo clean complete...
