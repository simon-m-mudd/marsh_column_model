# make with make -f Column_test.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=Column_test.cpp \
        ./class_headers/depo_particle.cpp \
        ./class_headers/depo_particle_info.cpp \
        ./class_headers/MuddColTrapping.cpp \
        ./class_headers/sediment_layer.cpp \
        ./class_headers/sediment_stack.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Column_test.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
