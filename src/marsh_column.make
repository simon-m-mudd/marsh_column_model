# make with make -f marsh_column.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=marsh_column.cpp \
        ./class_headers/depo_particle.cpp \
        ./class_headers/depo_particle_info.cpp \
        ./class_headers/MuddColTrapping.cpp \
        ./class_headers/sediment_layer.cpp \
        ./class_headers/sediment_stack.cpp \
        ./class_headers/LSDParameterParser.cpp \
        ./class_headers/LSDStatsTools.cpp \
        ./class_headers/marsh_util_fxns.cpp \


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=marsh_column.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

