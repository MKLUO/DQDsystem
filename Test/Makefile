CC			= g++
CFLAGS		= -std=c++11 -MMD -MP
LIBS		= 
INCLUDES	= -I../DQD
DBG			= -g -pg

DQDDIR		= ../DQD

SOURCES		= main.cpp $(DQDDIR)/*.cpp
OBJECTS		= $(SOURCES:.cpp=.o)
DEPS		= $(SOURCES:.cpp=.d) 

all: main $(SOURCES)

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) main.cpp $(LIBS) $(INCLUDES) -o main.exe $(DBG)

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $< $(LIBS) $(INCLUDES) -c -o $@ $(DBG)

clean:
	rm *.o *.d $(DQDDIR)/*.o $(DQDDIR)/*.d

-include $(DEPS)
