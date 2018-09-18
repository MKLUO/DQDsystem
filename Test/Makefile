CC			= g++
CFLAGS		= -std=c++11 -MMD -MP

DQDDIR		= ../DQD

LIBS		= 
INCLUDES	= -I$(DQDDIR)
DBG			= -g -pg
SOURCES		= $(DQDDIR)/MathUtilities.cpp $(DQDDIR)/DQD.cpp 
OBJECTS		= $(SOURCES:.cpp=.o)
DEPS		= $(SOURCES:.cpp=.d) 

all: main $(SOURCES)

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) main.cpp $(LIBS) $(INCLUDES) -o main.exe $(DBG)

MathUtilities.o: %.cpp %.h
	$(CC) $(CFLAGS) $< $(LIBS) $(INCLUDES) -c -o $@ $(DBG)

clean:
	rm *.o *.d $(DQDDIR)/*.o $(DQDDIR)/*.d

-include $(DEPS)
