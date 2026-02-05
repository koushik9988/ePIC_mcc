# Compiler and compiler flags
CXX = g++
CXXFLAGS =  -g -Wall -O2 -std=c++17 -I./include -I./linearalgebra -I/usr/include/hdf5/serial/ -I/usr/include/python3.10 

# Directories
SRCDIR = src
OBJDIR = object

# Library directories
LIB_DIRS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
LIBS = -lhdf5 -lhdf5_cpp -lpthread -lpython3.10

#LIB_DIRS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib
#LIBS = -lhdf5 -lhdf5_cpp -lpthread -lglut -lGLU -lGL

# Source files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)

# Object files
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Executable name
EXECUTABLE = ePIC++

# Default rule
all: welcome $(OBJDIR) compiling_field compiling_linalg $(EXECUTABLE)
	@echo "Compilation complete. Run './$(EXECUTABLE) ./inputfiles/inputfilename' to execute."

# Welcome message
welcome:
	@echo "1D Electrostatic  PIC (Particle-in-Cell) Simulation Program Compilation!"

# Create the object directory
$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Compile each source file into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "Compiling $< ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@


# Linking object files to create executable
$(EXECUTABLE): $(OBJECTS)
	@echo "Linking object files to create $(EXECUTABLE) ..."
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB_DIRS) $(LIBS)

# Clean rule to remove object files and executable
clean:
	@rm -f $(OBJDIR)/*.o $(EXECUTABLE)
	@echo "Cleanup complete."

# Phony targets
.PHONY: all welcome compiling_field compiling_linalg clean
