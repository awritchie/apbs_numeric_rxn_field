NAME := apbs_numeric_rxn_field
EXTRA_SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(EXTRA_SRC))

CXX := llvm-g++ 
CXXFLAGS := -O3 -Wall 
INCLUDE := -I/usr/local/Cellar/eigen/3.2.1/include/eigen3 
LDFLAGS :=
LDLIBS :=

$(NAME): $(OBJS)
	$(CXX) $(LDFLAGS) $(LDLIBS) -o $@ $^
$.o: %.c
	$(CXX) $(CXXFLAGS) -c $<
clean:
	rm -f $(NAME) $(OBJS)
print-%:
	@echo '$* =  $($*)'

.PHONY: clean
