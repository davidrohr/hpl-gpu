CXXFLAGS     = -fomit-frame-pointer -O3 -march=barcelona -Wall -Wimplicit \
	       -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar \
	       -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type \
	       -Wno-unused-function -pedantic -Wno-long-long -Wshadow \
	       -fno-threadsafe-statics -Wno-variadic-macros -std=c++98 \
	       -I../../../include -L../../../lib -ltbb -Wl,-rpath=../../../lib

HPL_dlaswp01T: HPL_dlaswp01T.cpp ../HPL_dlaswp01T.cpp matrix.h
	@$(CXX) $(CXXFLAGS) -o $@ $<
