gcc-11 -O3 -std=c99  -Werror -Wall -Wextra -Wpedantic -Wformat=2 -Wformat-overflow=2 -Wformat-truncation=2 -Wformat-security -Wnull-dereference -Wstack-protector -Wtrampolines -Walloca -Wvla -Warray-bounds=2 -Wimplicit-fallthrough=3 -Wconversion -Wshift-overflow=2 -Wcast-qual -Wstringop-overflow=4  -Warith-conversion -Wlogical-op -Wduplicated-cond -Wduplicated-branches -Wformat-signedness -Wshadow -Wstrict-overflow=4 -Wundef -Wstrict-prototypes -Wswitch-default -Wswitch-enum -Wstack-usage=1000000 -Wcast-align=strict -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fstack-clash-protection -fPIE  test.c -L./../build/ -l_nd_array -L/opt/homebrew/Cellar/openblas/0.3.20/lib -lopenblas -L/Users/murali/softwares/core/lib -lnetcdf -ltblis -lm
