import os

################################
tblis = True

CC = 'gcc-12'
includes = ['-I/Users/murali/softwares/core/include', '-I/opt/homebrew/Cellar/openblas/0.3.21/include'] ## requires netCDF and cblas headers
libs = ['']
CFLAGS = ''

## For gcc you can use below tags
#CFLAGS = '-O2 -std=c99 -Werror -Wall -Wextra -Wpedantic -Wformat=2 -Wformat-overflow=2 -Wformat-truncation=2 -Wformat-security -Wnull-dereference -Wstack-protector -Wtrampolines -Walloca -Wvla -Warray-bounds=2 -Wimplicit-fallthrough=3 -Wconversion -Wshift-overflow=2 -Wcast-qual -Wstringop-overflow=4  -Warith-conversion -Wlogical-op -Wduplicated-cond -Wduplicated-branches -Wformat-signedness -Wshadow -Wstrict-overflow=4 -Wundef -Wstrict-prototypes -Wswitch-default -Wswitch-enum -Wstack-usage=1000000 -Wcast-align=strict -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fstack-clash-protection '
CFLAGS = '-O3 -std=c99  -Wall -Werror  -Wextra -Wpedantic -Wformat=2 -march=native -fopenmp -fopenmp-simd'


#includes = ['-I/home/users/mnalabothula/softwares/tblis/library/include', '-I/opt/apps/resif/aion/2020b/epyc/software/OpenBLAS/0.3.12-GCC-10.2.0/include', '-I/opt/apps/resif/aion/2020b/epyc/software/netCDF/4.7.4-gompi-2020b/include' ]

#includes = ['-I/home/users/mnalabothula/softwares/iris_nd_array/include' , '-I/opt/apps/resif/iris/2020b/broadwell/software/imkl/2020.4.304-iimpi-2020b/mkl/include', '-I/opt/apps/resif/iris/2020b/broadwell/software/netCDF/4.7.4-iimpi-2020b/include']


#################################
### Compiling starts from here

if tblis:
    CFLAGS +='  -DCOMPILE_ND_TBLIS   '
    CC     +='  -DCOMPILE_ND_TBLIS   '

INCLUDE = ' '
LIBS = ' '

for i in includes:
    INCLUDE = INCLUDE + ' ' + i

for i in libs:
    LIBS = LIBS + ' ' + i

types = ['INT', 'FLOAT', 'DOUBLE', 'SINGLE_COMPLEX', 'DOUBLE_COMPLEX']

compile_c_files = ['alloc.c',  'array.c',  'linalg.c',  'nd_ulit.c',  'netcdf_io.c']

def remove_hash(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    with open(filename, "w") as f:
        for line in lines:
            if not line.strip().startswith("#"):
                f.write(line)


def generate_temp_header():
    with open('nd_array_src.h', "r") as f:
        lines = f.readlines()
    counter =0
    with open("generate_headers.h", "w") as f:
        for line in lines:
            if line.strip().startswith("/* PYTHON_HEADER_FILE START */"):
                counter = 1
            if counter >0:
                f.write(line)
# sed -i '/Baeldung/d' myfile.txt


compile_tag = CC + '  ' +  ' -c -fPIC  ' + '  ' + CFLAGS

asm_tag = CC + '  ' +  ' -S -fPIC  ' + '  ' + CFLAGS

generate_temp_header()

for i in types:
    print(CC + '  ' +  " -E generate_headers.h -DCOMPILE_ND_%s > nd_%s.h " %(i,i))
    os.system(CC + '  ' +  " -E generate_headers.h -DCOMPILE_ND_%s > nd_%s.h " %(i,i))
    remove_hash("nd_%s.h" %(i))

ar_tag = 'ar rcs lib_nd_array.a'
for i in types:
    for j in compile_c_files:
        ar_tag = ar_tag + ' ' + "nd_%s_%s.o" %(j.strip()[:-2],i)
        print(compile_tag + " " + INCLUDE + " " + " %s -DCOMPILE_ND_%s  -o nd_%s_%s.o " %(j,i,j.strip()[:-2],i))
        os.system(compile_tag + " " + INCLUDE + " " + " %s -DCOMPILE_ND_%s  -o nd_%s_%s.o " %(j,i,j.strip()[:-2],i))
        os.system(asm_tag + " " + INCLUDE + " " + " %s -DCOMPILE_ND_%s  -o nd_%s_%s.S " %(j,i,j.strip()[:-2],i))

# compile_tag = compile_tag + " " + "-L/opt/homebrew/Cellar/openblas/0.3.20/lib -lopenblas -L/Users/murali/softwares/core/lib -lnetcdf -lm "
## '-L/opt/homebrew/Cellar/openblas/0.3.20/lib/libopenblas.a -L/Users/murali/softwares/core/lib/libnetcdf.a -lm
os.system("rm generate_headers.h")
print(ar_tag)
os.system(ar_tag)
