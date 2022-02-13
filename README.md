SgInfo - Space Group Info
=========================

(c) 1994-96 Ralf W. Grosse-Kunstleve

Original SgInfo 1.01 from 1996. The only difference is the
new open source license.

See also:
  - http://cci.lbl.gov/sginfo/
  - https://github.com/rwgk/sglite  # No known bugs but also no documentation.
  - http://cctbx.sf.net/  # Look for sgtbx. No known bugs, but large.

# Build instructions

```
cd sginfo_1_01
clang -o sginfo sgclib.c sgfind.c sghkl.c sgio.c sgsi.c sginfo.c -lm
```

This also works with `gcc` instead of `clang` (and probably with any modern
C compiler).

To "install", simply copy the `sginfo` executable to a directory on your
command-line `PATH`.

You can use the following commands to build the sginfo **shared library**
using CMake:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```
