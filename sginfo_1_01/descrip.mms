cflags = $(cflags)/optimize=inline
rtlib = sys$library:vaxcrtl

sginfo.exe : sginfo.obj sgclib.obj sgio.obj sgfind.obj sghkl.obj sgsi.obj
        $(link)$(linkflags)/notrace sginfo.obj, sgclib.obj, sgio.obj, sgfind.obj, sghkl.obj, sgsi.obj, $(rtlib)/lib

sgtest.exe : sgtest.obj sgclib.obj sgio.obj sgfind.obj
        $(link)$(linkflags)/notrace sgtest.obj, sgclib.obj, sgio.obj, sgfind.obj, $(rtlib)/lib

sgquick.exe : sgquick.obj sgclib.obj sgio.obj
        $(link)$(linkflags)/notrace sgquick.obj, sgclib.obj, sgio.obj, $(rtlib)/lib

.c.obj :
        $(cc) $(cflags) $(mms$source)
