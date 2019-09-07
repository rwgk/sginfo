add_definitions(-DVERSION="1.01")

add_library(sginfo
  SHARED

  sgclib.c	
  sghkl.c		
  sgio.c		
  sgsi.c
  sgfind.c	
  sginfo.c	
  
  sginfo.h
  )


