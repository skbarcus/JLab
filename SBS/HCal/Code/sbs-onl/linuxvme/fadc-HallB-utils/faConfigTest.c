/*
 * File:
 *    faConfigTest.c
 *
 * Description:
 *    Program to test read of config file
 *
 *
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "fadc250Config.h"

extern void fadc250PrintGlobals();
int
main(int argc, char *argv[])
{
  char rol_usrConfig[256] = "";

  /* grab filename from arguments */
  if(argc > 1)
    {
      strcpy(rol_usrConfig,argv[1]);
      printf("rol = %s\n", rol_usrConfig);
    }

  fadc250InitGlobals();
  fadc250Config(rol_usrConfig);

  exit(0);
}
