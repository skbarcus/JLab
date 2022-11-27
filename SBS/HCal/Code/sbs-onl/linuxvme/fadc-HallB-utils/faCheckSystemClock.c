/*
 * File:
 *    faCheckSystemClock.c
 *
 * Description:
 *    return the status of the system clock
 *
 *
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "jvme.h"
#include "fadcLib.h"

int faTestSystemClock(int id, int pflag);
extern volatile struct fadc_struct *FAp[(FA_MAX_BOARDS+1)]; /* pointers to FADC memory map */
extern int fadcID[FA_MAX_BOARDS];                           /* array of slot numbers for FADCs */
extern int nfadc;


#define FADC_ADDR (3<<19)

int
main(int argc, char *argv[])
{

    int32_t status;
    int i,FA_SLOT;

    printf("\nJLAB fadc Lib Tests\n");
    printf("----------------------------\n");

    status = vmeOpenDefaultWindows();
    if(status != OK)
      goto CLOSE;

    vmeBusLock();
    /* Set the FADC structure pointer */
    faInit((unsigned int)(FADC_ADDR),(1<<19),18, 0x25);

    faGStatus(0);

    faTestSystemClock(0,1);

    vmeBusUnlock();

 CLOSE:
    vmeCloseDefaultWindows();

    exit(0);
}

int
faTestSystemClock(int id, int pflag)
{
  unsigned int rval=OK;

  if(id==0) id=fadcID[0];

  if((id<=0) || (id>21) || (FAp[id] == NULL))
    {
      printf("%s: ERROR : ADC in slot %d is not initialized \n",__FUNCTION__,
	     id);
      return ERROR;
    }

  /* Enable test mode */
  faTestSetSystemTestMode(id, 1);

  /* reset clock counter */
  faTestResetClock250Counter(id);

  int iwait = 0;
  /* Wait for the 20us internal timer */
  while(iwait++ < 50)
    {
      if(faTestGetClock250CounterStatus(id) == 0)
	break;
    }

  /* Counter should return 5000 if the system clock is 250Mhz */
  int expected = 5000, measured = 0, diff = 0;

  measured = faTestGetClock250Counter(id);
  diff = abs(expected - measured);

  if(diff < 5)
    rval = OK;
  else
    rval = ERROR;

  /* Disable test mode */
  faTestSetSystemTestMode(id, 0);

  if(pflag)
    {
      printf("%s: System Clock is %s\n",
	     __func__, (rval==OK) ? "Present" : "NOT PRESENT");
    }

  return rval;
}
