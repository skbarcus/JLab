/******************************************************************************
*
*  f1tdcLib.h  - Driver library header file for readout of the JLAB F1 TDC
*                using a VxWorks 5.4 or later based single board
*                computer.
*
*  Author: David Abbott
*          Jefferson Lab Data Acquisition Group
*          November 2003
*
* SVN: $Rev$
*
*/

#ifndef __F1TDCLIB__
#define __F1TDCLIB__

#define F1_BOARD_ID           0xf100
#define F1_MAX_BOARDS             20
#define F1_MAX_TDC_CHANNELS       64
#define F1_MAX_HITS_PER_CHANNEL    8
#define F1_MAX_TDC_CHIPS           8
#define F1_MAX_A32_MEM      0x400000	/* 4 Meg */
#define F1_MAX_A32MB_SIZE  0x1000000	/* 16 MB */
#define F1_VME_INT_LEVEL           2
#define F1_VME_INT_VEC          0xEE

/* Define structure for access to F1 config and control registers */
struct f1tdc_struct
{
  /* 0x000 */ volatile unsigned int csr;
  /* 0x004 */ volatile unsigned int ctrl;
  /* 0x008 */ volatile unsigned int ev_count;
  /* 0x00C */ volatile unsigned int ev_level;
  /* 0x010 */ volatile unsigned int intr;
  /* 0x014 */ volatile unsigned int adr32;
  /* 0x018 */ volatile unsigned int adr_mb;
  /* 0x01C */ volatile unsigned short stat[8];
  /* 0x02C */ volatile unsigned int version;
  /* 0x030          */ unsigned int blank[(0x800 - 0x30)>>2];
  /* 0x800 */ volatile unsigned int config[F1_MAX_TDC_CHIPS][16];
};

struct f1SDC_struct
{
  /* 0x00 */ volatile unsigned short csr;
  /* 0x02 */ volatile unsigned short ctrl;
};

/* Define Signal Distribution Card Bits */
#define F1_SDC_INTERNAL_CLK       0x1
#define F1_SDC_INTERNAL_SYNC      0x2
#define F1_SDC_BUSY_MASK         0x80

#define F1_SDC_RESET             0x1
#define F1_SDC_SYNC_RESET        0x2


/* Define CSR bits */
#define F1_CSR_ENABLE_INPUTS          0x1
#define F1_CSR_CONFIG_ERROR           0x2
#define F1_CSR_CONFIG_ACTIVE          0x4
#define F1_CSR_EVENT_LEVEL_FLAG       0x8
#define F1_CSR_ZERO_EVENTS_FLAG      0x10
#define F1_CSR_BERR_STATUS           0x20
#define F1_CSR_TOKEN_STATUS          0x40
#define F1_CSR_ERROR                 0x80
#define F1_CSR_ERROR_TDC0           0x100
#define F1_CSR_ERROR_TDC1           0x200
#define F1_CSR_ERROR_TDC2           0x400
#define F1_CSR_ERROR_TDC3           0x800
#define F1_CSR_ERROR_TDC4          0x1000
#define F1_CSR_ERROR_TDC5          0x2000
#define F1_CSR_ERROR_TDC6          0x4000
#define F1_CSR_ERROR_TDC7          0x8000
#define F1_CSR_NEXT_BUF_NO        0x10000
#define F1_CSR_BUF0_EMPTY         0x20000
#define F1_CSR_BUF1_EMPTY         0x40000
#define F1_CSR_EOB_FLAG           0x80000
#define F1_CSR_FILLER_FLAG       0x100000
#define F1_CSR_EXTFIFO_EMPTY     0x200000
#define F1_CSR_PULSE_SPARE_OUT 0x10000000
#define F1_CSR_TRIGGER         0x20000000
#define F1_CSR_SOFT_RESET      0x40000000
#define F1_CSR_HARD_RESET      0x80000000


/* Define Control Register bits */
#define F1_FB_SEL                  0x1
#define F1_REF_CLK_PCB             0x2
#define F1_REF_CLK_SEL             0x4
#define F1_EVENT_LEVEL_INT         0x8
#define F1_ERROR_INT              0x10
#define F1_ENABLE_SOFT_TRIG       0x20
#define F1_ENABLE_BERR            0x40
#define F1_ENABLE_MULTIBLOCK      0x80
#define F1_FIRST_BOARD           0x100
#define F1_LAST_BOARD            0x200
#define F1_ENABLE_HEADERS        0x400
#define F1_ENABLE_FP_BUSY        0x800
#define F1_ENABLE_SPARE_OUT     0x1000
#define F1_ENABLE_DATA_TDC0    0x10000
#define F1_ENABLE_DATA_TDC1    0x20000
#define F1_ENABLE_DATA_TDC2    0x40000
#define F1_ENABLE_DATA_TDC3    0x80000
#define F1_ENABLE_DATA_TDC4   0x100000
#define F1_ENABLE_DATA_TDC5   0x200000
#define F1_ENABLE_DATA_TDC6   0x400000
#define F1_ENABLE_DATA_TDC7   0x800000
#define F1_ENABLE_DATA_TDCMASK 0x00FF0000

/* Define Chip register Bits */
#define F1_CHIP_RES_LOCKED             0x1
#define F1_CHIP_HITFIFO_OVERFLOW       0x2
#define F1_CHIP_TRIGFIFO_OVERFLOW      0x4
#define F1_CHIP_OUTFIFO_OVERFLOW       0x8
#define F1_CHIP_EXTFIFO_FULL          0x10
#define F1_CHIP_EXTFIFO_ALMOST_FULL   0x20
#define F1_CHIP_EXTFIFO_EMPTY         0x40
#define F1_CHIP_INITIALIZED           0x80

#define F1_CHIP_CLEAR_STATUS        0x8000

#define F1_CHIP_ERROR_COND          0x1f1e


/* Define CONFIG REGISTER Bits */
#define F1_HIREZ_MODE           0x8000

/* Define Bit Masks */
#define F1_VERSION_MASK        0x00ff
#define F1_BOARD_MASK          0xff00
#define F1_CSR_MASK            0xf0ffffff
#define F1_CSR_ERROR_MASK      0x0000ff00
#define F1_CONTROL_MASK        0x00ff1fff
#define F1_CONTROL_SEL_MASK    0x7
#define F1_EVENT_COUNT_MASK    0xffff
#define F1_EVENT_LEVEL_MASK    0xffff
#define F1_INT_ENABLE_MASK     0x18
#define F1_INT_VEC_MASK        0xff
#define F1_INT_LEVEL_MASK      0x700
#define F1_SLOT_ID_MASK        0x1f0000
#define F1_MODULE_EMPTY_MASK   (F1_CSR_ZERO_EVENTS_FLAG | F1_CSR_BUF0_EMPTY | F1_CSR_BUF1_EMPTY | F1_CSR_EXTFIFO_EMPTY)

#define F1_SDC_MODE_MASK    0x10000
#define F1_SDC_MASK          0xfff0
#define F1_BDC_MASK          0x0008

#define F1_ALL_CHIPS          0xff
#define F1_OFFSET_MASK        0x3f3f
#define F1_ENABLE_EDGES       0x4040	/* Rising edges only - default */
#define F1_ENABLE_DUAL_EDGES  0xc0c0	/* Both rising and falling edges */
#define F1_DISABLE_EDGES_ODD  0xff3f
#define F1_DISABLE_EDGES_EVEN 0x3fff
#define F1_DISABLE_EDGES      0x3f3f

#define F1_A32_ENABLE        0x1
#define F1_AMB_ENABLE        0x1
#define F1_A32_ADDR_MASK     0xffc0
#define F1_AMB_MIN_MASK      0xffc0
#define F1_AMB_MAX_MASK      0xffc00000

#define F1_DATA_TDC_MASK     0x00ffffff
#define F1_DATA_FLAG_MASK    0x07000000
#define F1_DATA_SLOT_MASK    0xf8000000
#define F1_DATA_INVALID      0xf0000000

#define F1_DATA_TYPE_MASK    0x00800000
#define F1_DATA_CHIP_MASK    0x00380000
#define F1_DATA_CHAN_MASK    0x00070000
#define F1_DATA_TIME_MASK    0x0000ffff

#define F1_HEAD_DATA         0x00000000
#define F1_TAIL_DATA         0x00000007
#define F1_DUMMY_DATA        0xf1daffff	/* Slot 30  Data word, Res Not locked, TDC data 0xffff */

#define F1_HT_DATA_MASK      0x00800007
#define F1_HT_CHAN_MASK      0x00000007
#define F1_HT_CHIP_MASK      0x00000038
#define F1_HT_XOR_MASK       0x00000040
#define F1_HT_TRIG_MASK      0x0000ff80
#define F1_HT_EVENT_MASK     0x003f0000
#define F1_HT_TRIG_OVF_MASK  0x00400000


/* Define some macros */


/* Define Prototypes */
int f1Init(unsigned int addr, unsigned int addr_inc, int ntdc, int iFlag);
int f1CheckAddresses();
int f1Slot(unsigned int i);

int f1ConfigWrite(int id, int *config_data, int chipMask);
int f1SetConfig(int id, int iflag, int chipMask);
int f1ConfigRead(int id, unsigned int *config_data, int chipID);
int f1ConfigReadFile(char *filename);
void f1ConfigShow(int id, int chipMask);

void f1Status(int id, int sflag);
void f1GStatus(int sFlag);
void f1ChipStatus(int id, int pflag);
void f1ClearStatus(int id, unsigned int chipMask);

int f1ReadEvent(int id, volatile unsigned int *data, int nwrds, int rflag);
int f1PrintEvent(int id, int rflag);
int f1FlushEvent(int id);
int f1GPrintEvent(int rflag);

void f1Clear(int id);
void f1GClear();

void f1ClearStatus(int id, unsigned int chipMask);
void f1GClearStatus(unsigned int chipMask);

unsigned int f1ErrorStatus(int id, int sflag);
unsigned int f1GErrorStatus(int sflag);

int f1CheckLock(int id);
int f1GCheckLock(int pflag);

void f1Reset(int id, int iFlag);
int f1HardReset(int id);
void f1Trig(int id);
void f1GTrig();

int f1Dready(int id);
int f1Bready(int id);
unsigned int f1GBready();
int f1DataScan(int pflag);
unsigned int f1ScanMask();

int f1GetRez(int id);
int f1SetWindow(int id, int window, int latency, int chipMask);
void f1GSetWindow(int window, int latency, int chipMask);

unsigned int f1ReadCSR(int id);

int f1WriteControl(int id, unsigned int val);
void f1GWriteControl(unsigned int val);

int f1Enabled(int id);
int f1EnableData(int id, int chipMask);
void f1GEnableData(int chipMask);

int f1DisableData(int id);
int f1DisableChannel(int id, int input);

void f1EnableClk(int id, int cflag);
void f1DisableClk(int id, int cflag);

unsigned int f1EnableLetra(int id, int chipMask);
unsigned int f1DisableLetra(int id, int chipMask);

void f1EnableSoftTrig(int id);
void f1GEnableSoftTrig();
void f1DisableSoftTrig(int id);

void f1EnableBusError(int id);
void f1GEnableBusError();
void f1DisableBusError(int id);

int f1SetBlockLevel(int id, int level);
void f1GSetBlockLevel(int level);

void f1SetInputPort(int id, int fb);
void f1EnableMultiBlock();
void f1DisableMultiBlock();

unsigned int f1GetA32(int id);
unsigned int f1GetA32M();
void f1SetClockPeriod(int period);

STATUS f1SDC_Config(int mode);
void f1SDC_Sync();
STATUS f1BDC_Config();
void f1BDC_Sync();

#endif
