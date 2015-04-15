#if	WIN32

typedef unsigned __int64 TICKS;

#pragma warning(disable:4035)
inline TICKS GetClockTicks() {
	unsigned __int64 mmRet;
		/*_asm
			{
			_emit	0x0f
			_emit	0x31
			}*/
		asm("emit, 0x0f");
		asm("emit, 0x31");
		asm("mov dword ptr [mmRet+0],eax");
		asm("mov dword ptr [mmRet+4],edx");
	return mmRet;
}

#define	StartTimer()	__int64 t1__ = GetClockTicks()

#define	GetElapsedTicks()	(GetClockTicks() - t1__)

static double TicksToSecs(TICKS t)
	{
	return (__int64) t/2.5e9;
	}

#endif	// WIN32
