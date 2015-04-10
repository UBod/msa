//@@TODO reconcile /muscle with /muscle3.6

#include "muscle.h"
#include <stdio.h>
#ifdef	WIN32
#include <windows.h>	// for SetPriorityClass()
#include <io.h>			// for isatty()
#else
#include <unistd.h>		// for isatty()
#endif

const char *MUSCLE_LONG_VERSION	= "MUSCLE " SHORT_VERSION "."
#include "svnversion.h"
" ";

int g_argc;
char **g_argv;

int main(int argc, char **argv)
	{
#if	WIN32
// Multi-tasking does not work well in CPU-bound
// console apps running under Win32.
// Reducing the process priority allows GUI apps
// to run responsively in parallel.
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	g_argc = argc;
	g_argv = argv;

	SetNewHandler();
	SetStartTime();
	ProcessArgVect(argc - 1, argv + 1);
	SetParams();
	SetLogFile();

	//extern void TestSubFams(const char *);
	//TestSubFams(g_pstrInFileName);
	//return 0;

	if (g_bVersion)
		{
		printf("%s\n", MUSCLE_LONG_VERSION);
			throw EXIT_SUCCESS;
		}

	if (!g_bQuiet)
		Credits();

	if (MissingCommand() && isatty(0))
		{
		Usage();
			throw EXIT_SUCCESS;
		}

	if (g_bCatchExceptions)
		{
		try
			{
			//Run();
			}
		catch (...)
			{
			OnException();
			 	 throw EXIT_Except;
			}
		}
	else {
		//Run();
	}

	 	 throw EXIT_Success;
	}
