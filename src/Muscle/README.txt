MUSCLE v3.0 source code README
------------------------------

http://www.drive5.com/muscle

This version of MUSCLE was built and tested on two platforms:
Windows XP and Red Hat Linux 8.0.

On Windows, I used Microsoft Visual C++ .Net, which I find
to be the best C++ compile / edit / test environment I've
tried on any platform. The Microsoft project file is
muscle.vcproj.

The Linux make file is Makefile. This is a very simple-minded
make file (because I am a Linux development novice), so should
be easy to understand. By default, it uses shared libraries,
but I found this to give problems when copying between
different Linux versions. The fix was to use the linker
flag -lm static (commented out), which gives a much bigger
but more portable binary. The posted binary was linked with
static libraries.

The source code was not written to be maintained by anyone
but me, so the usual apologies and caveats apply.

Bob Edgar,
January 2004
