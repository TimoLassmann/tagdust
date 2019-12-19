AC_INIT([tldevel], [1.0])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL

AM_SILENT_RULES([no])



LT_INIT
AC_CONFIG_MACRO_DIRS([m4])

m4_include([m4/ax_check_zlib.m4])

AX_CHECK_ZLIB(
[
CPPFLAGS="$CPPFLAGS -I${ZLIB_HOME}/include"
LDFLAGS="$LDFLAGS -L${ZLIB_HOME}/lib"
LIBS="-lz $LIBS"
AC_DEFINE([HAVE_LIBZ], [1],[Define to 1 if you have `z' library (-lz)])
AM_CONDITIONAL([ENABLE_ZLIB], [test "$zlib_cv_zlib_h" = "yes"])
]
,
AC_MSG_WARN([Zlib not found!!!]))

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])
AC_CHECK_LIB([pthread], [pthread_create])

#AC_CHECK_LIB([hts],[fai_build])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([string.h])
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADERS([limits.h])
AC_CHECK_HEADERS([stddef.h])
AC_CHECK_HEADERS([time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T

AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_C_INLINE
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime memset pow sqrt getpid])



AC_CHECK_FUNCS([strdup])
AC_CHECK_FUNCS([strtol])



# non-standard stuff


AC_CANONICAL_HOST
case $host_os in
darwin* )
DEAD_STRIP="-Wl,-dead_strip"
;;
*)
DEAD_STRIP="-Wl,--gc-sections -Wl,--as-needed"
;;
esac
AC_SUBST(DEAD_STRIP)

AC_ARG_ENABLE(debugging,
[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging="no")

case $enable_debugging in
yes)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging. ]);;
1)  AC_DEFINE(DEBUGLEVEL, 1,[Defines debugging level 1.]);;
2)  AC_DEFINE(DEBUGLEVEL, 2,[Defines debugging level 2.]);;
3)  AC_DEFINE(DEBUGLEVEL, 3,[Defines debugging level 3.]);;
no)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging.]);;
*)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac


if test "$enable_debugging" != "no"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS="-ggdb -pedantic -Wall -std=gnu11 -fdata-sections -ffunction-sections"
else
CFLAGS="-O3 -pedantic -funroll-loops -Wall -std=gnu11 -fdata-sections -ffunction-sections"
DEBUG=0
fi


AX_VALGRIND_DFLT([memcheck], [on])
AX_VALGRIND_CHECK


#AC_CANONICAL_HOST

AC_SUBST([DATE]) DATE="`date`"

AC_SUBST([USER]) USER="`whoami`@`hostname`"

AC_SUBST([UNAME_INFO]) UNAME_INFO=`uname -n`

AC_SUBST([COMPILER]) COMPILER=`$CC --version | head -n 1`



ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"


AC_CONFIG_FILES([
Makefile
dev/Makefile
])

AC_CONFIG_HEADERS([config.h])

AC_OUTPUT

AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([General Information:])
AC_MSG_NOTICE([               Version: $VERSION])
AC_MSG_NOTICE([         Configured on: $DATE])
AC_MSG_NOTICE([         Configured by: $USER])
AC_MSG_NOTICE([           Host system: $host_cpu-$host_vendor-$host_os])
AC_MSG_NOTICE([     Uname information: $UNAME_INFO])
AC_MSG_NOTICE([                CFLAGS: $CFLAGS])
AC_MSG_NOTICE([               CPFLAGS: $CPPFLAGS])
AC_MSG_NOTICE([                  LIBS: $LIBS])
AC_MSG_NOTICE([               LDFLAGS: $LDFLAGS])
AC_MSG_NOTICE([    Installation point: $prefix])
AC_MSG_NOTICE([---------------------------------------------------------])


