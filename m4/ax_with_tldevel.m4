# SYNOPSIS
#
#   AX_WITH_TLDEVEL
#
# DESCRIPTION
#
#   This macro checks whether tldevel is installed nearby. This macro is
#   in essence a copy of the ax_with_htslib library. 
#
#   The following output variables are set by this macro:
#
#     TLDEVELDIR              Directory containing HTSlib source tree
#
#   The following shell variables may be defined:
#
#     ax_cv_tldevel        Set to "yes" if HTSlib was found
#     ax_cv_tldevel_which  Set to "source", "install", or "none"
#
# LICENSE
#
#   Copyright (C) 2015,2017 Genome Research Ltd
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.
#serial 1

AC_DEFUN([AX_WITH_TLDEVEL], [
AC_ARG_WITH([tldevel],
[AS_HELP_STRING([--with-tldevel=DIR],
    [use the tldevel source tree or installation in DIR])
dnl Not indented, to avoid extra whitespace outwith AS_HELP_STRING()
AS_HELP_STRING([--with-tldevel=system],
    [use only a system tldevel installation])],
  [with_tldevel=system], [with_tldevel=search])

AC_ARG_ENABLE([configure-tldevel],
  [AS_HELP_STRING([--enable-configure-tldevel],
     [run configure for tldevel as well @<:@default=only_in_subdir@:>@])],
  [], [enable_configure_tldevel=only_in_subdir])

AC_MSG_NOTICE([ search: $with_tldevel])

case $with_tldevel in
yes|search)
  AC_MSG_CHECKING([location of tldevel source tree])
  case $srcdir in
    .) srcp= ;;
    *) srcp=$srcdir/ ;;
  esac
  found=
  for dir in ${srcp}tldevel* -- ${srcp}tldevel
  do
    if test "$dir" = "--"; then
      test -n "$found" && break
    elif test -f "$dir/tldevel.c" && test -f "$dir/tldevel.h"; then
      found="${found}1"
      TLDEVELDIR=$dir
    fi
  done
  if test -z "$found"; then
    AC_MSG_RESULT([tldevel not found])
    ax_cv_tldevel_which=system
  elif test "$found" = 1; then
    AC_MSG_RESULT([$TLDEVELDIR])
    ax_cv_tldevel_which=source
    if test "x$enable_configure_tldevel" = "xonly_in_subdir" ; then
      case $TLDEVELDIR in
        "${srcp}tldevel"*) enable_configure_tldevel=yes ;;
        *) ;;
      esac
    fi
  else
    AC_MSG_RESULT([several directories found])
    AC_MSG_ERROR([use --with-tldevel=DIR to select which tldevel to use])
  fi
  ;;
no) ax_cv_tldevel_which=none ;;
system) ax_cv_tldevel_which=system ;;
*)
    TLDEVELDIR=$with_tldevel
  if test -f "$TLDEVELDIR/tldevel.c" && test -f "$TLDEVELDIR/tldevel.h"; then
    ax_cv_tldevel_which=source
  else
    ax_cv_tldevel_which=install
  fi
  ;;
esac


case $ax_cv_tldevel_which in
source)
  ax_cv_tldevel=yes
#CPPFLAGS="$CPPFLAGS -I${TLDEVELDIR}"
  #LDFLAGS="$LDFLAGS -L${TLDEVELDIR}"
  TLDEVEL_CPPFLAGS='-I$(top_srcdir)'
  TLDEVEL_CPPFLAGS+="/$TLDEVELDIR"
  TLDEVEL_LDFLAGS='-L${top_srcdir}'
  TLDEVEL_LDFLAGS+="/$TLDEVELDIR"
  TLDEVEL_LIB='${top_srcdir}'
  TLDEVEL_LIB+="/$TLDEVELDIR/libtldevel.la"
  #TLDEVEL_LIB="../${srcdir}/$TLDEVELDIR/libtldevel.la"
  if test "x$enable_configure_tldevel" = "xyes"; then
    # We can't use a literal, because $HTSDIR is user-provided and variable
    #ac_configure_args="${ac_configure_args} "
    AC_CONFIG_SUBDIRS($TLDEVELDIR)
  fi
  ;;
system)
  AC_CHECK_HEADER([tldevel.h],
   [AC_CHECK_LIB(tldevel, tldevel_version, [ax_cv_tldevel=yes], [ax_cv_tldevel=no])],
  [ax_cv_tldevel=no], [;])
  ax_cv_tldevel_which=install
  TLDEVELDIR=
  #TLDEVEL_CPPFLAGS=
  #TLDEVEL_LDFLAGS=
  TLDEVEL_LIB=-ltldevel
  ;;
install)
  ax_saved_CPPFLAGS=$CPPFLAGS
  ax_saved_LDFLAGS=$LDFLAGS
  TLDEVEL_CPPFLAGS="-I$TLDEVELDIR/include"
  TLDEVEL_LDFLAGS="-L$TLDEVELDIR/lib"
  CPPFLAGS="$CPPFLAGS $TLDEVEL_CPPFLAGS"
  LDFLAGS="$LDFLAGS $TLDEVEL_LDFLAGS"
  AC_CHECK_HEADER([tldevel.h],
    [AC_CHECK_LIB(tldevel, tldevel_version, [ax_cv_tldevel=yes], [ax_cv_tldevel=no])],
    [ax_cv_tldevel=no], [;])
  TLDEVELDIR=
  CPPFLAGS=$ax_saved_CPPFLAGS
  LDFLAGS=$ax_saved_LDFLAGS
  ;;
none)
  ax_cv_tldevel=no
  ;;
esac
AC_SUBST([TLDEVELDIR])
AC_SUBST([TLDEVEL_CPPFLAGS])
AC_SUBST([TLDEVEL_LDFLAGS])
AC_SUBST([TLDEVEL_LIB])


if test "$ax_cv_tldevel" != yes; then
AC_MSG_ERROR([tldevelfiles not found
Samtools uses HTSlib to parse bioinformatics file formats etc.  Building it
requires an unpacked HTSlib source tree (which will be built in conjunction
with samtools) or a previously-installed HTSlib.  In either case you may
need to configure --with-htslib=DIR to locate the appropriate HTSlib.
FAILED.  You must supply an HTSlib in order to build samtools successfully.])
fi


])
