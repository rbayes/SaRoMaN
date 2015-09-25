dnl Process this file with autoconf to produce a configure script.
AC_PREREQ(2.53)
AC_INIT(lib/version.cpp)
AC_CANONICAL_TARGET([])

dnl Initi automake with the name of the package and the version
AM_INIT_AUTOMAKE(mind,"v0r1p0")
dnl AC_CONFIG_HEADERS(config.h)
AM_CONFIG_HEADER(config.h)

dnl things required by automake

AC_PROG_MAKE_SET([])


dnl Checks for programs.
AC_LANG([C++])
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX


AC_MSG_CHECKING(whether current system is supported)
case "$host_os" in
*darwin*)
  LDFLAGS="-flat_namespace"
  AC_DISABLE_SHARED	  
esac


AC_PROG_LIBTOOL

AC_MSG_CHECKING([for math library])
AC_CHECK_LIB(m,main)	

AC_MSG_CHECKING([for g2c library])
AC_CHECK_LIB(g2c,main)	


AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_CHECK_TOOL(AR, ar, :)
AC_HEADER_STDC


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

AC_SUBST(CXXFLAGS)

CXXFLAGS="-g -Wall $(recpack-config --include) $(clhep-config --include) $(bhep-config --include) -I$(root-config --incdir)" 

AC_CHECK_HEADERS([zlib.h])
AC_CHECK_HEADERS([particle.h])
AC_CHECK_HEADERS([Volume.h])

#linker flags 
AC_SUBST(LDFLAGS)
AC_SUBST(AM_LDFLAGS)

dnl Use standard <sstream> if found; 
AC_CHECK_HEADERS([sstream])

dnl Use standard <ostream> if found; 
AC_CHECK_HEADERS([ostream])


dnl change path if your libraries are somewhere else
AC_SUBST(LIBS)

prog_list="root-config recpack-config bhep-config clhep-config"
for prog in $prog_list; do
  prog=`which $prog 2> /dev/null`
  test ! -z $prog || AC_MSG_ERROR($prog not found in paths='$PATH': Check your installation)
  LDFLAGS="$LDFLAGS -L$($prog --prefix)/lib"
  LIBS="$LIBS $($prog --libs)"
  test $prog = "root-config" && LIBS="$LIBS $($prog --glibs)"
done

AC_MSG_CHECKING([for bhep library])
AC_CHECK_LIB(bhep, main)

AC_MSG_CHECKING([for CLHEP library])
AC_CHECK_LIB(CLHEP, main)

AC_MSG_CHECKING([for recpack library])
AC_CHECK_LIB(recpack, main)

AC_MSG_CHECKING([for root library])
AC_CHECK_LIB(root, main)



# LIBS="$LIBS $(root-config --libs)"

# Configuration inquiry program/script:
AC_CONFIG_FILES([mind-config], [chmod +x mind-config])

echo "host               = $host_os "
echo "CXX                = $CXX"
echo "CXXFLAGS           = $CXXFLAGS"
echo "CPPFLAGS           = $CPPFLAGS"
echo "LIBS               = $LIBS"
echo "LDFLAGS            = $LDFLAGS"
echo ""


AC_CONFIG_FILES([
	Makefile 
	mind/Makefile
	base/Makefile
	examples/Makefile
	lib/Makefile
        util/Makefile
	doc/Makefile
	doc/Doxyfile
	])
AC_OUTPUT
