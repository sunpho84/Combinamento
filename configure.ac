AC_INIT([Combinatoria],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_CC

AC_CHECK_LIB([m],[cos])

#silent automake
AM_SILENT_RULES([yes])

AC_CONFIG_FILES(Makefile)

AC_OUTPUT
