#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"
libGint_ver="EXP"
#spglib_sha256="c295dbea7d2fc9e50639aa14331fef277878c35f00ef0766e688bfbb7b17d44c"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libGint" ] && rm "${BUILDDIR}/setup_libGint"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libGint" in
  __INSTALL__)
    echo "==================== Installing libGint ===================="
    pkg_install_dir="${INSTALLDIR}/libGint-${libGint_ver}"
    install_lock_file="$pkg_install_dir/install_successful"

#    if verify_checksums "${install_lock_file}"; then
#      echo "libGint-${libGint_ver} is already installed, skipping it."
#    else
#      if [ -f libGint-${libGint_ver}.tar.gz ]; then
#        echo "libGint-${libGint_ver}.tar.gz is found"
#      else
#        download_pkg_from_cp2k_org "${spglib_sha256}" "spglib-${spglib_ver}.tar.gz"
#      fi
#      echo "Installing from scratch into ${pkg_install_dir}"
#      rm -rf spglib-${spglib_ver} "${pkg_install_dir}"
#      tar -xzf spglib-${spglib_ver}.tar.gz

#     OR

      rm -rf libGint-${libGint_ver}
      cp -r /work4/scd/scarf1152/libGint libGint-${libGint_ver}
#

      cd libGint-${libGint_ver}

#      mkdir build
#      cd build
#      cmake \
#        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
#        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
#        -DCMAKE_VERBOSE_MAKEFILE=ON \
#        -DSPGLIB_SHARED_LIBS=OFF \
#        -DSPGLIB_USE_OMP=ON \
#        -DSPGLIB_WITH_TESTS=OFF \
#        .. > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      echo "PKG INSTALL DIR IS " $pkg_install_dir
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install PREFIX="${pkg_install_dir}" > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
#    fi

    LIBGINT_CFLAGS="-I${pkg_install_dir}/include"
    LIBGINT_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    if [ -d "${pkg_install_dir}/lib64" ]; then
      ln -sf lib64 ${pkg_install_dir}/lib
      cd ${pkg_install_dir}
    fi
    ;;

  __SYSTEM__)
    echo "==================== Finding libGint from system paths ===================="
    echo "choose install next time. Or path. "
    exit 1
#    check_command pkg-config --modversion libGint
#    add_include_from_paths LIBGINT_CFLAGS "include/libGint.mod" $INCLUDE_PATHS
#    add_lib_from_paths SPGLIB_LDFLAGS "libspglib.*" $LIB_PATHS
    ;;

  __DONTUSE__) ;;
# pass

  *)
    echo "==================== Linking libGint to user paths ===================="
    echo "$with_libGint"
    echo "$with_libgint"
    pkg_install_dir="$with_libGint"
    check_dir "$pkg_install_dir/lib"
    check_dir "$pkg_install_dir/include"
    LIBGINT_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBGINT_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_libGint" != "__DONTUSE__" ]; then
  LIBGINT_LIBS="-lcp2kGint"
  if [ "$with_libGint" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libGint"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libGint"
export LIBGINT_VER="${libGint_ver}"
export LIBGINT_CFLAGS="IF_CUDA(-I${pkg_install_dir}/include ${LIBGINT_CFLAGS}|)"
export LIBGINT_LDFLAGS="IF_CUDA(${LIBGINT_LDFLAGS}|)"
export CP_DFLAGS="\${CP_DFLAGS} IF_CUDA(-D__LIBGINT|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_CUDA(${LIBGINT_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_CUDA(${LIBGINT_LDFLAGS}|)"
export CP_LIBS=" IF_CUDA(${LIBGINT_LIBS}|) \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_libGint" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_libGint"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libGint"
