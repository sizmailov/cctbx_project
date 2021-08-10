#include "simtbx/kokkos/kokkos_utils.h"

namespace simtbx { namespace Kokkos {
  void
  transfer_double2kokkos(vector_cudareal_t &dst, const double *src, const size_t length) {
    if (dst.span() < length) {
      resize(dst, length);
    }

    transfer_X2kokkos(dst, src, length);
  }

} // Kokkos
} // simtbx
